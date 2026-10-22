"""Generate the PR summary HTML from data.json and explain_*.json"""
import collections
import glob
import html
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
# data.json is written by extract.py
data = json.load(open(os.path.join(HERE, "data.json")))
explain = {}
for path in sorted(glob.glob(os.path.join(HERE, "explain_*.json"))):
    explain.update(json.load(open(path)))

E = html.escape
GH = {"beamFoam": "https://github.com/solids4foam/beamFoam/commit/",
      "moorFV": "https://github.com/solids4foam/moorFV/commit/"}

SCRIPT_NAMES = {"Allrun", "Allclean", "Alltest", "makeCases", "Allwmake"}
SCRIPT_EXT = (".py", ".sh", ".gp", ".plt", ".gnuplot")
DOC_PREFIX = ("codexLogs/", "documentation/")


def new_path(p):
    m = re.match(r"(.*)\{(.*) => (.*)\}(.*)", p)
    if m:
        return re.sub("//+", "/", m.group(1) + m.group(3) + m.group(4))
    return p.split(" => ")[-1]


def kind(p):
    q = new_path(p)
    name = q.split("/")[-1]
    if q.startswith(DOC_PREFIX) or name.lower().endswith((".md", ".pdf", ".txt", ".html", ".png")):
        return "doc"
    if q.startswith("tutorials/"):
        if name in SCRIPT_NAMES or name.endswith(SCRIPT_EXT):
            return "script"
        return "case"
    return "code"


CASE_SUBDIRS = {"0", "0.orig", "constant", "system", "postProcessing", "plots", "references", "referenceResults"}


def case_dir(p):
    parts = new_path(p).split("/")
    out = []
    for part in parts[:-1]:
        if part in CASE_SUBDIRS or re.fullmatch(r"[0-9.e+-]+", part):
            break
        out.append(part)
    return "/".join(out) or "tutorials"


def ranges(f, limit=10):
    if f["add"] == "-":
        return "binary"
    out = []
    for os_, oc, ns, nc in f["hunks"]:
        if nc == 0:
            out.append(f"−{oc} at {ns}")
        elif nc == 1:
            out.append(f"{ns}")
        else:
            out.append(f"{ns}–{ns + nc - 1}")
    if not out:
        return "rename only" if "=>" in f["path"] else ""
    extra = len(out) - limit
    text = ", ".join(out[:limit])
    return text + (f" (+{extra} more)" if extra > 0 else "")


def num(x):
    return 0 if x == "-" else int(x)


def pm(a, d):
    if a == "-":
        return '<span class="bin">bin</span>'
    return f'<span class="add">+{a}</span> <span class="del">−{d}</span>'


def short(p):
    q = E(p)
    for pre in ("src/wireBunchingModels/", "src/sDoFRGBFvBeam/"):
        q = q.replace(pre, f'<span class="pre">{pre}</span>', 1)
    return q


def file_table(files, max_code=30):
    code = [f for f in files if kind(f["path"]) in ("code", "script")]
    docs = [f for f in files if kind(f["path"]) == "doc"]
    cases = collections.OrderedDict()
    for f in files:
        if kind(f["path"]) == "case":
            c = cases.setdefault(case_dir(f["path"]), [0, 0, 0])
            c[0] += 1; c[1] += num(f["add"]); c[2] += num(f["del"])
    rows = []
    shown = code[:max_code]
    for f in shown:
        rows.append(f'<tr><td class="path">{short(f["path"])}</td><td class="pm">{pm(f["add"], f["del"])}</td>'
                    f'<td class="lines">{E(ranges(f))}</td></tr>')
    if len(code) > max_code:
        rest = code[max_code:]
        dirs = collections.Counter("/".join(new_path(f["path"]).split("/")[:4]) for f in rest)
        a = sum(num(f["add"]) for f in rest); d = sum(num(f["del"]) for f in rest)
        rows.append(f'<tr><td class="path">… {len(rest)} more code files in '
                    + ", ".join(E(k) + f" ({v})" for k, v in dirs.most_common(6))
                    + f'</td><td class="pm">{pm(str(a), str(d))}</td><td class="lines"></td></tr>')
    for f in docs:
        rows.append(f'<tr class="docrow"><td class="path">{short(f["path"])}</td><td class="pm">{pm(f["add"], f["del"])}</td>'
                    f'<td class="lines">document</td></tr>')
    for k, (n, a, d) in cases.items():
        rows.append(f'<tr class="caserow"><td class="path">{E(k)}/ <span class="muted">(case files)</span></td>'
                    f'<td class="pm">{pm(str(a), str(d))}</td><td class="lines">{n} file{"s" if n != 1 else ""}</td></tr>')
    if not rows:
        return '<p class="muted small">No files changed (empty commit).</p>'
    return ('<table class="files"><tr><th>File</th><th class="pmh">Lines</th><th>Changed line ranges (after the commit)</th></tr>'
            + "".join(rows) + "</table>")


def card(repo, c):
    ex = explain.get(c["sha"], {})
    summary = ex.get("summary", "")
    details = ex.get("details", [])
    tag = f'<span class="repo {repo}">{repo}</span>'
    body = f'<p class="summary">{E(summary)}</p>' if summary else ""
    if details:
        body += "<ul class='det'>" + "".join(f"<li>{E(x)}</li>" for x in details) + "</ul>"
    return (f'<div class="commit"><div class="ch">{tag}<a class="sha" href="{GH[repo]}{c["sha"]}">{c["sha"]}</a>'
            f'<span class="date">{c["date"]}</span><span class="who">{E(c["author"])}</span></div>'
            f'<div class="subj">{E(c["subject"])}</div>{body}{file_table(c["files"])}</div>')


# ---------------------------------------------------------------- eras
commits = [("beamFoam", c) for c in data["beamFoam"]["commits"]] + [("moorFV", c) for c in data["moorFV"]["commits"]]
commits.sort(key=lambda rc: rc[1]["time"])

ERAS = [
    ("era1", "2.1 June – September 2025: beamFoam refactoring, function objects, code removal and tutorials",
     "Work by Seevani Bali, with Philip Cardiff and Amirhossein Taran, on the <code>openfoam-v2306-seevani</code> line. It is "
     "merged into this branch, but not yet into <code>openfoam-v2306</code>, so it shows up in the PR diff.",
     lambda c: c["date"] < "2025-11-01"),
    ("era2", "2.2 November – December 2025: momentum contributions, spring boundary condition, build fixes",
     "Also part of the <code>openfoam-v2306-seevani</code> work. The first moorFV commit (local changes) falls in this period.",
     lambda c: "2025-11-01" <= c["date"] < "2026-01-01"),
    ("era3", "2.3 February – August 2026: first rigid-body coupling attempt (Colm Mc Alister, with Codex)",
     "Adds rigid-body unknowns to the BlockEigen system and moorFV state passing. The review on 30 September 2026 "
     "(<code>rigidBodyCouplingPRSummary.pdf</code>) found that this path is still a partitioned scheme. The descriptions below "
     "say what the code does.",
     lambda c: "2026-01-01" <= c["date"] < "2026-09-30" or c["sha"] in {"778e43d"}),
    ("era4", "2.4 30 September – 1 October 2026: review and monolithic coupling, Phases 0–4 (Colm Mc Alister, with Claude)",
     "Review of the earlier path, the development plan, then beamFoam solving the body monolithically "
     "(<code>rigidBodyEnd</code>), moorFV's <code>beamFoamCoupled</code> solver, and validation in interFoam.",
     lambda c: c["date"] >= "2026-09-30" and c["sha"] not in {"778e43d"}),
]


def era_of(c):
    for key, _, _, test in ERAS:
        if test(c):
            return key


by_era = collections.defaultdict(list)
for repo, c in commits:
    by_era[era_of(c)].append((repo, c))

# ---------------------------------------------------------------- net tables
def net_rows(repo):
    rows = []
    dele = collections.defaultdict(lambda: [0, 0])
    for a, d, p in data[repo]["net"]:
        k = kind(p)
        if k == "case":
            continue
        root = {"beamFoam": "/Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV/src/beamFoam",
                "moorFV": "/Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV"}[repo]
        if not os.path.lexists(os.path.join(root, new_path(p))):
            parts = new_path(p).split("/")
            key = "/".join(parts[:min(len(parts) - 1, 4 if parts[0] == "src" else 3)])
            dele[key][0] += 1; dele[key][1] += num(d)
            continue
        rows.append((p, a, d, k))
    return rows, dele


def net_table(repo):
    rows, dele = net_rows(repo)
    groups = collections.OrderedDict()
    for p, a, d, k in rows:
        q = new_path(p)
        g = "/".join(q.split("/")[:3]) if q.startswith("src/") else q.split("/")[0] if "/" in q else "(top level)"
        groups.setdefault(g, []).append((p, a, d, k))
    out = ['<table class="files"><tr><th>File</th><th class="pmh">Lines</th><th>Type</th></tr>']
    for g, items in groups.items():
        out.append(f'<tr class="grp"><td colspan="3">{E(g)}</td></tr>')
        for p, a, d, k in items:
            out.append(f'<tr><td class="path">{short(p)}</td><td class="pm">{pm(a, d)}</td><td class="lines">{k}</td></tr>')
    if dele:
        out.append('<tr class="grp"><td colspan="3">Deleted files, grouped by directory</td></tr>')
        for g, (n, d) in sorted(dele.items()):
            out.append(f'<tr><td class="path">{E(g)}/</td><td class="pm"><span class="del">−{d}</span></td>'
                       f'<td class="lines">{n} files</td></tr>')
    out.append("</table>")
    return "".join(out)


def case_summary(repo):
    cases = collections.OrderedDict()
    for a, d, p in data[repo]["net"]:
        if kind(p) == "case":
            c = cases.setdefault(case_dir(p), [0, 0, 0])
            c[0] += 1; c[1] += num(a); c[2] += num(d)
    out = ['<table class="files"><tr><th>Tutorial case</th><th class="pmh">Lines</th><th>Files</th></tr>']
    for k, (n, a, d) in sorted(cases.items()):
        out.append(f'<tr><td class="path">{E(k)}/</td><td class="pm">{pm(str(a), str(d))}</td><td class="lines">{n}</td></tr>')
    out.append("</table>")
    return "".join(out)


def totals(repo):
    a = sum(num(x[0]) for x in data[repo]["net"]); d = sum(num(x[1]) for x in data[repo]["net"])
    code = [x for x in data[repo]["net"] if kind(x[2]) in ("code", "script")]
    ca = sum(num(x[0]) for x in code); cd = sum(num(x[1]) for x in code)
    return len(data[repo]["net"]), a, d, len(code), ca, cd


# ---------------------------------------------------------------- page
parts = open(os.path.join(HERE, "head.html")).read()
bt, mt = totals("beamFoam"), totals("moorFV")
nb, nm = len(data["beamFoam"]["commits"]), len(data["moorFV"]["commits"])
parts = parts.replace("{{BF_TOTALS}}", f"{bt[0]} files, +{bt[1]:,} / −{bt[2]:,} lines; code and scripts {bt[3]} files, +{bt[4]:,} / −{bt[5]:,}")
parts = parts.replace("{{MF_TOTALS}}", f"{mt[0]} files, +{mt[1]:,} / −{mt[2]:,} lines; code {mt[3]} files, +{mt[4]:,} / −{mt[5]:,}")
parts = parts.replace("{{BF_BASE}}", E(data["beamFoam"]["mergeBase"])).replace("{{MF_BASE}}", E(data["moorFV"]["mergeBase"]))
parts = parts.replace("{{NB}}", str(nb)).replace("{{NM}}", str(nm))

timeline = ['<table><tr><th>Period</th><th class="n">beamFoam commits</th><th class="n">moorFV commits</th><th>Authors</th></tr>']
for key, title, intro, _ in ERAS:
    items = by_era[key]
    authors = sorted({c["author"] for _, c in items})
    timeline.append(f'<tr><td><a href="#{key}">{E(title.split(":")[0])}: {E(title.split(": ", 1)[1])}</a></td>'
                    f'<td class="n">{sum(r == "beamFoam" for r, _ in items)}</td><td class="n">{sum(r == "moorFV" for r, _ in items)}</td>'
                    f'<td>{E(", ".join(authors))}</td></tr>')
timeline.append("</table>")
parts = parts.replace("{{TIMELINE}}", "".join(timeline))
parts = parts.replace("{{BF_NET}}", net_table("beamFoam")).replace("{{BF_CASES}}", case_summary("beamFoam"))
parts = parts.replace("{{MF_NET}}", net_table("moorFV"))

log = ['<h2 style="border:none;margin-bottom:0">2. Chronological commit log</h2>']
for key, title, intro, _ in ERAS:
    log.append(f'<h2 id="{key}">{E(title)}</h2><p>{intro}</p>')
    log.extend(card(r, c) for r, c in by_era[key])
parts = parts.replace("{{LOG}}", "".join(log))

merges = ['<table><tr><th style="width:10%">Merge</th><th style="width:12%">Date</th><th style="width:16%">Author</th><th>Subject</th></tr>']
for m in data["beamFoam"]["merges"]:
    sha, date, who, subj = m.split("|", 3)
    merges.append(f"<tr><td><code>{sha}</code></td><td>{date}</td><td>{E(who)}</td><td>{E(subj)}</td></tr>")
merges.append("</table>")
parts = parts.replace("{{MERGES}}", "".join(merges))
parts = parts.replace("{{TAIL}}", open(os.path.join(HERE, "tail.html")).read())

missing = [c["sha"] for _, c in commits if c["sha"] not in explain]
open(os.path.join(HERE, "prSummary.html"), "w").write(parts)
print("written; commits without explanation:", missing)
