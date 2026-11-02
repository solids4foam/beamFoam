import json, os, re, subprocess, sys
def git(repo, *a):
    return subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True, errors="replace").stdout

def extract(repo, base, name):
    mb = git(repo, "merge-base", "HEAD", base).strip()
    shas = git(repo, "log", "--reverse", "--no-merges", "--format=%H", f"{mb}..HEAD").split()
    merges = git(repo, "log", "--reverse", "--merges", "--format=%h|%ad|%an|%s", "--date=short", f"{mb}..HEAD").splitlines()
    out = []
    for sha in shas:
        meta = git(repo, "show", "-s", "--format=%h%x00%ad%x00%an%x00%s%x00%at%x00%b", "--date=short", sha).split("\x00")
        numstat = git(repo, "show", "--format=", "--numstat", "-M", sha)
        files = []
        for line in numstat.splitlines():
            a, d, path = line.split("\t", 2)
            files.append({"path": path, "add": a, "del": d, "hunks": []})
        diff = git(repo, "show", "--format=", "-U0", "-M", "--no-color", sha)
        cur = None
        for line in diff.splitlines():
            if line.startswith("diff --git"):
                m = re.match(r"diff --git a/(.*) b/(.*)", line); cur = m.group(2)
            m = re.match(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@", line)
            if m and cur:
                os_, oc, ns, nc = int(m.group(1)), int(m.group(2) or 1), int(m.group(3)), int(m.group(4) or 1)
                for f in files:
                    p = f["path"]
                    if p == cur or ("=>" in p and p.endswith(cur.split("/")[-1])):
                        f["hunks"].append([os_, oc, ns, nc]); break
        out.append({"sha": meta[0], "date": meta[1], "author": meta[2], "subject": meta[3], "time": int(meta[4]), "body": meta[5].strip(), "files": files})
    net = git(repo, "diff", "--numstat", "-M", mb, "HEAD")
    return {"repo": name, "mergeBase": git(repo, "show", "-s", "--format=%h %ad %s", "--date=short", mb).strip(),
            "commits": out, "merges": merges, "net": [l.split("\t", 2) for l in net.splitlines()]}

data = {"beamFoam": extract("/Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV/src/beamFoam", "origin/openfoam-v2306", "beamFoam"),
        "moorFV": extract("/Volumes/OpenFOAM/colmmcalister-v2312/run/moorFV", "origin/main", "moorFV")}
json.dump(data, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data.json"), "w"))
for k, v in data.items():
    nfiles = sum(len(c["files"]) for c in v["commits"])
    print(k, v["mergeBase"], len(v["commits"]), "non-merge commits,", len(v["merges"]), "merges,", nfiles, "file entries,", len(v["net"]), "net files")
