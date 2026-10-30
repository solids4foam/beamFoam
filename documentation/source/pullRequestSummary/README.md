# Source of documentation/pullRequestSummary.pdf

- `extract.py` reads both git histories (beamFoam since `origin/openfoam-v2306`, moorFV
  since `origin/main`) and writes `data.json`: commits, files, +/- counts and hunk ranges
- `explain_*.json`: a hand-checked description of each commit (keyed by short SHA)
- `gen.py` builds `prSummary.html` from `data.json`, the explanations, `head.html` and
  `tail.html` (open issues and the Phase 5 plan)

    python3 extract.py && python3 gen.py
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new \
        --disable-gpu --no-pdf-header-footer \
        --print-to-pdf=../../pullRequestSummary.pdf "file://$PWD/prSummary.html"

New commits need an entry in one of the `explain_*.json` files; `gen.py` prints any
commit without one.
