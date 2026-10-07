#!/usr/bin/env python3
"""
Write Mooring/VTK/mdv2.vtk.series, so ParaView reads the MoorDyn line files
(mdv2_NNNN.vtk) at their simulation times instead of as steps 0, 1, 2, ...

Each file's time is read from its header ("MoorDyn v2 vtk output time=...").
Can be run at any time, also while the case is still running. Then in
ParaView open case.foam and Mooring/VTK/mdv2.vtk.series together; they share
the time axis.
"""

import glob
import json
import os
import re

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
VTK_DIR = os.path.join(CASE_DIR, "Mooring", "VTK")
PREFIX = "mdv2"

files = []
for path in sorted(glob.glob(os.path.join(VTK_DIR, f"{PREFIX}_*.vtk"))):
    with open(path) as handle:
        handle.readline()
        match = re.search(r"time=\s*([-+0-9.eE]+)", handle.readline())
    if match:
        files.append({"name": os.path.basename(path), "time": float(match.group(1))})
    else:
        print(f"No time in the header of {os.path.basename(path)}: skipped")

if not files:
    raise SystemExit(f"No {PREFIX}_*.vtk files with a time in {VTK_DIR}")

series = os.path.join(VTK_DIR, f"{PREFIX}.vtk.series")
with open(series, "w") as handle:
    json.dump({"file-series-version": "1.0", "files": files}, handle, indent=1)

print(
    f"Wrote {os.path.relpath(series, CASE_DIR)}: {len(files)} files, "
    f"t = {files[0]['time']} to {files[-1]['time']} s"
)
