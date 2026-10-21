"""Regenerate every figure: runs each script in figures/src/ (or only those named on the command line).

    python figures/build_all.py            # all
    python figures/build_all.py ch26 ch27  # scripts whose name starts with ch26 or ch27
"""

import runpy
import sys
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

prefixes = sys.argv[1:]
scripts = sorted((HERE / "src").glob("*.py"))
if prefixes:
    scripts = [s for s in scripts if s.stem.startswith(tuple(prefixes))]

failed = []
for script in scripts:
    try:
        runpy.run_path(str(script), run_name="__main__")
    except Exception:
        traceback.print_exc()
        failed.append(script.name)

print(f"{len(scripts) - len(failed)}/{len(scripts)} figure scripts succeeded")
if failed:
    print("failed:", ", ".join(failed))
    sys.exit(1)
