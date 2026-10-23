"""Rerun the local assessment in an isolated, temporary source copy."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess
import tempfile

audit = Path(__file__).resolve().parent
root = audit.parents[1]
scratch = Path(tempfile.mkdtemp(prefix="multippi-audit-"))
package = scratch / "multiPPI"
package.mkdir()
for name in ["R", "src", "tests", "man", "DESCRIPTION", "NAMESPACE", "LICENSE"]:
    source = root / name
    if source.is_dir():
        shutil.copytree(source, package / name,
                        ignore=shutil.ignore_patterns("*.o", "*.so", "*.dll"))
    elif source.exists():
        shutil.copy2(source, package / name)

manifest = json.loads((audit / "source-manifest.json").read_text())
changed = [name for name, digest in manifest.items()
           if not (root / name).exists()
           or hashlib.sha256((root / name).read_bytes()).hexdigest() != digest]
(scratch / "source-differences.json").write_text(json.dumps(changed, indent=2) + "\n")
print(f"Scratch and results: {scratch}", flush=True)
print(f"Files differing from the assessed manifest: {len(changed)}", flush=True)
env = dict(os.environ, LC_ALL="C", NOT_CRAN="true", AUDIT_PROBE="")
statuses = {}
for script in ["run-suite.R", "probes.R", "integration.R", "benchmark.R"]:
    with (scratch / f"{script}.log").open("w") as log:
        completed = subprocess.run(
            ["Rscript", "--vanilla", str(audit / script), str(package), str(scratch)],
            cwd=scratch, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    statuses[script] = completed.returncode
    print(f"{script}: exit {completed.returncode}", flush=True)
(scratch / "exit-status.json").write_text(json.dumps(statuses, indent=2) + "\n")
print("Inspect the captured errors and numerical results; probes capture expected failures.")
