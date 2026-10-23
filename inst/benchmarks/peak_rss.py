"""Run isolated PPI benchmark processes and retain macOS/Linux wait4 peak RSS."""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

root = Path(__file__).resolve().parents[2]
out = root / "audit/2026-10-23-ppi"
paths = [root / "R/ppi_design.R", root / "R/ppi_fit.R",
         root / "inst/benchmarks/benchmark_ppi.R"]
hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in paths}
receipts = []
for targets in (400, 1000, 50000):
    if os.statvfs(root).f_bavail * os.statvfs(root).f_frsize < 20 * 1024**3:
        raise RuntimeError("Less than 20 GiB free; benchmark paused")
    log = out / f"implementation-final-benchmark-{targets}.log"
    with log.open("w") as handle:
        start = time.monotonic()
        proc = subprocess.Popen(["Rscript", "inst/benchmarks/benchmark_ppi.R", str(targets)],
                                cwd=root, stdout=handle, stderr=subprocess.STDOUT)
        _, status, usage = os.wait4(proc.pid, 0)
        proc.returncode = os.waitstatus_to_exitcode(status)
    receipt = dict(targets=targets, elapsed_process_seconds=time.monotonic()-start,
                   exit_code=proc.returncode, peak_rss_bytes=usage.ru_maxrss *
                   (1 if platform.system() == "Darwin" else 1024),
                   log=str(log.relative_to(root)))
    receipts.append(receipt)
    print(json.dumps(receipt), flush=True)
    if proc.returncode:
        raise RuntimeError(f"Benchmark failed: {log}")
assert hashes == {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in paths}, "Source changed during benchmark"
(out / "implementation-final-benchmark-rss.json").write_text(json.dumps(
    dict(platform=platform.platform(), source_sha256=hashes,
         measurement="wait4 child ru_maxrss; whole process including input generation and three fits",
         results=receipts), indent=2) + "\n")
