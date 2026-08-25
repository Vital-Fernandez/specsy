from pathlib import Path
import os
import subprocess
import sys
import time

Z_LIST = ['MWC', 'MW', 'LMC', 'SMC', 'IZw18', 'XMP', 'Z0']

script_file = Path('pySB99_cleanUppySB99_cleanUp.py')
env = {**os.environ, 'MPLBACKEND': 'Agg'}

results = {}
for Z in Z_LIST:

    print(f'--- {Z} ---', flush=True)
    st = time.time()
    run = subprocess.run([sys.executable, str(script_file), Z], capture_output=True, text=True, env=env)
    elapsed = time.time() - st

    results[Z] = run.returncode
    if run.returncode == 0:
        print(f'{Z}: OK ({elapsed:.0f} s)')
    else:
        print(f'{Z}: FAILED (exit {run.returncode})')
        print(run.stdout[-2000:])
        print(run.stderr[-2000:])

print('\n--- summary ---')
for Z, code in results.items():
    print(f'{Z:<8} {"OK" if code == 0 else f"FAILED ({code})"}')