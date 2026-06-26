import numpy as np
import subprocess
import sys
import tomlkit
from pathlib import Path

script_dir  = Path(__file__).parent
output_toml = script_dir / 'pySB99_checkpoint_orig.toml'
tmp_dir     = Path('/tmp')
prefix      = 'orig'
file_prefix = f'pySB99_{prefix}_'

print('Running original...')
subprocess.run(
    [sys.executable, 'pySB99_test.py'],
    check=True,
    cwd=script_dir,
)

# ── Auto-discover checkpoint files ─────────────────────────────────────────────
checkpoint_files = sorted(tmp_dir.glob(f'{file_prefix}*.npy'))
if not checkpoint_files:
    print(f'No checkpoint files found in {tmp_dir} with prefix "{file_prefix}"')
    sys.exit(1)
print(f'Found {len(checkpoint_files)} checkpoint files')

# ── Load and truncate ──────────────────────────────────────────────────────────
def load_checkpoint(path, max_values=5):
    arr = np.load(path, allow_pickle=True).flatten()
    if len(arr) > max_values:
        return arr[:max_values].tolist()
    return arr.tolist()

# ── Build and write TOML ───────────────────────────────────────────────────────
doc = tomlkit.document()
doc.add(tomlkit.comment('pySB99 original script checkpoint values'))
doc.add(tomlkit.comment('Arrays longer than 5 elements are truncated to first 5 values'))
doc.add(tomlkit.nl())

for fpath in checkpoint_files:
    name = fpath.stem.replace(file_prefix, '')
    values = load_checkpoint(fpath)
    doc.add(name, values)

output_toml.write_text(tomlkit.dumps(doc))
print(f'Checkpoints saved to {output_toml}')