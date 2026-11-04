"""Run every gallery calculation and CSV postprocessing before its first renderer."""
import re, subprocess, tempfile
from pathlib import Path
from invocation import compiler_cmd
root = Path.cwd()
sources = sorted((root / 'examples').glob('*.lms'))
assert len(sources) == 17, len(sources)
with tempfile.TemporaryDirectory(prefix='example-data-', dir=root / 'build') as directory:
    for source in sources:
        text = source.read_text()
        assert 'r.stdout' not in text, source
        renderer = re.search(r'^\s*\w+\s*=\s*_lmcode_ cli="(?:lmplot|lmcad|lmdraw)"', text, re.M)
        assert renderer, source
        prefix = text[:renderer.start()].replace(' install="github:lmctl-ai/lmfluid@latest"', '')
        scratch = Path(directory) / source.stem
        scratch.mkdir()
        program = scratch / source.name
        program.write_text(prefix)
        result = subprocess.run(compiler_cmd() + [str(program)], cwd=scratch,
                                capture_output=True, text=True, timeout=600)
        assert result.returncode == 0, (source.name, result.stdout, result.stderr)
print(f'PASS {len(sources)} offline gallery calculations and CSV postprocessing')
