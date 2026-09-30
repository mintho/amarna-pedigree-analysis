"""Inventory the portable calculation bundle, excluding local publication assets."""

from pathlib import Path

EXCLUDED_PARTS = {'.git', 'reproduced', '__pycache__', '.DS_Store', 'build', '.venv'}
LOCAL_ONLY = (
    'Data/data/source_cache', 'Figures', 'Supplementary_Source',
    'Supplementary_Material.pdf', 'Results/Manuscript_Tables/en',
    'Results/Manuscript_Tables/de',
)


def retained_files(root):
    root = Path(root)
    paths = []
    for path in root.rglob('*'):
        relative = path.relative_to(root)
        if not path.is_file() or any(part in EXCLUDED_PARTS for part in relative.parts):
            continue
        if any(relative == Path(prefix) or Path(prefix) in relative.parents for prefix in LOCAL_ONLY):
            continue
        paths.append(path)
    return sorted(paths)
