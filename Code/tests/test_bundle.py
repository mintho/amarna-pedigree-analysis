from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools'))
from repository_files import retained_files
import reproduce


class BundleTests(unittest.TestCase):
    def test_inventory_excludes_local_pages_and_publication_assets(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            retained = ('Code/src/example.py', 'Data/data/example.json',
                        'Results/summary.csv', 'README.md', 'LICENSE')
            excluded = ('Data/data/source_cache/page.html', 'reproduced/source_cache/page.html',
                        'Figures/colour/figure.png', 'Supplementary_Material.pdf',
                        'Supplementary_Source/main.tex', 'Results/Manuscript_Tables/en/row.tex',
                        'Results/Manuscript_Tables/de/row.tex', '.git/config',
                        'Code/src/__pycache__/example.pyc')
            for relative in retained + excluded:
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch()
            self.assertEqual({str(p.relative_to(root)) for p in retained_files(root)}, set(retained))

    def test_source_retrieval_uses_ignored_output_cache(self):
        with patch.object(reproduce, 'run') as run:
            reproduce.priors()
        args = run.call_args_list[0].args
        cache = args[args.index('--cache_dir') + 1]
        self.assertEqual(cache, reproduce.OUT / 'source_cache/duesseldorf_archive')
        self.assertEqual(run.call_args_list[1].args[0], 'build_color_str_mapping.py')

    def test_retained_inputs_and_calculation_code_have_no_html(self):
        self.assertFalse(any(p.suffix.lower() in {'.html', '.htm', '.pdf', '.png'}
                             for p in retained_files(ROOT)))
