import sys
sys.dont_write_bytecode = True
from pathlib import Path
import runpy
import unittest


def main():
    root = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(
        str(root/"tests"), pattern="test_*.py", top_level_dir=str(root))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        return 1
    # Joint exact supercell audit of freshly published, re-opened and
    # off-database evidence.  `supercell_cases.py --verify` is the same
    # workflow behind a script entry point.
    namespace = runpy.run_path(str(root/"supercell_cases.py"),
                               run_name="supercell_cases")
    namespace["workflow"]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
