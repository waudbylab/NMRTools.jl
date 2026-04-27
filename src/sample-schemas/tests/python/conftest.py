import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MIGRATION_CODE_DIR = os.path.join(REPO_ROOT, "migration-code")
FIXTURES_DIR = os.path.join(REPO_ROOT, "tests", "fixtures")
PATCH_PATH = os.path.join(REPO_ROOT, "current", "patch.json")

if MIGRATION_CODE_DIR not in sys.path:
    sys.path.insert(0, MIGRATION_CODE_DIR)
