"""
Pytest configuration and fixtures for mwprop tests
"""
import pytest
import os
import sys
import shutil

# Force non-interactive matplotlib backend for tests
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Add src to path so we can import mwprop modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


@pytest.fixture(scope="session", autouse=True)
def setup_output_directory():
    """
    Create output directory needed for dmdsm tests
    Uses os.getcwd() to match what dmdsm.py expects
    """
    output_dir = os.path.join(os.getcwd(), 'output_ne2025p')
    os.makedirs(output_dir, exist_ok=True)
    yield
    # Cleanup after tests
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir, ignore_errors=True)


@pytest.fixture(autouse=True)
def mpl_rcparams():
    """
    Apply stable matplotlib defaults for image-based tests.
    """
    original = plt.rcParams.copy()
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 100,
        "font.size": 10,
    })
    yield
    plt.rcParams.update(original)
