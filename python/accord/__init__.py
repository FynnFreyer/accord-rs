import sys

from ._internal import Calculator, data

# necessary for imports, see discussion at: https://github.com/PyO3/pyo3/issues/759
sys.modules["accord.data"] = data
sys.modules["accord.data.indel"] = data.indel
sys.modules["accord.data.stats"] = data.stats

__all__ = ["Calculator", "data"]
