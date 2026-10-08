"""pandas >= 3 compatibility shim for TIR-Learner 3 inside the EDTA image.

TIR-Learner's check_TIR_TSD.py reads `family = x[0]` from a DataFrame row
whose labels are column names. pandas 2.x fell back to the position there
(with a FutureWarning); pandas 3 raises KeyError: 0, and EDTA_raw.pl --type
tir then finds no TIR-Learner output. TIR-Learner pins pandas 2.2.2, the
EDTA 2.3.0 image ships pandas 3.

dna_te_windows.py puts this directory first on PYTHONPATH for EDTA_raw.pl
only, so Python imports this file at startup (sitecustomize) in every
process EDTA starts. It restores the pandas 2 fallback for exactly that
case -- an int key, a Series whose index is not integer-typed, and the label
not found -- and changes nothing else. A no-op with pandas < 3 or without
pandas."""

try:
    import pandas as _pd
    _major = int(_pd.__version__.split(".")[0])
except Exception:  # no pandas, or an odd version string: leave everything alone
    _pd, _major = None, 0

if _pd is not None and _major >= 3 and not getattr(_pd.Series.__getitem__, "_positional_shim", False):
    _orig_getitem = _pd.Series.__getitem__

    def _getitem(self, key):
        try:
            return _orig_getitem(self, key)
        except KeyError:
            if (type(key) is int and not _pd.api.types.is_integer_dtype(self.index.dtype)
                    and -len(self) <= key < len(self)):
                return self.iloc[key]
            raise

    _getitem._positional_shim = True
    _pd.Series.__getitem__ = _getitem
