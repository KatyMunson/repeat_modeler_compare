"""Compatibility shims for TIR-Learner 3 inside the EDTA 2.3.0 image.

1. pandas >= 3 positional fallback (below).
2. swifter: no dask on string frames, thread scheduler (after it).
3. apply(axis=1) on an empty DataFrame (last).

--- 1. pandas >= 3 ---

TIR-Learner's check_TIR_TSD.py reads `family = x[0]` from a DataFrame row
whose labels are column names. pandas 2.x fell back to the position there
(with a FutureWarning); pandas 3 raises KeyError: 0, and EDTA_raw.pl --type
tir then finds no TIR-Learner output. TIR-Learner pins pandas 2.2.2, the
EDTA 2.3.0 image ships pandas 3.

This file restores the pandas 2 fallback for exactly that case -- an int
key, a Series whose index is not integer-typed, and the label not found --
and changes nothing else. A no-op with pandas < 3 or without pandas.

--- 2. swifter ---
get_fasta_sequence.py runs df.swifter.apply inside a multiprocessing.Pool
worker. swifter (1.4.0, pinned by TIR-Learner) parallelises an apply it
estimates to be slow with dask's *process* scheduler, and a pool worker is
daemonic and may not start processes: "AssertionError: daemonic processes
are not allowed to have children". It only happens in windows with many
TIRvish hits -- the TIR-richest windows, which must not be skipped. When
swifter is imported, its objects get the thread scheduler instead; the
apply is a string slice, so nothing is lost.

Why dask ran at all: swifter keeps an apply on a frame with strings in
pandas (no dask) by testing for dtype "object". TIR-Learner's frames hold
sequences and ids; under its pinned pandas 2.2.2 those were object, so it
never used dask. pandas 3 stores them as "str", the test misses them, and
any apply swifter estimates at > 1 s (large candidate tables, a busy node)
went to dask -- a path TIR-Learner was never run on, and dask threads
followed by multiprocessing's fork can deadlock. So a swifter object on a
frame with any string column gets an infinite dask threshold: the plain
pandas apply TIR-Learner was written against.

--- 3. empty frames ---
DataFrame.apply(f, axis=1) on zero rows returns a copy of the (empty)
frame, not a Series, when f fails on an empty row -- every TIR-Learner
lambda does. `df["check_N_TIR1_check"] = df.apply(...)` then raises
"Cannot set a DataFrame with multiple columns to the single column" (same
in pandas 2). It happens when TIRvish finds candidates and the N / TA
filters drop all of them: a window with nothing left to report fails
instead of reporting nothing. Here such an apply returns an empty Series
(the column assignment works, the frame stays empty, TIR-Learner's own
"no candidates" branch follows). Only zero-row frames, axis=1, no
result_type: nothing that has data changes.

dna_te_windows.py puts this directory first on PYTHONPATH for EDTA_raw.pl
only, so Python imports this file at startup (sitecustomize) in every
process EDTA starts."""

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


# --- 2. swifter: thread scheduler, patched when swifter.swifter is imported
import sys as _sys


def _has_strings(obj):
    """True for a pandas object with a string column (object or str dtype)."""
    try:
        from pandas.api.types import is_string_dtype
        dtypes = obj.dtypes if hasattr(obj, "columns") else [obj.dtype]
        return any(is_string_dtype(t) for t in dtypes)
    except Exception:
        return False


def _patch_swifter(module):
    cls = getattr(module, "_SwifterObject", None)
    if cls is None or getattr(cls.__init__, "_thread_shim", False):
        return
    _orig_init = cls.__init__

    def __init__(self, *args, **kwargs):
        _orig_init(self, *args, **kwargs)
        if getattr(self, "_scheduler", None) == "processes":
            self._scheduler = "threads"
        if _has_strings(getattr(self, "_obj", None)):
            self._dask_threshold = float("inf")

    __init__._thread_shim = True
    cls.__init__ = __init__


class _SwifterFinder:
    """Meta-path hook: load swifter.swifter normally, then patch it."""

    def find_spec(self, name, path=None, target=None):
        if name != "swifter.swifter":
            return None
        for finder in _sys.meta_path:
            if finder is self or not hasattr(finder, "find_spec"):
                continue
            spec = finder.find_spec(name, path, target)
            if spec is not None and spec.loader is not None and hasattr(spec.loader, "exec_module"):
                orig_exec = spec.loader.exec_module

                def exec_module(module, _orig=orig_exec):
                    _orig(module)
                    _patch_swifter(module)

                spec.loader.exec_module = exec_module
                return spec
        return None


if "swifter.swifter" in _sys.modules:
    _patch_swifter(_sys.modules["swifter.swifter"])
else:
    _sys.meta_path.insert(0, _SwifterFinder())


# --- 3. apply(axis=1) on an empty DataFrame: an empty Series, not a frame copy
if _pd is not None and not getattr(_pd.DataFrame.apply, "_empty_shim", False):
    _orig_apply = _pd.DataFrame.apply

    def _apply(self, func, axis=0, *args, **kwargs):
        res = _orig_apply(self, func, axis, *args, **kwargs)
        result_type = kwargs.get("result_type", args[1] if len(args) > 1 else None)
        if (len(self) == 0 and axis in (1, "columns") and result_type is None
                and isinstance(res, _pd.DataFrame)):
            return _pd.Series(index=self.index, dtype=object)
        return res

    _apply._empty_shim = True
    _pd.DataFrame.apply = _apply
