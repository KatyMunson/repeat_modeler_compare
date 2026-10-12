"""Compatibility shims for TIR-Learner 3 inside the EDTA 2.3.0 image.

1. pandas >= 3 positional fallback (below).
2. swifter: no dask on string frames, thread scheduler (after it).
3. apply(axis=1) on an empty DataFrame.
4. TIR-Learner 3.0.x crashes when a step finds nothing (last).

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

--- 4. nothing found ---
The image runs TIR-Learner 3.0.2 (/usr/local/share/TIR-Learner3.0/bin),
which has no guard for an empty step result; 3.0.7 added most of them. A
window with nothing to find failed instead of reporting nothing:
- TIRvish finds no hits: get_fasta_pieces_SeqIO concatenates only Nones,
  "All objects passed were None" (before GRF has even run);
- GRF writes no candidate.fasta anywhere: the same, in run_GRF;
- every sequence is shorter than TIR-Learner's minimum (max TIR length +
  500 bp, 5.5 kb under EDTA): TIRvish and GRF get no files, "No objects to
  concatenate";
- neither TIRvish nor GRF leaves a candidate: concat of two Nones;
- any step that stores None (TIRvish or GRF with nothing left after its
  filters): the checkpoint writer calls None.to_csv;
- candidates found but none passes the TIR/TSD check: check_TIR_TSD returns
  an empty frame without its "source" column and combine_all's sort fails,
  "KeyError: 'source'".
Each is wrapped so the empty case does what the code meant: an empty frame
passes through, "no candidates at all" exits 0 with a NOTICE (EDTA_raw.pl
then writes an empty TIR file, as for a TIR-Learner run that finds
nothing), and combine_all ignores empty module results, and the checkpoint writer skips
None entries. Inputs with data
go to the original functions unchanged. Only modules loaded from a
TIR-Learner directory are touched.

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


def _patch_tirlearner(module):
    """Wrap TIR-Learner 3.0.x functions that crash on an empty step (4.)."""
    import functools
    pd = __import__("pandas")
    short = module.__name__.rsplit(".", 1)[-1]

    def wrap(name, make):
        f = getattr(module, name, None)
        if f is not None and not getattr(f, "_empty_shim", False):
            g = functools.wraps(f)(make(f))
            g._empty_shim = True
            setattr(module, name, g)

    def empty(df):
        return df is None or len(df) == 0

    def _nothing(e):
        # pd.concat of only Nones, or of an empty list
        return any(m in str(e) for m in ("All objects passed were None", "No objects to concatenate"))

    if short == "get_fasta_sequence":
        def make(f):
            def get_fasta_pieces_SeqIO(genome_file, df_in, *args, **kwargs):
                if df_in is not None and len(df_in) == 0:
                    out = df_in.copy()
                    out["seq"] = pd.Series(index=out.index, dtype=object)
                    return out
                return f(genome_file, df_in, *args, **kwargs)
            return get_fasta_pieces_SeqIO
        wrap("get_fasta_pieces_SeqIO", make)
    elif short == "run_GRF":
        def make(f):
            def get_GRF_result_df_para(*args, **kwargs):
                try:
                    return f(*args, **kwargs)
                except ValueError as e:
                    if not _nothing(e):
                        raise
                    print("NOTICE: GRF found no candidates.")
                    return pd.DataFrame({"id": pd.Series(dtype=object), "seq": pd.Series(dtype=object),
                                         "len": pd.Series(dtype=int)})
            return get_GRF_result_df_para
        wrap("get_GRF_result_df_para", make)
    elif short == "run_TIRvish":
        def make(f):
            def run_TIRvish_py_para(*args, **kwargs):
                try:
                    return f(*args, **kwargs)
                except ValueError as e:
                    if not _nothing(e):
                        raise
                    print("NOTICE: TIRvish had no sequence to search (all shorter than TIR-Learner's minimum).")
                    cols = ["seqid", "start", "end", "TIR1_start", "TIR1_end", "TIR2_start", "TIR2_end", "id"]
                    return pd.DataFrame({c: pd.Series(dtype=object if c in ("seqid", "id") else int)
                                         for c in cols})
            return run_TIRvish_py_para
        wrap("run_TIRvish_py_para", make)
    elif short == "process_de_novo_result":
        def make(f):
            def combine_de_novo_result(inst, *args, **kwargs):
                if empty(inst.get("TIRvish", None)) and empty(inst.get("GRF", None)):
                    print("NOTICE: No TIR candidates found by TIRvish or GRF; nothing to classify.", flush=True)
                    _sys.exit(0)
                return f(inst, *args, **kwargs)
            return combine_de_novo_result
        wrap("combine_de_novo_result", make)
    elif short == "main":
        cls = getattr(module, "TIRLearner", None)
        f = getattr(cls, "save_checkpoint_file", None)
        if f is not None and not getattr(f, "_empty_shim", False):
            def save_checkpoint_file(self, *args, **kwargs):
                # a step that found nothing stores None; 3.0.2 calls None.to_csv
                d = self.working_df_dict
                held = {k: v for k, v in d.items() if v is None}
                for k in held:
                    del d[k]
                try:
                    return f(self, *args, **kwargs)
                finally:
                    d.update(held)
            save_checkpoint_file = functools.wraps(f)(save_checkpoint_file)
            save_checkpoint_file._empty_shim = True
            cls.save_checkpoint_file = save_checkpoint_file
    elif short == "post_processing":
        def make(f):
            def combine_all(df_list, *args, **kwargs):
                kept = [d for d in df_list if not empty(d)]
                return f(kept, *args, **kwargs) if kept else pd.DataFrame()
            return combine_all
        wrap("combine_all", make)


_TIRL_MODULES = {"main", "get_fasta_sequence", "run_TIRvish", "run_GRF", "process_de_novo_result", "post_processing"}


def _is_tirlearner(spec):
    return "TIR-Learner" in (getattr(spec, "origin", None) or "")


class _PatchFinder:
    """Meta-path hook: load a module normally, then patch it (swifter.swifter;
    TIR-Learner's own modules, imported as bare names or as bin.<name>)."""

    def find_spec(self, name, path=None, target=None):
        short = name.rsplit(".", 1)[-1]
        if name != "swifter.swifter" and short not in _TIRL_MODULES:
            return None
        for finder in _sys.meta_path:
            if finder is self or not hasattr(finder, "find_spec"):
                continue
            spec = finder.find_spec(name, path, target)
            if spec is None:
                continue
            if name == "swifter.swifter":
                patch = _patch_swifter
            elif _is_tirlearner(spec):
                patch = _patch_tirlearner
            else:
                return None
            if spec.loader is None or not hasattr(spec.loader, "exec_module"):
                return None
            orig_exec = spec.loader.exec_module

            def exec_module(module, _orig=orig_exec, _patch=patch):
                _orig(module)
                _patch(module)

            spec.loader.exec_module = exec_module
            return spec
        return None


if "swifter.swifter" in _sys.modules:
    _patch_swifter(_sys.modules["swifter.swifter"])
_sys.meta_path.insert(0, _PatchFinder())


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
