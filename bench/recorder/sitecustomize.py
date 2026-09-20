"""Record every aligner call, in EVERY process.

uLTRA does its alignment in worker processes spawned by `pc.Managers`, and on
macOS/Python 3.8+ `spawn` means those children re-import everything -- so
patching the parent records nothing. `sitecustomize` is imported during each
interpreter's site initialisation, including spawned children, which is why the
hook lives here.

Activated by putting this directory on PYTHONPATH and setting:
    ULTRA_REC_DIR   directory to write <pid>.jsonl into
    ULTRA_REC_LIMIT optional per-process cap
"""
import atexit
import json
import os

_dir = os.environ.get("ULTRA_REC_DIR")
if _dir:
    os.makedirs(_dir, exist_ok=True)
    _limit = int(os.environ.get("ULTRA_REC_LIMIT", "0"))
    _fh = open(os.path.join(_dir, f"{os.getpid()}.jsonl"), "w")
    _n = [0]

    def _rec(kind, payload):
        if _limit and _n[0] >= _limit:
            return
        _n[0] += 1
        _fh.write(json.dumps({"kind": kind, **payload}, sort_keys=True) + "\n")

    @atexit.register
    def _close():
        try:
            _fh.close()
        except Exception:
            pass

    def _install():
        import edlib
        if getattr(edlib, "_ultra_recorded", False):
            return
        _align = edlib.align

        def align(query, target, **kw):
            r = _align(query, target, **kw)
            _rec("edlib.align", {
                "query": query, "target": target,
                "mode": kw.get("mode"), "task": kw.get("task"), "k": kw.get("k", -1),
                "editDistance": r.get("editDistance"),
                "locations": r.get("locations"),
                "cigar": r.get("cigar"),
            })
            return r

        edlib.align = align
        edlib._ultra_recorded = True

        # parasail is wrapped at uLTRA's own function, so the _16 -> _32
        # saturation fallback and the derived alignment strings are captured
        # as one record rather than as a raw library call.
        try:
            from modules import help_functions as hf
        except Exception:
            return
        if getattr(hf, "_ultra_recorded", False):
            return
        _pa = hf.parasail_alignment

        def parasail_alignment(s1, s2, **kw):
            r = _pa(s1, s2, **kw)
            read_aln, ref_aln, cigar_string, cigar_tuples, score = r
            _rec("parasail_alignment", {
                "s1": s1, "s2": s2, "kw": dict(kw),
                "read_aln": read_aln, "ref_aln": ref_aln,
                "cigar": cigar_string, "score": score,
            })
            return r

        hf.parasail_alignment = parasail_alignment
        hf._ultra_recorded = True

    def _install_chaining():
        """Record the colinear solver's inputs and outputs.

        These are the chaining core -- the bulk of uLTRA's own 59% of wall
        clock. Mems are namedtuples; they serialise as their field order
        (x, y, c, d, val, j, exon_part_id) so a replay can rebuild them.
        """
        from modules import colinear_solver as cs
        if getattr(cs, "_ultra_recorded", False):
            return

        def _mem(m):
            return [m.x, m.y, m.c, m.d, m.val, m.j, m.exon_part_id]

        def _sols(solutions):
            return [[_mem(m) for m in sol] for sol in solutions]

        _rc = cs.read_coverage
        def read_coverage(mems, max_intron):
            out = _rc(mems, max_intron)
            _rec("read_coverage", {
                "mems": [_mem(m) for m in mems], "max_intron": max_intron,
                "solutions": _sols(out[0]), "value": out[1],
            })
            return out
        cs.read_coverage = read_coverage

        _nl = cs.n_logn_read_coverage
        def n_logn_read_coverage(mems):
            out = _nl(mems)
            _rec("n_logn_read_coverage", {
                "mems": [_mem(m) for m in mems],
                "solutions": _sols(out[0]), "value": out[1],
            })
            return out
        cs.n_logn_read_coverage = n_logn_read_coverage

        def _mam(m):
            return [m.x, m.y, m.c, m.d, m.val, m.j, m.min_segment_length, m.mam_id, m.ref_chr_id]

        _rcm = cs.read_coverage_mam_score
        def read_coverage_mam_score(mams, overlap_threshold=20):
            out = _rcm(mams, overlap_threshold)
            _rec("read_coverage_mam_score", {
                "mams": [_mam(m) for m in mams], "overlap_threshold": overlap_threshold,
                "solution": [_mam(m) for m in out[0]], "value": out[1], "unique": out[2],
            })
            return out
        cs.read_coverage_mam_score = read_coverage_mam_score

        _nlm = cs.n_logn_read_coverage_mams
        def n_logn_read_coverage_mams(mams, overlap_threshold=5):
            out = _nlm(mams, overlap_threshold)
            _rec("n_logn_read_coverage_mams", {
                "mams": [_mam(m) for m in mams], "overlap_threshold": overlap_threshold,
                "solution": [_mam(m) for m in out[0]], "value": out[1], "unique": out[2],
            })
            return out
        cs.n_logn_read_coverage_mams = n_logn_read_coverage_mams

        import modules.align as am
        am.colinear_solver.read_coverage = read_coverage
        am.colinear_solver.n_logn_read_coverage = n_logn_read_coverage
        try:
            from modules import classify_read_with_mams as cw
            cw.colinear_solver.read_coverage_mam_score = read_coverage_mam_score
            cw.colinear_solver.n_logn_read_coverage_mams = n_logn_read_coverage_mams
        except Exception:
            pass
        cs._ultra_recorded = True

    def _install_mamplumbing():
        from modules import classify_read_with_mams as cw
        if getattr(cw, "_ultra_recorded_plumb", False):
            return
        _g = cw.get_unique_exon_and_flank_locations
        def get_unique_exon_and_flank_locations(solution, parts_to_segments):
            out = _g(solution, parts_to_segments)
            shl, psh, fhl, pfh, c2e, fps, lps = out
            _rec("get_unique_exon_and_flank_locations", {
                "solution": [[m.x, m.y, m.c, m.d, m.val, m.j, m.exon_part_id] for m in solution],
                "segment_hit_locations": [list(t) for t in shl],
                "flank_hit_locations": [list(t) for t in fhl],
                "partial_segment_hit_locations": {repr(k): v for k, v in psh.items()},
                "partial_flank_hit_locations": {repr(k): v for k, v in pfh.items()},
                "first_part_stop": fps, "last_part_start": lps,
            })
            return out
        cw.get_unique_exon_and_flank_locations = get_unique_exon_and_flank_locations

        _a = cw.add_segment_to_mam
        def add_segment_to_mam(read_seq, ref_chr_id, exon_seq, e_start, e_stop, segm_id,
                               mam_instance, min_acc, annot_label):
            before = len(mam_instance)
            _a(read_seq, ref_chr_id, exon_seq, e_start, e_stop, segm_id,
               mam_instance, min_acc, annot_label)
            added = mam_instance[before:]
            _rec("add_segment_to_mam", {
                "read_seq": read_seq, "ref_chr_id": ref_chr_id, "exon_seq": exon_seq,
                "e_start": e_start, "e_stop": e_stop, "segm_id": str(segm_id),
                "min_acc": min_acc, "annot_label": annot_label,
                "added": [[m.x, m.y, m.c, m.d, m.val, m.j, m.min_segment_length,
                           m.mam_id, m.ref_chr_id] for m in added],
            })
        cw.add_segment_to_mam = add_segment_to_mam
        cw._ultra_recorded_plumb = True

    # modules/ is only importable once uLTRA has set up sys.path, so defer the
    # install until the first import of the package rather than doing it here.
    import importlib.abc
    import importlib.machinery
    import sys

    class _Hook(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname in ("modules.help_functions", "modules.align", "modules.classify_read_with_mams"):
                spec = importlib.machinery.PathFinder.find_spec(fullname, path)
                if spec and spec.loader:
                    orig_exec = spec.loader.exec_module
                    want = fullname

                    def exec_module(module, _orig=orig_exec, _w=want):
                        _orig(module)
                        if _w == "modules.help_functions":
                            _install()
                        elif _w == "modules.classify_read_with_mams":
                            _install_mamplumbing()
                        else:
                            _install_chaining()

                    spec.loader.exec_module = exec_module
                return spec
            return None

    sys.meta_path.insert(0, _Hook())
