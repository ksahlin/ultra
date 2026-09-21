//! Build uLTRA's annotation index: the Rust equivalent of
//! `modules/create_augmented_gene.py` plus `uLTRA`'s `prep_seqs`.
//!
//! This is a transliteration, not a redesign. Where the reference does
//! something surprising it is reproduced and the reason is noted, because the
//! contract is the index's semantic content (PORTING.md Part 5) and a tidier
//! algorithm is a different algorithm.
//!
//! The one deliberate departure is Finding 22, marked DIVERGENCE below.

use std::collections::{BTreeMap, BTreeSet, HashMap};

/// (chr_id, start, stop). The reference keys everything by
/// `array("L",[chr,start,stop]).tobytes()`; we keep the triple and only
/// serialise when we must. Finding 20: that byte form is u64 little-endian on
/// every platform uLTRA runs on.
pub type Key = (u64, u64, u64);

#[derive(Default, Debug, PartialEq)]
pub struct Index {
    pub chr_to_id: BTreeMap<String, u64>,
    pub id_to_chr: BTreeMap<u64, String>,
    pub refs_lengths: BTreeMap<String, u64>,
    /// The reference names in FASTA order. The SAM `@SQ` header is written in
    /// this order -- the reference builds it from `list(refs_lengths.keys())`,
    /// i.e. Python dict insertion order, which is the order `load_reference`
    /// read the fasta. Sorting instead would coincide on SIRV and differ on any
    /// assembly with mixed-form contig names.
    pub ref_order: Vec<String>,
    pub refs_id_lengths: BTreeMap<u64, u64>,

    pub exon_choordinates_to_id: BTreeSet<Key>,
    pub flank_choordinates: BTreeSet<Key>,

    pub segment_id_to_choordinates: BTreeMap<Key, (u64, u64)>,
    pub segment_to_ref: BTreeMap<Key, u64>,
    pub segment_to_gene: BTreeMap<Key, BTreeSet<String>>,
    pub parts_to_segments: BTreeMap<Key, Vec<Key>>,
    pub gene_to_small_segments: BTreeMap<String, Vec<Key>>,

    pub max_intron_chr: BTreeMap<u64, u64>,
    pub all_splice_sites_annotations: BTreeMap<u64, BTreeSet<u64>>,
    pub all_splice_pairs_annotations: BTreeMap<u64, BTreeMap<(u64, u64), BTreeSet<String>>>,
    pub splices_to_transcripts: BTreeMap<u64, BTreeMap<Vec<(u64, u64)>, BTreeSet<String>>>,
    pub transcripts_to_splices: BTreeMap<u64, BTreeMap<String, Vec<(u64, u64)>>>,

    pub ref_part_sequences: BTreeMap<Key, String>,
    pub ref_segment_sequences: BTreeMap<Key, String>,
    pub ref_exon_sequences: BTreeMap<Key, String>,
    pub ref_flank_sequences: BTreeMap<Key, String>,
}

pub struct Params {
    pub flank_size: u64,
    pub small_exon_threshold: u64,
    pub min_segm: u64,
}

/// `add_to_chr_mapping`: ids are assigned 1.. in order of first appearance in
/// the (seqid, start)-sorted exon iteration, so they depend on the annotation's
/// contig ordering and nothing else.
fn chr_id_for(name: &str, chr_to_id: &mut BTreeMap<String, u64>, id_to_chr: &mut BTreeMap<u64, String>) -> u64 {
    if let Some(&id) = chr_to_id.get(name) {
        return id;
    }
    let id = chr_to_id.len() as u64 + 1;
    chr_to_id.insert(name.to_string(), id);
    id_to_chr.insert(id, name.to_string());
    id
}

/// Per-part state accumulated during the exon sweep.
#[derive(Default)]
struct PartAcc {
    canonical_pos: BTreeSet<u64>,
    active_genes: BTreeSet<String>,
    /// (position, is_start) -> exon ids
    pos_to_exon_ids: BTreeMap<(u64, bool), BTreeSet<u32>>,
    choord: (u64, u64),
}

pub fn build(gtf: &crate::gtf::Gtf, refs_lengths: &BTreeMap<String, u64>, p: &Params) -> Index {
    let mut ix = Index::default();
    ix.refs_lengths = refs_lengths.clone();

    let mut exon_id_to_choord: HashMap<u32, (u64, u64)> = HashMap::new();
    // parts keyed by (chr_id, part_counter), kept in creation order
    let mut parts: Vec<((u64, u64), PartAcc)> = Vec::new();
    let mut part_index: HashMap<(u64, u64), usize> = HashMap::new();
    let mut part_counter: u64 = 0;

    let mut prev_chr_id: u64 = 0;
    let mut prev_chr_name = String::new();
    let mut active_start: u64 = 0;
    let mut active_stop: u64 = 0;
    let mut active_genes: BTreeSet<String> = BTreeSet::new();
    let mut cur_chr_name = String::new();
    let mut cur_chr_id: u64 = 0;

    macro_rules! part {
        ($cid:expr) => {{
            let k = ($cid, part_counter);
            let idx = *part_index.entry(k).or_insert_with(|| {
                parts.push((k, PartAcc::default()));
                parts.len() - 1
            });
            &mut parts[idx].1
        }};
    }

    for (i, ex) in gtf.exons.iter().enumerate() {
        let chr_id = chr_id_for(&ex.seqid, &mut ix.chr_to_id, &mut ix.id_to_chr);
        cur_chr_name = ex.seqid.clone();
        cur_chr_id = chr_id;
        exon_id_to_choord.insert(ex.id, (ex.start, ex.stop));
        ix.exon_choordinates_to_id.insert((chr_id, ex.start, ex.stop));
        let genes: BTreeSet<String> = ex.gene_ids.iter().cloned().collect();

        if i == 0 {
            prev_chr_id = chr_id;
            prev_chr_name = ex.seqid.clone();
            active_start = ex.start;
            active_stop = ex.stop;
            active_genes = genes.clone();
            if ex.start > 0 {
                // The reference computes max(0, exon.start - 2*flank_size)
                // from the 1-BASED exon.start, while ex.start is 0-based.
                ix.flank_choordinates.insert((
                    chr_id,
                    (ex.start + 1).saturating_sub(2 * p.flank_size),
                    ex.start,
                ));
            }
            let pa = part!(chr_id);
            pa.canonical_pos.insert(ex.start);
            pa.canonical_pos.insert(ex.stop);
            pa.active_genes.extend(genes.iter().cloned());
            // NOTE: the reference does NOT populate pos_to_exon_ids here, only
            // in the three later branches. Reproduced; see PORTING.md
            // Finding 24.
            continue;
        }

        if chr_id != prev_chr_id {
            // close the previous chromosome's part
            if let Some(&idx) = part_index.get(&(prev_chr_id, part_counter)) {
                parts[idx].1.choord = (active_start, active_stop);
            }
            let chr_length = *refs_lengths
                .get(&prev_chr_name)
                .unwrap_or(&(active_stop + 2 * p.flank_size));
            if active_stop < chr_length {
                ix.flank_choordinates.insert((
                    prev_chr_id,
                    active_stop,
                    std::cmp::min(chr_length, active_stop + 2 * p.flank_size),
                ));
            }
            prev_chr_id = chr_id;
            prev_chr_name = ex.seqid.clone();
            active_start = ex.start;
            active_stop = ex.stop;
            if ex.start > 0 {
                // The reference computes max(0, exon.start - 2*flank_size)
                // from the 1-BASED exon.start, while ex.start is 0-based.
                ix.flank_choordinates.insert((
                    chr_id,
                    (ex.start + 1).saturating_sub(2 * p.flank_size),
                    ex.start,
                ));
            }
            part_counter += 1;
            {
                let pa = part!(chr_id);
                pa.canonical_pos.insert(ex.start);
                pa.canonical_pos.insert(ex.stop);
                pa.active_genes.extend(genes.iter().cloned());
                pa.pos_to_exon_ids.entry((ex.start, true)).or_default().insert(ex.id);
                pa.pos_to_exon_ids.entry((ex.stop, false)).or_default().insert(ex.id);
            }
            active_genes = genes.clone();
        } else if ex.start > active_stop + 20 {
            if let Some(&idx) = part_index.get(&(chr_id, part_counter)) {
                parts[idx].1.choord = (active_start, active_stop);
            }
            // 2*flank_size when this exon shares no gene with the open part
            let segment_size = if genes.is_disjoint(&active_genes) {
                2 * p.flank_size
            } else {
                p.flank_size
            };
            // The reference tests `exon.start - active_stop > 2*segment_size`
            // where exon.start is 1-BASED, while ex.start here is 0-based
            // (exon.start - 1). So the 1-based value is ex.start + 1.
            if (ex.start + 1).saturating_sub(active_stop) > 2 * segment_size {
                let chr_length = *refs_lengths
                    .get(&ex.seqid)
                    .unwrap_or(&(active_stop + segment_size));
                ix.flank_choordinates.insert((
                    chr_id,
                    active_stop,
                    std::cmp::min(chr_length, active_stop + segment_size),
                ));
                ix.flank_choordinates.insert((
                    chr_id,
                    (ex.start + 1).saturating_sub(segment_size),
                    ex.start,
                ));
            } else {
                ix.flank_choordinates.insert((chr_id, active_stop, ex.start));
            }
            active_genes = genes.clone();
            active_start = ex.start;
            active_stop = ex.stop;
            part_counter += 1;
            let pa = part!(chr_id);
            pa.canonical_pos.insert(ex.start);
            pa.canonical_pos.insert(ex.stop);
            pa.active_genes.extend(genes.iter().cloned());
            pa.pos_to_exon_ids.entry((ex.start, true)).or_default().insert(ex.id);
            pa.pos_to_exon_ids.entry((ex.stop, false)).or_default().insert(ex.id);
        } else {
            active_stop = std::cmp::max(active_stop, ex.stop);
            active_genes.extend(genes.iter().cloned());
            let pa = part!(chr_id);
            pa.canonical_pos.insert(ex.start);
            pa.canonical_pos.insert(ex.stop);
            pa.active_genes.extend(genes.iter().cloned());
            pa.pos_to_exon_ids.entry((ex.start, true)).or_default().insert(ex.id);
            pa.pos_to_exon_ids.entry((ex.stop, false)).or_default().insert(ex.id);
        }
    }

    // the very last flank, and the last part's coordinates
    let chr_length = *refs_lengths
        .get(&cur_chr_name)
        .unwrap_or(&(active_stop + 2 * p.flank_size));
    if active_stop < chr_length {
        ix.flank_choordinates.insert((
            cur_chr_id,
            active_stop,
            std::cmp::min(chr_length, active_stop + 2 * p.flank_size),
        ));
    }
    if let Some(&idx) = part_index.get(&(cur_chr_id, part_counter)) {
        parts[idx].1.choord = (active_start, active_stop);
    }

    canonical_segments(&mut ix, &parts, &exon_id_to_choord, p);
    transcript_splices(&mut ix, gtf);
    ix
}

/// `get_canonical_segments`.
fn canonical_segments(
    ix: &mut Index,
    parts: &[((u64, u64), PartAcc)],
    exon_id_to_choord: &HashMap<u32, (u64, u64)>,
    p: &Params,
) {
    for ((chr_id, _pid), pa) in parts.iter() {
        let chr_id = *chr_id;
        let (astart, astop) = pa.choord;
        let part_name: Key = (chr_id, astart, astop);
        let genes = &pa.active_genes;
        let sorted_pos: Vec<u64> = pa.canonical_pos.iter().cloned().collect();
        if sorted_pos.len() < 2 {
            continue;
        }
        let pos_tuples: Vec<(u64, u64)> = sorted_pos
            .windows(2)
            .map(|w| (w[0], w[1]))
            .collect();

        let mut open_starts: BTreeSet<u32> = BTreeSet::new();

        let mut add_segment = |ix: &mut Index, a: u64, b: u64| {
            let name: Key = (chr_id, a, b);
            ix.segment_id_to_choordinates.insert(name, (a, b));
            ix.segment_to_ref.insert(name, chr_id);
            ix.parts_to_segments.entry(part_name).or_default().push(name);
            ix.segment_to_gene.insert(name, genes.clone());
            if b - a <= p.small_exon_threshold {
                for g in genes.iter() {
                    ix.gene_to_small_segments.entry(g.clone()).or_default().push(name);
                }
            }
        };

        for (i, &(p1, p2)) in pos_tuples.iter().enumerate() {
            if let Some(s) = pa.pos_to_exon_ids.get(&(p1, true)) {
                open_starts.extend(s.iter().cloned());
            }
            if let Some(s) = pa.pos_to_exon_ids.get(&(p1, false)) {
                for e in s {
                    open_starts.remove(e);
                }
            }

            if p2 - p1 >= p.min_segm {
                add_segment(ix, p1, p2);
                if let Some(s) = pa.pos_to_exon_ids.get(&(p2, false)) {
                    for e in s {
                        open_starts.remove(e);
                    }
                }
                continue;
            }

            // the segment is too short: try to extend backwards, then forwards
            if i > 0 {
                let mut k = 1usize;
                while i >= k {
                    let cand = pos_tuples[i - k].0;
                    if p2 - cand >= p.min_segm {
                        let name: Key = (chr_id, cand, p2);
                        if !ix.segment_id_to_choordinates.contains_key(&name) {
                            add_segment(ix, cand, p2);
                        }
                        break;
                    }
                    k += 1;
                }
            }
            if i < pos_tuples.len() - 1 {
                let mut k = 1usize;
                while i + k <= pos_tuples.len() - 1 {
                    let cand = pos_tuples[i + k].1;
                    if cand - p1 >= p.min_segm {
                        // NOTE: the reference guards this with
                        //     if (chr_id, p1, cand) not in segment_id_to_choordinates
                        // -- a TUPLE tested against a dict keyed by BYTES, so the
                        // guard is always true and never dedups. Measured: it
                        // produces no duplicates on either corpus anyway, so the
                        // faithful behaviour is simply "always add".
                        // PORTING.md Finding 24.
                        add_segment(ix, p1, cand);
                        break;
                    }
                    k += 1;
                }
            }
            if pos_tuples.len() == 1 {
                add_segment(ix, p1, p2);
            }

            // whole exons spanning this too-short gap
            let mut spanning: BTreeSet<u32> = BTreeSet::new();
            if let Some(s) = pa.pos_to_exon_ids.get(&(p1, true)) {
                spanning.extend(s.iter().cloned());
            }
            if let Some(s) = pa.pos_to_exon_ids.get(&(p2, false)) {
                spanning.extend(s.iter().cloned());
            }
            if spanning.is_empty() {
                spanning = open_starts.clone();
            }
            for e_id in spanning.iter() {
                let (es, ee) = match exon_id_to_choord.get(e_id) {
                    Some(v) => *v,
                    None => continue,
                };
                let name: Key = (chr_id, es, ee);
                if ix.segment_id_to_choordinates.contains_key(&name) {
                    continue;
                }
                ix.segment_id_to_choordinates.insert(name, (es, ee));
                ix.segment_to_ref.insert(name, chr_id);
                ix.parts_to_segments.entry(part_name).or_default().push(name);
                ix.segment_to_gene.insert(name, genes.clone());
                if ee - es <= p.small_exon_threshold {
                    // DIVERGENCE, PORTING.md Finding 22. The reference writes
                    //     add_items(gene_to_small_segments[gene_id], ...)
                    // with a bare `gene_id` left bound by an earlier, unrelated
                    // loop -- it crashes with UnboundLocalError when nothing
                    // bound it, and files the segment under an arbitrary gene
                    // when something did. We implement the evident intent,
                    // matching the two sibling call sites above.
                    for g in genes.iter() {
                        ix.gene_to_small_segments.entry(g.clone()).or_default().push(name);
                    }
                }
            }

            if let Some(s) = pa.pos_to_exon_ids.get(&(p2, false)) {
                for e in s {
                    open_starts.remove(e);
                }
            }
        }
    }
}

fn transcript_splices(ix: &mut Index, gtf: &crate::gtf::Gtf) {
    for t in gtf.transcripts.iter() {
        let chr_id = match ix.chr_to_id.get(&t.seqid) {
            Some(&c) => c,
            // A transcript on a contig with no exons never reaches chr_to_id.
            // The reference would KeyError here; we skip, which cannot differ
            // on any annotation where every transcript has exons.
            None => continue,
        };
        let splices: Vec<(u64, u64)> = t
            .exons
            .windows(2)
            .map(|w| (w[0].1, w[1].0))
            .collect();
        for &(i1, i2) in splices.iter() {
            let gap = i2.saturating_sub(i1);
            let e = ix.max_intron_chr.entry(chr_id).or_insert(0);
            if gap > *e {
                *e = gap;
            }
        }
        ix.splices_to_transcripts
            .entry(chr_id)
            .or_default()
            .entry(splices.clone())
            .or_default()
            .insert(t.id.clone());
        for &(s1, s2) in splices.iter() {
            ix.all_splice_pairs_annotations
                .entry(chr_id)
                .or_default()
                .entry((s1, s2))
                .or_default()
                .insert(t.id.clone());
            let sites = ix.all_splice_sites_annotations.entry(chr_id).or_default();
            sites.insert(s1);
            sites.insert(s2);
        }
        if let (Some(first), Some(last)) = (t.exons.first(), t.exons.last()) {
            let sites = ix.all_splice_sites_annotations.entry(chr_id).or_default();
            sites.insert(first.0);
            sites.insert(last.1);
        }
    }
    // transcripts_to_splices is the inversion of splices_to_transcripts
    for (chr_id, m) in ix.splices_to_transcripts.iter() {
        let dest = ix.transcripts_to_splices.entry(*chr_id).or_default();
        for (splices, tids) in m.iter() {
            for tid in tids {
                dest.insert(tid.clone(), splices.clone());
            }
        }
    }
}
