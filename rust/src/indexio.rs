//! On-disk index format for the port.
//!
//! The reference stores 20 Python pickles plus a gffutils sqlite database.
//! The port cannot reproduce those (PORTING.md Finding 3) and does not try: it
//! writes one binary file, `ultra.idx`, with an explicit little-endian layout.
//!
//! Why hand-rolled rather than a serialisation crate: the format is a dozen
//! maps of fixed-shape keys, the encoding is the same u64 little-endian the
//! reference's `array("L")` already uses on every platform uLTRA runs on
//! (Finding 20), and it keeps the dependency list at the two C libraries we
//! actually need.
//!
//! Layout: a magic and version, then each structure as a u64 count followed by
//! its entries. Strings are u32 length + bytes. Keys are three u64s.

use crate::index::{Index, Key};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufReader, BufWriter, Read, Write};

const MAGIC: &[u8; 8] = b"uLTRAidx";
const VERSION: u32 = 1;

fn wu64<W: Write>(w: &mut W, v: u64) -> std::io::Result<()> { w.write_all(&v.to_le_bytes()) }
fn wi64<W: Write>(w: &mut W, v: i64) -> std::io::Result<()> { w.write_all(&v.to_le_bytes()) }
fn wstr<W: Write>(w: &mut W, s: &str) -> std::io::Result<()> {
    w.write_all(&(s.len() as u32).to_le_bytes())?;
    w.write_all(s.as_bytes())
}
fn wkey<W: Write>(w: &mut W, k: &Key) -> std::io::Result<()> {
    wu64(w, k.0)?; wu64(w, k.1)?; wu64(w, k.2)
}

fn ru64<R: Read>(r: &mut R) -> std::io::Result<u64> {
    let mut b = [0u8; 8]; r.read_exact(&mut b)?; Ok(u64::from_le_bytes(b))
}
fn ri64<R: Read>(r: &mut R) -> std::io::Result<i64> {
    let mut b = [0u8; 8]; r.read_exact(&mut b)?; Ok(i64::from_le_bytes(b))
}
fn rstr<R: Read>(r: &mut R) -> std::io::Result<String> {
    let mut b = [0u8; 4]; r.read_exact(&mut b)?;
    let n = u32::from_le_bytes(b) as usize;
    let mut v = vec![0u8; n]; r.read_exact(&mut v)?;
    Ok(String::from_utf8_lossy(&v).into_owned())
}
fn rkey<R: Read>(r: &mut R) -> std::io::Result<Key> {
    Ok((ru64(r)?, ru64(r)?, ru64(r)?))
}

fn wseqmap<W: Write>(w: &mut W, m: &BTreeMap<Key, String>) -> std::io::Result<()> {
    wu64(w, m.len() as u64)?;
    for (k, v) in m { wkey(w, k)?; wstr(w, v)?; }
    Ok(())
}
fn rseqmap<R: Read>(r: &mut R) -> std::io::Result<BTreeMap<Key, String>> {
    let n = ru64(r)?;
    let mut m = BTreeMap::new();
    for _ in 0..n { let k = rkey(r)?; let v = rstr(r)?; m.insert(k, v); }
    Ok(m)
}

pub fn write(ix: &Index, path: &std::path::Path) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    w.write_all(MAGIC)?;
    w.write_all(&VERSION.to_le_bytes())?;

    wu64(&mut w, ix.chr_to_id.len() as u64)?;
    for (name, id) in &ix.chr_to_id { wstr(&mut w, name)?; wu64(&mut w, *id)?; }

    // refs_lengths, written in FASTA order so the @SQ header survives.
    // ref_order is set by the caller, not by index::build, so it can be short
    // or empty; drive the write off refs_lengths and use ref_order only to
    // order it, or an index built without it loses every reference silently.
    let mut order: Vec<&String> = Vec::with_capacity(ix.refs_lengths.len());
    for name in &ix.ref_order {
        if ix.refs_lengths.contains_key(name) { order.push(name); }
    }
    let seen: std::collections::BTreeSet<&String> = order.iter().cloned().collect();
    for name in ix.refs_lengths.keys() {
        if !seen.contains(name) { order.push(name); }
    }
    wu64(&mut w, order.len() as u64)?;
    for name in order {
        wstr(&mut w, name)?;
        wu64(&mut w, ix.refs_lengths[name])?;
    }

    wu64(&mut w, ix.refs_id_lengths.len() as u64)?;
    for (id, l) in &ix.refs_id_lengths { wu64(&mut w, *id)?; wu64(&mut w, *l)?; }

    wu64(&mut w, ix.max_intron_chr.len() as u64)?;
    for (c, m) in &ix.max_intron_chr { wu64(&mut w, *c)?; wu64(&mut w, *m)?; }

    wu64(&mut w, ix.parts_to_segments.len() as u64)?;
    for (k, v) in &ix.parts_to_segments {
        wkey(&mut w, k)?; wu64(&mut w, v.len() as u64)?;
        for s in v { wkey(&mut w, s)?; }
    }

    wu64(&mut w, ix.segment_to_gene.len() as u64)?;
    for (k, g) in &ix.segment_to_gene {
        wkey(&mut w, k)?; wu64(&mut w, g.len() as u64)?;
        for name in g { wstr(&mut w, name)?; }
    }

    wu64(&mut w, ix.gene_to_small_segments.len() as u64)?;
    for (g, v) in &ix.gene_to_small_segments {
        wstr(&mut w, g)?; wu64(&mut w, v.len() as u64)?;
        for s in v { wkey(&mut w, s)?; }
    }

    wu64(&mut w, ix.all_splice_sites_annotations.len() as u64)?;
    for (c, sites) in &ix.all_splice_sites_annotations {
        wu64(&mut w, *c)?; wu64(&mut w, sites.len() as u64)?;
        for s in sites { wu64(&mut w, *s)?; }
    }

    wu64(&mut w, ix.all_splice_pairs_annotations.len() as u64)?;
    for (c, pairs) in &ix.all_splice_pairs_annotations {
        wu64(&mut w, *c)?; wu64(&mut w, pairs.len() as u64)?;
        for ((a, b), trs) in pairs {
            wu64(&mut w, *a)?; wu64(&mut w, *b)?; wu64(&mut w, trs.len() as u64)?;
            for t in trs { wstr(&mut w, t)?; }
        }
    }

    wu64(&mut w, ix.transcripts_to_splices.len() as u64)?;
    for (c, m) in &ix.transcripts_to_splices {
        wu64(&mut w, *c)?; wu64(&mut w, m.len() as u64)?;
        for (tid, sp) in m {
            wstr(&mut w, tid)?; wu64(&mut w, sp.len() as u64)?;
            for (a, b) in sp { wu64(&mut w, *a)?; wu64(&mut w, *b)?; }
        }
    }

    wseqmap(&mut w, &ix.ref_part_sequences)?;
    wseqmap(&mut w, &ix.ref_segment_sequences)?;
    wseqmap(&mut w, &ix.ref_exon_sequences)?;
    wseqmap(&mut w, &ix.ref_flank_sequences)?;
    w.flush()
}

pub fn read(path: &std::path::Path) -> std::io::Result<Index> {
    let mut r = BufReader::new(std::fs::File::open(path)?);
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
            "not a uLTRA index (bad magic); rebuild it with `uLTRA index`"));
    }
    let mut vb = [0u8; 4]; r.read_exact(&mut vb)?;
    let v = u32::from_le_bytes(vb);
    if v != VERSION {
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidData,
            format!("index format version {v}, this build expects {VERSION}; rebuild it")));
    }

    let mut ix = Index::default();
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let name = rstr(&mut r)?; let id = ru64(&mut r)?;
        ix.id_to_chr.insert(id, name.clone());
        ix.chr_to_id.insert(name, id);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let name = rstr(&mut r)?; let l = ru64(&mut r)?;
        ix.ref_order.push(name.clone());
        ix.refs_lengths.insert(name, l);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n { let id = ru64(&mut r)?; let l = ru64(&mut r)?; ix.refs_id_lengths.insert(id, l); }
    let n = ru64(&mut r)?;
    for _ in 0..n { let c = ru64(&mut r)?; let m = ru64(&mut r)?; ix.max_intron_chr.insert(c, m); }

    let n = ru64(&mut r)?;
    for _ in 0..n {
        let k = rkey(&mut r)?; let m = ru64(&mut r)?;
        let mut v = Vec::with_capacity(m as usize);
        for _ in 0..m { v.push(rkey(&mut r)?); }
        ix.parts_to_segments.insert(k, v);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let k = rkey(&mut r)?; let m = ru64(&mut r)?;
        let mut g = BTreeSet::new();
        for _ in 0..m { g.insert(rstr(&mut r)?); }
        ix.segment_to_gene.insert(k, g);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let g = rstr(&mut r)?; let m = ru64(&mut r)?;
        let mut v = Vec::with_capacity(m as usize);
        for _ in 0..m { v.push(rkey(&mut r)?); }
        ix.gene_to_small_segments.insert(g, v);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let c = ru64(&mut r)?; let m = ru64(&mut r)?;
        let mut s = BTreeSet::new();
        for _ in 0..m { s.insert(ru64(&mut r)?); }
        ix.all_splice_sites_annotations.insert(c, s);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let c = ru64(&mut r)?; let m = ru64(&mut r)?;
        let mut pairs = BTreeMap::new();
        for _ in 0..m {
            let a = ru64(&mut r)?; let b = ru64(&mut r)?; let t = ru64(&mut r)?;
            let mut trs = BTreeSet::new();
            for _ in 0..t { trs.insert(rstr(&mut r)?); }
            pairs.insert((a, b), trs);
        }
        ix.all_splice_pairs_annotations.insert(c, pairs);
    }
    let n = ru64(&mut r)?;
    for _ in 0..n {
        let c = ru64(&mut r)?; let m = ru64(&mut r)?;
        let mut tm = BTreeMap::new();
        for _ in 0..m {
            let tid = rstr(&mut r)?; let k = ru64(&mut r)?;
            let mut sp = Vec::with_capacity(k as usize);
            for _ in 0..k { sp.push((ru64(&mut r)?, ru64(&mut r)?)); }
            tm.insert(tid, sp);
        }
        ix.transcripts_to_splices.insert(c, tm);
    }

    ix.ref_part_sequences = rseqmap(&mut r)?;
    ix.ref_segment_sequences = rseqmap(&mut r)?;
    ix.ref_exon_sequences = rseqmap(&mut r)?;
    ix.ref_flank_sequences = rseqmap(&mut r)?;
    let _ = (wi64::<Vec<u8>>, ri64::<&[u8]>);
    Ok(ix)
}
