//! The on-disk index must round-trip losslessly.
//!
//! The port's index format exists because the reference's cannot be reproduced
//! (PORTING.md Finding 3) -- so nothing external validates it, and a silent
//! truncation here would look like an alignment bug much later. This builds a
//! real index from the repository's own fixtures, writes it, reads it back and
//! compares every structure.

#[path = "../src/fasta.rs"]
mod fasta;
#[path = "../src/gtf.rs"]
mod gtf;
#[path = "../src/index.rs"]
mod index;
#[path = "../src/indexio.rs"]
mod indexio;

#[test]
fn index_survives_a_write_and_read() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let refs = fasta::read(root.join("test/SIRV_genes.fasta").to_str().unwrap()).unwrap();
    let mut refs_lengths = std::collections::BTreeMap::new();
    for r in &refs {
        refs_lengths.insert(r.name.clone(), r.seq.len() as u64);
    }
    let parsed = gtf::parse(root.join("test/SIRV_genes_C_170612a.gtf").to_str().unwrap()).unwrap();
    let mut ix = index::build(&parsed, &refs_lengths, &index::Params {
        flank_size: 1000, small_exon_threshold: 200, min_segm: 25,
    });
    // a couple of sequence maps so the string paths are exercised too
    for (i, k) in ix.segment_id_to_choordinates.keys().cloned().enumerate().take(50) {
        ix.ref_segment_sequences.insert(k, "ACGT".repeat(i % 7 + 1));
    }
    for k in ix.flank_choordinates.iter().cloned().take(20) {
        ix.ref_flank_sequences.insert(k, "TTTTGGGG".into());
    }

    let dir = std::env::temp_dir().join(format!("ultra-idx-rt-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("ultra.idx");
    indexio::write(&ix, &p).unwrap();
    let back = indexio::read(&p).unwrap();
    let _ = std::fs::remove_dir_all(&dir);

    assert_eq!(back.chr_to_id, ix.chr_to_id, "chr_to_id");
    assert_eq!(back.id_to_chr, ix.id_to_chr, "id_to_chr");
    assert_eq!(back.refs_lengths, ix.refs_lengths, "refs_lengths");
    assert!(!back.refs_lengths.is_empty(), "refs_lengths round-tripped empty");
    // ref_order is populated by the caller, not by index::build. This test
    // deliberately leaves it unset, which is how the silent loss of every
    // reference length was caught; the reader must still recover them all.
    assert_eq!(
        back.ref_order.iter().cloned().collect::<std::collections::BTreeSet<_>>(),
        ix.refs_lengths.keys().cloned().collect::<std::collections::BTreeSet<_>>(),
        "ref_order must cover every reference"
    );
    assert_eq!(back.refs_id_lengths, ix.refs_id_lengths, "refs_id_lengths");
    assert_eq!(back.max_intron_chr, ix.max_intron_chr, "max_intron_chr");
    assert_eq!(back.parts_to_segments, ix.parts_to_segments, "parts_to_segments");
    assert_eq!(back.segment_to_gene, ix.segment_to_gene, "segment_to_gene");
    assert_eq!(back.gene_to_small_segments, ix.gene_to_small_segments, "gene_to_small_segments");
    assert_eq!(back.all_splice_sites_annotations, ix.all_splice_sites_annotations, "splice sites");
    assert_eq!(back.all_splice_pairs_annotations, ix.all_splice_pairs_annotations, "splice pairs");
    assert_eq!(back.transcripts_to_splices, ix.transcripts_to_splices, "transcripts_to_splices");
    assert_eq!(back.ref_segment_sequences, ix.ref_segment_sequences, "segment sequences");
    assert_eq!(back.ref_flank_sequences, ix.ref_flank_sequences, "flank sequences");

    assert!(!back.parts_to_segments.is_empty(), "round-tripped an empty index");
    println!(
        "index round-trip: {} parts, {} segments, {} genes, {} chromosomes",
        back.parts_to_segments.len(), back.segment_to_gene.len(),
        back.gene_to_small_segments.len(), back.chr_to_id.len()
    );
}

#[test]
fn a_foreign_file_is_rejected_rather_than_misread() {
    let dir = std::env::temp_dir().join(format!("ultra-idx-bad-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("ultra.idx");
    std::fs::write(&p, b"this is not an index at all, it is a text file").unwrap();
    let e = indexio::read(&p);
    let _ = std::fs::remove_dir_all(&dir);
    assert!(e.is_err(), "a non-index file must be refused, not parsed as garbage");
    let msg = format!("{}", e.unwrap_err());
    assert!(msg.contains("rebuild"), "the error should tell the user what to do, got: {msg}");
}
