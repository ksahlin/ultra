
from collections import defaultdict

from modules.help_functions import readfq

def reverse_mapping(d):
    d_dict = defaultdict(list)
    for k,v in d.items():
        for i in v:
            d_dict[i].append(k)
    return d_dict

def dd_set(): # top level function declaration needed for multiprocessing
    return defaultdict(set)

def create_graph_from_exon_parts(db, min_mem): 
    """
        We need to link parts --> exons and exons --> transcripts
    """
    # print(dir(db))
    # genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_to_gene = {} # exon_id : [gene_id ]
    gene_to_small_exons = {} # gene_id : [exon_id ]

    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dd_set)
    all_splice_pairs_annotations = defaultdict(dd_set)
    all_splice_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)

    parts_to_exons = defaultdict(dd_set)
    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        # if exon.seqid.isdigit() or exon.seqid == 'X' or exon.seqid == 'Y':
        #     chr_id = 'chr'+ exon.seqid
        # elif exon.seqid == 'MT':
        #     chr_id = 'chrM'
        # else:
        chr_id = exon.seqid

        exons_to_ref[exon.id] = chr_id
        # print(dir(exon))
        # print(exon.attributes["gene_id"])
        exon_to_gene[exon.id] = exon.attributes["gene_id"]
        exon_id_to_choordinates[exon.id] = (exon.start - 1, exon.stop)
        # creating the augmentation
        if i == 0: # initialization
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)

        if chr_id != prev_seq_id: # switching chromosomes
            parts_to_exons[prev_seq_id][(active_start, active_stop)] = active_exons
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)           

        if exon.start - 1 > active_stop:
            parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
            active_exons = set()
            active_exons.add(exon.id)

            active_start = exon.start - 1
            active_stop = exon.stop

        else:
            active_exons.add(exon.id)    
            active_stop = max(active_stop, exon.stop)

        assert active_start <= exon.start - 1

    parts_to_exons[chr_id][(active_start, active_stop)] = active_exons


    for transcript in db.features_of_type('transcript', order_by='seqid'): #db.children(gene, featuretype='transcript', order_by='start'):
        # if transcript.seqid.isdigit() or transcript.seqid == 'X' or transcript.seqid == 'Y':
        #     chr_id = 'chr'+ transcript.seqid
        # elif transcript.seqid == 'MT':
        #     chr_id = 'chrM'
        # else:
        chr_id = transcript.seqid

        # print("here", transcript.id,  str(chr_id))  
        transcript_exons = []
        for exon in db.children(transcript, featuretype='exon', order_by='start'):
            transcript_exons.append( (exon.start-1, exon.stop) )
        internal_transcript_splices = [ (e1[1],e2[0]) for e1, e2 in zip(transcript_exons[:-1],transcript_exons[1:])]
        
        # internal transcript splices
        splices_to_transcripts[chr_id][ tuple(internal_transcript_splices)].add(transcript.id)
        for site1, site2 in internal_transcript_splices:
            all_splice_pairs_annotations[str(chr_id)][(site1, site2)].add( transcript.id )
            # if site2 == 3105:
            #     print('LOOOL',chr_id)
            #     ccc.add(chr_id)

            all_splice_sites_annotations[str(chr_id)].add(site1)
            all_splice_sites_annotations[str(chr_id)].add(site2)
        
        # add start and end splice to all_splice_sites_annotations 
        if transcript_exons:
            all_splice_sites_annotations[str(chr_id)].add(transcript_exons[0][0])
            all_splice_sites_annotations[str(chr_id)].add(transcript_exons[-1][-1])
        else:
            print("Something is wrong with transcript annotation: {0} on gene: {1}, and could not be added. Check that the gene ID and transcript ID is not the same!".format(transcript.id, transcript.seqid))
            # sys.exit()

    # transcripts_to_splices = reverse_mapping(splices_to_transcripts)
    transcripts_to_splices = defaultdict(dict)
    for chr_id, chr_sp_sites_dict in splices_to_transcripts.items():
        for unique_sp_sites, tr_ids in chr_sp_sites_dict.items():
            for tr_id in tr_ids:
                transcripts_to_splices[chr_id][tr_id] = unique_sp_sites


    for gene in db.features_of_type('gene', order_by='seqid'):
        gene_to_small_exons[gene.id] = []
        for exon in  db.children(gene, featuretype='exon', order_by='start'):
            if exon.stop - exon.start < 50:
                gene_to_small_exons[gene.id].append(exon.id)

    # print(gene_to_small_exons)
    # for g in gene_to_small_exons:
    #     for e in gene_to_small_exons[g]:
    #         print(exon_id_to_choordinates[e])
    # sys.exit()
    return  exons_to_ref, parts_to_exons, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, exon_id_to_choordinates, \
            exon_to_gene, gene_to_small_exons



def get_part_sequences_from_choordinates(parts_to_exons, refs):
    segments = {}
    for chr_id in parts_to_exons:
        if chr_id not in refs:
            continue
        else:
            parts_instance = parts_to_exons[chr_id]
            # chromosome = genes_to_ref[chr_id]
            segments[chr_id] = {}
            for part in parts_instance:
                start,stop = part[0], part[1]
                seq = refs[chr_id][start : stop] 
                segments[chr_id][part] = seq
    return segments

def get_exon_sequences_from_choordinates(exon_id_to_choordinates, exons_to_ref, refs):
    exon_sequences = defaultdict(dict)
    for exon_id in exon_id_to_choordinates:
        start,stop = exon_id_to_choordinates[exon_id]
        chr_id = exons_to_ref[exon_id]
        if chr_id not in refs:
            continue
        else:
            seq = refs[chr_id][start : stop] 
            exon_sequences[chr_id][(start,stop)] = seq
    # print(segments)
    return exon_sequences



# Functions for masking overly abundant kmers 

import itertools
from collections import defaultdict, deque
def kmer_counter(ref_part_sequences, kmer_size):
    count = defaultdict(int)
    position_count = defaultdict(list)
    it_count = 0
    for chr_id  in ref_part_sequences:
        for part, seq in ref_part_sequences[chr_id].items():
            read_kmers = deque([seq[i:i+kmer_size] for i in range(len(seq) - kmer_size + 1)])
            for i, kmer in enumerate(read_kmers):
                count[kmer] += 1
                it_count += 1

                if it_count % 10000000 == 0:
                    print(int(it_count/10000000.0), " million kmers processed.")
                    # clean every 1M kmers
                    if it_count % 20000000 == 0:
                        for kmer in list(count.keys()):
                            if count[kmer] == 1:
                                del count[kmer]
                        print(len(count),"kmers stored")
    return count


def mask_refs(ref_part_sequences, to_mask, kmer_size):
    mask_counter = 0
    tot_counter = 0
    for chr_id  in ref_part_sequences:
        for part, seq in ref_part_sequences[chr_id].items():
            tot_counter += 1
            read_kmers = [seq[i:i+kmer_size] for i in range(len(seq) - kmer_size + 1 )]
            seq_masked = []
            has_been_modified = False
            if read_kmers:
                # separately treat first kmer in every seq
                if read_kmers[0] in to_mask:
                    seq_masked.append(read_kmers[0][:-1]+'N')
                    has_been_modified = True
                else:
                    seq_masked.append(read_kmers[0])

                for i, kmer in enumerate(read_kmers[1:]):
                    if kmer in to_mask:
                        seq_masked.append("N") # mask the last base
                        has_been_modified = True
                    else:
                        seq_masked.append(kmer[-1]) # add the true last base
                
                if has_been_modified:
                    seq_masked = "".join([s for s in seq_masked])
                    # print("masking", seq, "to", seq_masked) 
                    ref_part_sequences[chr_id][part] = seq_masked
                    mask_counter += 1
    print(mask_counter, "{0} out of {1} sequences has been modified in masking step.".format(mask_counter, tot_counter))

def mask_abundant_kmers(ref_part_sequences, kmer_size, mask_threshold):
    DBG = kmer_counter(ref_part_sequences, kmer_size)
    n = float(len(DBG))
    print(n, "Unique kmers in reference part sequences with abundance > 1")
    to_mask = set()
    for i, (kmer, abundance) in enumerate(sorted(DBG.items(), key=lambda x: x[1], reverse=True)):
        if abundance >= mask_threshold:
            to_mask.add(kmer) 
            print(kmer, abundance)
        else:
            break

    mask_refs(ref_part_sequences, to_mask, kmer_size)
    # print(len(to_mask), "kemrs masked.")
    # print(sorted([DBG[kmer] for kmer in DBG], reverse=True))





