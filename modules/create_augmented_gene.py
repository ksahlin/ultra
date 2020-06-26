
from collections import defaultdict

import intervaltree

from modules.help_functions import readfq

def reverse_mapping(d):
    d_dict = defaultdict(list)
    for k,v in d.items():
        for i in v:
            d_dict[i].append(k)
    return d_dict

def dd_set(): # top level function declaration needed for multiprocessing
    return defaultdict(set)

def create_graph_from_exon_parts(db, flank_size, small_exon_threshold): 
    """
        We need to link parts --> exons and exons --> transcripts
    """
    # print(dir(db))
    # genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_to_gene = {} # exon_id : [gene_id ]

    flanks_to_gene2 = defaultdict(dict)  
    total_flanks2 = 0

    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dd_set)
    all_splice_pairs_annotations = defaultdict(dd_set)
    all_splice_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)
    part_intervals = defaultdict(intervaltree.IntervalTree)
    parts_to_exons = defaultdict(dd_set)
    max_intron_chr = defaultdict(int)
    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):

        # if i > 0:
        #     if  exon.stop - active_stop > 100000:
        #         check so that  its the  same gene
        #         print("intron size:", exon.stop - active_stop, exon.seqid,prev_seq_id)
        # print(exon.seqid, exon.start, exon.stop, exon.featuretype,dir(exon))
        # import sys
        # sys.exit()
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
            # adding very first flank on chromosome
            # print(max(0, exon.start - 1000), exon.start - 1)
            flanks_to_gene2[chr_id][(max(0, exon.start - flank_size), exon.start - 1)] = "flank_{0}".format(i)
            total_flanks2 += 1

        if chr_id != prev_seq_id: # switching chromosomes
            parts_to_exons[prev_seq_id][(active_start, active_stop)] = active_exons
            part_intervals[prev_seq_id].addi(active_start, active_stop, None)
            
            # adding very last flank on chromosome
            flanks_to_gene2[prev_seq_id][(max(0, active_stop), active_stop + flank_size)] = "flank_{0}".format(i)
            total_flanks2 += 1
            # print(max(0, active_stop), active_stop + 1000)
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)     


        if exon.start - 1 > active_stop:
            parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
            part_intervals[prev_seq_id].addi(active_start, active_stop, None)
            if exon.start - active_stop > 2*flank_size:
                flanks_to_gene2[chr_id][(max(0, active_stop), active_stop + flank_size)] = "flank_{0}_1".format(i)
                flanks_to_gene2[chr_id][(max(0, exon.start - flank_size), exon.start - 1)] = "flank_{0}_2".format(i)
                total_flanks2 += 2
                # print(max(0, active_stop), active_stop + 1000)
                # print(max(0, exon.start - 1000), exon.start - 1)

            else: # add the whole intron
                flanks_to_gene2[chr_id][(max(0, active_stop), exon.start - 1)] = "flank_{0}".format(i)
                total_flanks2 += 1

            active_exons = set()
            active_exons.add(exon.id)

            active_start = exon.start - 1
            active_stop = exon.stop

        else:
            active_exons.add(exon.id)    
            active_stop = max(active_stop, exon.stop)

        assert active_start <= exon.start - 1

    print("total_flanks2:", total_flanks2)

    parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
    part_intervals[prev_seq_id].addi(active_start, active_stop, None)


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
        
        for i1,i2 in internal_transcript_splices:
            if i2 - i1 > max_intron_chr[chr_id]:
                max_intron_chr[chr_id] = i2 - i1
                # print("intron size:", i2 - i1, chr_id)

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

    flanks_to_gene = defaultdict(dict)   
    gene_to_small_exons = {} # gene_id : [exon_id ]
    flanks_not_overlapping = 0
    total_flanks = 0
    for gene in db.features_of_type('gene', order_by='seqid'):
        gene_to_small_exons[gene.id] = []
        exons_list = [exon for exon in  db.children(gene, featuretype='exon', order_by='start')]
        chr_id = gene.seqid
        if exons_list:
            ovl = part_intervals[chr_id].overlaps(max(0, exons_list[0].start - flank_size), exons_list[0].start - 1)
            if not ovl:
                flanks_to_gene[chr_id][(max(0, exons_list[0].start - flank_size), exons_list[0].start - 1)] = gene.id
                flanks_not_overlapping +=1
            total_flanks +=1            

            ovl = part_intervals[chr_id].overlaps(exons_list[-1].stop, exons_list[-1].stop + flank_size)
            if not ovl:
                flanks_to_gene[chr_id][(exons_list[-1].stop, exons_list[-1].stop + flank_size)] = gene.id
                flanks_not_overlapping +=1
            total_flanks +=1            

            for exon in exons_list:
                if exon.stop - exon.start < small_exon_threshold:
                    gene_to_small_exons[gene.id].append(exon.id)

    # for chr_id in flanks_to_gene:
    #     print(chr_id, flanks_to_gene[chr_id])
    print("total_flanks:", total_flanks)
    print("flanks_not_overlapping:", flanks_not_overlapping)

    # print(gene_to_small_exons)
    # for g in gene_to_small_exons:
    #     for e in gene_to_small_exons[g]:
    #         print(exon_id_to_choordinates[e])
    # sys.exit()
    return  exons_to_ref, parts_to_exons, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, exon_id_to_choordinates, \
            exon_to_gene, gene_to_small_exons, flanks_to_gene2, max_intron_chr



def get_part_sequences_from_choordinates(parts_to_exons, flanks_to_gene, refs):
    segments = {}
    tot_flanks = 0
    tot_parts = 0

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
                tot_parts += stop - start

            flank_instances = flanks_to_gene[chr_id]
            for flank in flank_instances:
                start,stop = flank[0], flank[1]
                seq = refs[chr_id][start : stop] 
                segments[chr_id][flank] = seq
                tot_flanks += stop - start

    print("Total parts size:", tot_parts)
    print("Total flanks size:", tot_flanks)

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


def get_flank_sequences_from_choordinates(flanks_to_gene, refs):
    flank_sequences = {}

    for chr_id in flanks_to_gene:
        if chr_id not in refs:
            continue
        else:
            flank_sequences[chr_id] = {}
            for flank in flanks_to_gene[chr_id]:
                start, stop = flank[0], flank[1]
                seq = refs[chr_id][start : stop] 
                flank_sequences[chr_id][flank] = seq

    return flank_sequences


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





