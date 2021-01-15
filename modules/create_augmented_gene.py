
from collections import defaultdict
from array import array
from struct import *
# import intervaltree

from modules.help_functions import readfq

def reverse_mapping(d):
    d_dict = defaultdict(list)
    for k,v in d.items():
        for i in v:
            d_dict[i].append(k)
    return d_dict

def dd_set(): # top level function declaration needed for multiprocessing
    return defaultdict(set)

def dd_tuple(): # top level function declaration needed for multiprocessing
    return defaultdict(tuple)

def add_items(array, chr_id, p1, p2):
    array.append(chr_id)
    array.append(p1)
    array.append(p2)

def add_tiling(p1, p2, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                 tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                            tiling_parts_to_segments, tiling_gene_to_small_segments):
    part_name = array("L", [chr_id, active_start, active_stop]).tobytes()
    k = 0
    while p2 - p1 - k > min_segment_size:
        # tiling_segment_name = "segm_{0}_{1}_{2}".format(chr_id, p1+k, p1 + k + min_segment_size)
        # tiling_segment_name = tiling_segment_id
        tiling_segment_key = array("L", [chr_id, p1+k, p1 + k + min_segment_size]).tobytes()
        tiling_segment_id_to_choordinates[tiling_segment_key] = (p1+k, p1 + k + min_segment_size)
        tiling_segment_to_ref[tiling_segment_key] = chr_id
        # tiling_parts_to_segments[part_name].append(tiling_segment_key)
        add_items(tiling_parts_to_segments[part_name], chr_id, p1+k, p1 + k + min_segment_size)
        tiling_segment_to_gene[tiling_segment_key] =  active_gene_ids 
        # total_unique_segment_counter += min_segment_size
        for gene_id in active_gene_ids:
            # tiling_gene_to_small_segments[gene_id].append(tiling_segment_key)
            add_items(tiling_gene_to_small_segments[gene_id], chr_id, p1+k, p1 + k + min_segment_size)
        k += min_segment_size

    # tiling_segment_name = "segm_{0}_{1}_{2}".format(chr_id, p1 + k - min_segment_size, p2)
    # tiling_segment_name = tiling_segment_id
    tiling_segment_key = array("L", [chr_id, max(0, p1 + k - min_segment_size), p2]).tobytes()
    tiling_segment_id_to_choordinates[tiling_segment_key] = (max(0, p1 + k - min_segment_size), p2)
    tiling_segment_to_ref[tiling_segment_key] = chr_id
    # tiling_parts_to_segments[part_name].append(tiling_segment_key)
    add_items(tiling_parts_to_segments[part_name], chr_id, max(0, p1 + k - min_segment_size), p2)
    tiling_segment_to_gene[tiling_segment_key] =  active_gene_ids 
    # total_unique_segment_counter += p2 - (p1 + k - min_segment_size)
    for gene_id in active_gene_ids:
        # tiling_gene_to_small_segments[gene_id].append(tiling_segment_key)
        add_items(tiling_gene_to_small_segments[gene_id], chr_id, max(0, p1 + k - min_segment_size), p2)

    # tiling_segment_id += 1
    # return tiling_segment_id # need to return cause immutable



def get_canonical_segments(part_to_canonical_pos, part_count_to_choord, part_to_active_gene, pos_to_exon_ids, exon_id_to_choordinates, small_segment_threshold, min_segment_size):
    # parts_to_exons, exon_to_gene, exon_id_to_choordinates,
    segment_id_to_choordinates = {}
    segment_to_gene = {} 
    segment_to_ref = {}
    parts_to_segments = defaultdict(lambda :array("L")) #defaultdict(dd_array)
    gene_to_small_segments = defaultdict(lambda :array("L")) #defaultdict(lambda :array("b"))
    total_unique_segment_counter = 0
    total_segments_bad = 0
    bad = 0

    # for all the poorly fitting reads to annotated start stop exon sites we employ a tiling segment structure and 
    # align reads in a second pass
    # M2 = defaultdict(lambda: defaultdict(lambda :array("I")))
    tiling_segment_id_to_choordinates = {}
    tiling_segment_to_gene = {} 
    tiling_segment_to_ref = {}
    tiling_parts_to_segments = defaultdict(lambda :array("L")) #defaultdict(dd_array)
    tiling_gene_to_small_segments = defaultdict(lambda :array("L")) #= defaultdict(lambda :array("b"))


    tiling_structures = [tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref, tiling_parts_to_segments, tiling_gene_to_small_segments]
    # tiling_segment_id = 0
    # segment_id = 0
    for (chr_id, part_id) in part_to_canonical_pos:
        active_start, active_stop = part_count_to_choord[(chr_id, part_id)]
        part_name = array("L", [chr_id, active_start, active_stop]).tobytes()
        active_gene_ids = part_to_active_gene[(chr_id, part_id)]
        sorted_pos = sorted(part_to_canonical_pos[(chr_id, part_id)])

        open_starts_e_ids = set() #list( pos_to_exon_ids[(chr_id, part_id)][p, is_start] for p, is_start in pos_to_exon_ids[(chr_id, part_id)] if p < ) # exons spanning over the small segment
        pos_tuples = [(p1, p2) for p1, p2 in zip(sorted_pos[:-1], sorted_pos[1:])]
        for i, (p1, p2) in enumerate(pos_tuples):
            open_starts_e_ids.update(pos_to_exon_ids[(chr_id, part_id)][p1, True]) # add the exons that start at this point 
            open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p1, False]) # remove the ones that ended here
            add_tiling(p1, p2, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                        tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                        tiling_parts_to_segments, tiling_gene_to_small_segments)
            if p2 - p1 >= min_segment_size:
                # print("here good", p2 - p1)
                # segment_name = "segm_{0}_{1}_{2}".format(chr_id,p1,p2)
                # segment_name = segment_id
                segment_name = array("L", [chr_id, p1, p2]).tobytes()
                # temp_already_added_choord.add((chr_id, p1, p2))
                segment_id_to_choordinates[segment_name] = (p1, p2)
                segment_to_ref[segment_name] = chr_id
                # parts_to_segments[part_name].append(segment_name)
                add_items(parts_to_segments[part_name], chr_id, p1, p2)
                segment_to_gene[segment_name] =  active_gene_ids 
                total_unique_segment_counter += p2 - p1
                # print("Std:", segment_name, (chr_id, p1, p2))
                open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p2, False])
                # add small segments
                if p2 - p1 <= small_segment_threshold:
                    for gene_id in active_gene_ids:
                        # gene_to_small_segments[gene_id].append(segment_name)
                        add_items(gene_to_small_segments[gene_id], chr_id, p1, p2)
                # segment_id += 1

            else:
                # add_earlier = False
                if i > 0: 
                    k = 1
                    while i - k >= 0:
                        if p2 - pos_tuples[i-k][0] >= min_segment_size: # or i - k == 0: # last case is i - k == 0 and we still have not added anything
                            # segment_name = "segm_{0}_{1}_{2}".format(chr_id, pos_tuples[i-k][0],p2)
                            # exon_name = "exon_{0}_{1}_{2}".format(chr_id, pos_tuples[i-k][0], p2) 
                            # segment_name = segment_id
                            segment_name = array("L", [chr_id, pos_tuples[i-k][0], p2]).tobytes()

                            if segment_name not in segment_id_to_choordinates: #(chr_id, pos_tuples[i-k][0],p2) not in temp_already_added_choord: #segment_id_to_choordinates and exon_name not in segment_id_to_choordinates:
                                # print("Added 1:", k, p2 - p1, p2 - pos_tuples[i-k][0])
                                total_segments_bad += p2 - pos_tuples[i-k][0]

                                segment_id_to_choordinates[segment_name] = (pos_tuples[i-k][0], p2)
                                segment_to_ref[segment_name] = chr_id
                                # parts_to_segments[part_name].append(segment_name)
                                add_items(parts_to_segments[part_name], chr_id, pos_tuples[i-k][0], p2)
                                segment_to_gene[segment_name] =  active_gene_ids 
                                # temp_already_added_choord.add( (chr_id, pos_tuples[i-k][0],p2) )
                                # print("while before:", segment_name,(chr_id, pos_tuples[i-k][0],p2))
                                # segment_id += 1

                                # add_earlier = True
                                # add small segments
                                if p2 - pos_tuples[i-k][0] <= small_segment_threshold:
                                    for gene_id in active_gene_ids:
                                        # gene_to_small_segments[gene_id].append(segment_name)
                                        add_items(gene_to_small_segments[gene_id], chr_id, pos_tuples[i-k][0], p2)
                            break
                        else:
                            pass
                            # print("DID not make it 1:", k, p2 - pos_tuples[i-k][0], p1, p2)

                        k += 1

                # add_later = False
                if i < len(pos_tuples) - 1:
                    k = 1
                    while i + k <= len(pos_tuples) - 1: 
                        if pos_tuples[i+k][1] - p1 >= min_segment_size: # or i + k == len(pos_tuples) - 1:  # last case is i - k == 0 and we still have not added anything
                            # segment_name = "segm_{0}_{1}_{2}".format(chr_id,p1,pos_tuples[i+k][1])
                            # exon_name = "exon_{0}_{1}_{2}".format(chr_id,p1,pos_tuples[i+k][1]) 
                            # segment_name = segment_id
                            segment_name = array("L", [chr_id, p1, pos_tuples[i+k][1]]).tobytes()

                            if (chr_id, p1, pos_tuples[i+k][1]) not in segment_id_to_choordinates: # not in temp_already_added_choord: # segment_name not in segment_id_to_choordinates and exon_name not in segment_id_to_choordinates: 
                                # print("Added 2:", k, p2 - p1, pos_tuples[i+k][1] - p1)
                                total_segments_bad += pos_tuples[i+k][1] - p1

                                segment_id_to_choordinates[segment_name] = (p1, pos_tuples[i+k][1])
                                segment_to_ref[segment_name] = chr_id
                                # parts_to_segments[part_name].append(segment_name)
                                add_items(parts_to_segments[part_name], chr_id, p1, pos_tuples[i+k][1])
                                segment_to_gene[segment_name] =  active_gene_ids 
                                # temp_already_added_choord.add( (chr_id, p1, pos_tuples[i+k][1]) )
                                # print("while after:", segment_name, (chr_id, p1, pos_tuples[i+k][1]))

                                # segment_id += 1

                                # add_later = True
                                # add small segments
                                if pos_tuples[i+k][1] - p1 <= small_segment_threshold:
                                    for gene_id in active_gene_ids:
                                        # gene_to_small_segments[gene_id].append(segment_name)
                                        add_items(gene_to_small_segments[gene_id], chr_id, p1, pos_tuples[i+k][1])
                            break
                        else:
                            pass
                            # print("DID not make it 2:", k, pos_tuples[i+k][1] - p1, p1, p2)
                        k += 1

                if len(pos_tuples) == 1:
                    # segment_name = "segm_{0}_{1}_{2}".format(chr_id,p1,p2)
                    # segment_name = segment_id
                    segment_name = array("L", [chr_id, p1, p2]).tobytes()

                    segment_id_to_choordinates[segment_name] = (p1, p2)
                    segment_to_ref[segment_name] = chr_id
                    # parts_to_segments[part_name].append(segment_name)
                    add_items(parts_to_segments[part_name], chr_id, p1, p2)
                    segment_to_gene[segment_name] =  active_gene_ids 
                    total_unique_segment_counter += p2 - p1
                    # temp_already_added_choord.add( (chr_id, p1, p2) )
                    # print("Third:", segment_name,  (chr_id, p1, p2))
                    # segment_id += 1

                    # add small segments
                    if p2 - p1 <= small_segment_threshold:
                        for gene_id in active_gene_ids:
                            # gene_to_small_segments[gene_id].append(segment_name)
                            add_items(gene_to_small_segments[gene_id], chr_id, p1, p2)


                # either part consists only of a small exon (e7 in fig in paper) or 
                # one of the adjacent segment cutpoints (immediately before and/or after) are too close
                # which also makes them not pass the minimum segment size
                # in this case simply add the exons
                # if not add_later or not add_earlier: 
                relevant_starts = list(pos_to_exon_ids[(chr_id, part_id)][p1, True])
                relevant_ends = list(pos_to_exon_ids[(chr_id, part_id)][p2, False])
                all_segm_spanning = set(relevant_starts + relevant_ends)
                if not all_segm_spanning:
                    all_segm_spanning = open_starts_e_ids
                bad += p2 - p1
                for e_id in all_segm_spanning:
                    e_start, e_stop = exon_id_to_choordinates[e_id]
                    add_tiling(e_start, e_stop, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                                tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                                tiling_parts_to_segments, tiling_gene_to_small_segments)

                    # exon_name = "exon_{0}_{1}_{2}".format(chr_id,e_start,e_stop) 
                    # segment_name = segment_id
                    segment_name = array("L", [chr_id, e_start, e_stop]).tobytes()

                    if segment_name not in segment_id_to_choordinates: #(chr_id, e_start, e_stop) not in temp_already_added_choord: #"segm_{0}_{1}_{2}".format(chr_id,e_start,e_stop) not in segment_id_to_choordinates:                
                        segment_id_to_choordinates[segment_name] = (e_start, e_stop) 
                        segment_to_ref[segment_name] = chr_id
                        # parts_to_segments[part_name].append(segment_name)
                        add_items(parts_to_segments[part_name], chr_id, e_start, e_stop)
                        segment_to_gene[segment_name] = active_gene_ids
                        total_segments_bad += e_stop - e_start
                        # temp_already_added_choord.add( (chr_id, e_start, e_stop) )
                        # print("Fourth:", segment_name,  (chr_id, e_start, e_stop))

                        # segment_id += 1

                        if e_stop - e_start <= small_segment_threshold:
                            # gene_to_small_segments[gene_id].append(segment_name)
                            add_items(gene_to_small_segments[gene_id], chr_id, e_start, e_stop)
                    else:
                        pass
                        # print("LOOOOOOL", e_stop - e_start)


                open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p2, False])


    print("total_unique_segment_counter", total_unique_segment_counter)
    print("total_segments_bad", total_segments_bad)
    print("bad", bad)
    return parts_to_segments, segment_to_gene, segment_id_to_choordinates, segment_to_ref, gene_to_small_segments, tiling_structures 


def add_to_chr_mapping(chr_name, chr_to_id, id_to_chr):
    hash_id = len(chr_to_id)
    if chr_name in chr_to_id:
        return chr_to_id[chr_name]
    else:
        hash_id += 1
        chr_to_id[chr_name] = hash_id
        id_to_chr[hash_id] = chr_name
        return hash_id

        
def create_graph_from_exon_parts(db, flank_size, small_exon_threshold, min_segment_size, refs_lengths): 
    """
        We need to link parts --> exons and exons --> transcripts
    """
    # print(dir(db))
    # print([c for c in db.features_of_type('chr', order_by='seqid')])
    # print( [c for c in db.all_features()])
    # sys.exit()
    chr_to_id = {}
    id_to_chr = {}
    # genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    # exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    # exon_to_gene = {} # exon_id : [gene_id ]

    flank_ids = set()  
    total_flanks2 = 0
    total_flank_size = 0

    exon_id_to_choordinates = {}
    splices_to_transcripts = defaultdict(dd_set)
    all_splice_pairs_annotations = defaultdict(dd_set)
    all_splice_sites_annotations = defaultdict(set)
    # annotated_transcripts = defaultdict(set)
    # part_intervals = defaultdict(intervaltree.IntervalTree)
    parts_to_exons = defaultdict(dd_set)
    max_intron_chr = defaultdict(int)
    pos_to_exon_ids = defaultdict(dd_set)
    part_counter = 0
    part_to_canonical_pos = defaultdict(set)
    part_count_to_choord = defaultdict(dd_tuple)
    part_to_active_gene = defaultdict(set)
    exon_ids = set() #defaultdict(dd_set)

    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        chr_name = exon.seqid
        # if chr_name in refs_lengths:
        #     chr_length = refs_lengths[chr_name]
        #     if exon.start - 1 < 0 or exon.stop > chr_length:
        #         print("Warning. Annotation on chromosome {0} outside chromosome length [0,{1}]. Annotation was: ({2}, {3})".format(chr_name,chr_length,exon.start - 1, exon.stop))
        # if chr_name == "211000022278436":
        #     print("chromosome length [0,{1}]. Annotation was: ({2}, {3})".format(chr_name,chr_length,exon.start - 1, exon.stop))
        chr_id = add_to_chr_mapping(chr_name, chr_to_id, id_to_chr)
        # exons_to_ref[exon.id] = chr_id
        # print(dir(exon))
        # print(exon.attributes["gene_id"])
        exon_gene_ids = exon.attributes["gene_id"] # is a list of strings
        # exon_to_gene[exon.id] = exon_gene_ids
        exon_id_to_choordinates[exon.id] = (exon.start - 1, exon.stop)
        exon_name = array("L", [chr_id, exon.start - 1, exon.stop]).tobytes()
        exon_ids.add(exon_name)
        # creating the augmentation
        if i == 0: # initialization
            prev_chr_id = chr_id
            prev_chr_name = chr_name
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)    
            active_gene_ids = set(exon_gene_ids) 
            # adding very first flank on chromosome
            # print(max(0, exon.start - 1000), exon.start - 1)
            # to avoid flank intervals of length 0
            if exon.start - 1 > 0:
                flank_name =  array("L", [chr_id, max(0, exon.start - 2*flank_size), exon.start - 1]).tobytes()
                flank_ids.add(flank_name)
                total_flanks2 += 1
                total_flank_size += 2*flank_size

            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(set(exon_gene_ids))

            continue
            # print(2*flank_size)

        if chr_id != prev_chr_id: # switching chromosomes
            parts_to_exons[prev_chr_id][(active_start, active_stop)] = active_exons
            part_count_to_choord[(prev_chr_id,part_counter)] = (active_start, active_stop)
           
            # adding the very last flank on previous chromosome
            # chr_length = refs_lengths[prev_chr_name]
            chr_length = refs_lengths.get(prev_chr_name, active_stop + 2*flank_size) # if key 'prev_chr_name' is not present in refs_lengths (i.e. in reference fasta), it means that we will not use the annotations to this chromosome anyway so value does not matter..
            # to avoid flank intervals of length 0
            if active_stop < chr_length:
                flank_name = array("L", [prev_chr_id, max(0, active_stop), min(chr_length, active_stop + 2*flank_size)]).tobytes()
                flank_ids.add(flank_name)
                total_flanks2 += 1
                total_flank_size += 2*flank_size
                # print(2*flank_size)

            # print(max(0, active_stop), active_stop + 1000)
            prev_chr_id = chr_id
            prev_chr_name = chr_name
            active_start = exon.start - 1
            active_stop = exon.stop

            # adding the very first flank on new chromosome
            # to avoid flank intervals of length 0
            if exon.start - 1 > 0:
                flank_name = array("L", [chr_id, max(0, exon.start - 2*flank_size), exon.start - 1]).tobytes()
                flank_ids.add(flank_name)
                total_flanks2 += 1
                total_flank_size += 2*flank_size

            part_counter +=1
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(set(exon_gene_ids))
            pos_to_exon_ids[(chr_id, part_counter)][(exon.start - 1, True)].add(exon.id)
            pos_to_exon_ids[(chr_id, part_counter)][ (exon.stop, False)].add(exon.id) 

            active_exons = set() 
            active_exons.add(exon.id)     
            active_gene_ids = set(exon_gene_ids) 


        elif exon.start - 1 > active_stop + 20:
            parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
            part_count_to_choord[(chr_id,part_counter)] = (active_start, active_stop)
            segment_size = 2*flank_size if set(exon_gene_ids).isdisjoint(active_gene_ids) else flank_size
            # part_intervalst(segment_size)
            if exon.start - active_stop > 2*segment_size:
                # chr_length = refs_lengths[chr_name]
                chr_length = refs_lengths.get(chr_name, active_stop + segment_size) # if key 'chr_name' is not present in refs_lengths (i.e. in reference fasta), it means that we will not use the annotations to this chromosome anyway so value does not matter..
                flank_name = array("L", [chr_id, max(0, active_stop), min(chr_length, active_stop + segment_size)]).tobytes()
                flank_ids.add(flank_name)
                flank_name = array("L", [chr_id, max(0, exon.start - segment_size), exon.start - 1]).tobytes()
                flank_ids.add(flank_name) 
                total_flanks2 += 2
                total_flank_size += 2*segment_size
                # print(max(0, active_stop), active_stop + 1000)
                # print(max(0, exon.start - 1000), exon.start - 1)

            else: # add the whole intron
                flank_name = array("L", [chr_id, max(0, active_stop), exon.start - 1]).tobytes()
                flank_ids.add(flank_name)
                total_flanks2 += 1
                total_flank_size += (exon.start - 1 - max(0, active_stop))

            active_exons = set()
            active_exons.add(exon.id)
            active_gene_ids = set(exon_gene_ids) 

            active_start = exon.start - 1
            active_stop = exon.stop
            part_counter +=1
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(set(exon_gene_ids))
            pos_to_exon_ids[(chr_id, part_counter)][(exon.start - 1, True)].add(exon.id)
            pos_to_exon_ids[(chr_id, part_counter)][ (exon.stop, False)].add(exon.id)

        else:
            active_exons.add(exon.id)    
            active_stop = max(active_stop, exon.stop)
            active_gene_ids.update(exon_gene_ids) 
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(exon_gene_ids)
            pos_to_exon_ids[(chr_id, part_counter)][(exon.start - 1, True)].add(exon.id)
            pos_to_exon_ids[(chr_id, part_counter)][ (exon.stop, False)].add(exon.id)

        assert active_start <= exon.start - 1

    # addig the very last flank at the last chromosome in the annotation
    # chr_length = refs_lengths[chr_name]
    chr_length = refs_lengths.get(chr_name, active_stop + 2*flank_size) # if key 'chr_name' is not present in refs_lengths (i.e. in reference fasta), it means that we will not use the annotations to this chromosome anyway so value does not matter..
    # to avoid flank intervals of length 0
    if active_stop < chr_length:
        flank_name = array("L", [chr_id, max(0, active_stop), min(chr_length, active_stop + 2*flank_size)]).tobytes()
        flank_ids.add(flank_name)
        total_flanks2 += 1
        total_flank_size += 2*flank_size

    # print("NR EXONS + COMPL:", len(exon_to_gene))
    print("total_flanks2:", total_flanks2)
    print("total_flank_size", total_flank_size)
    # print(flank_ids)
    parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
    part_count_to_choord[(chr_id,part_counter)] = (active_start, active_stop)

    parts_to_segments, segment_to_gene, \
    segment_id_to_choordinates, segment_to_ref, \
    gene_to_small_segments, tiling_structures  = get_canonical_segments(part_to_canonical_pos, part_count_to_choord, part_to_active_gene, pos_to_exon_ids, exon_id_to_choordinates, small_exon_threshold, min_segment_size)

    print("total parts size:", sum( [stop - start for chrrr in parts_to_exons for start,stop in parts_to_exons[chrrr] ]))
    print("total exons size:", sum( [stop - start for start, stop in exon_id_to_choordinates.values() ]))
    # print(chr_to_id)
    # print(id_to_chr)

    min_intron = 2**32

    for transcript in db.features_of_type('transcript', order_by='seqid'): #db.children(gene, featuretype='transcript', order_by='start'):
        # if transcript.seqid.isdigit() or transcript.seqid == 'X' or transcript.seqid == 'Y':
        #     chr_name = 'chr'+ transcript.seqid
        # elif transcript.seqid == 'MT':
        #     chr_name = 'chrM'
        # else:
        chr_name = transcript.seqid
        chr_id = chr_to_id[chr_name]
        # print("here", transcript.id,  str(chr_id))  
        transcript_exons = []
        for exon in db.children(transcript, featuretype='exon', order_by='start'):
            transcript_exons.append( (exon.start-1, exon.stop) )
        internal_transcript_splices = [ (e1[1],e2[0]) for e1, e2 in zip(transcript_exons[:-1],transcript_exons[1:])]
        
        for i1,i2 in internal_transcript_splices:
            if i2 - i1 > max_intron_chr[chr_id]:
                max_intron_chr[chr_id] = i2 - i1
                # print("intron size:", i2 - i1, chr_id)

            if i2 - i1 < min_intron:
                min_intron = i2 - i1
        # internal transcript splices
        splices_to_transcripts[chr_id][ tuple(internal_transcript_splices)].add(transcript.id)
        for site1, site2 in internal_transcript_splices:
            all_splice_pairs_annotations[chr_id][(site1, site2)].add( transcript.id )
            # if site2 == 3105:
            #     print('LOOOL',chr_id)
            #     ccc.add(chr_id)

            all_splice_sites_annotations[chr_id].add(site1)
            all_splice_sites_annotations[chr_id].add(site2)
        
        # add start and end splice to all_splice_sites_annotations 
        if transcript_exons:
            all_splice_sites_annotations[chr_id].add(transcript_exons[0][0])
            all_splice_sites_annotations[chr_id].add(transcript_exons[-1][-1])
        else:
            print("Something is wrong with transcript annotation: {0} on gene: {1}, and could not be added. Check that the gene ID and transcript ID is not the same!".format(transcript.id, transcript.seqid))
            # sys.exit()
    print("min_intron:", min_intron)
    # sys.exit()
    # transcripts_to_splices = reverse_mapping(splices_to_transcripts)
    transcripts_to_splices = defaultdict(dict)
    for chr_id, chr_sp_sites_dict in splices_to_transcripts.items():
        for unique_sp_sites, tr_ids in chr_sp_sites_dict.items():
            for tr_id in tr_ids:
                transcripts_to_splices[chr_id][tr_id] = unique_sp_sites

    # print(chr_to_id)
    # for f in flank_ids:
    #     chr_id, start, stop = unpack('LLL',f)
    #     if chr_id == 7:
    #         print(chr_id, start, stop)

    return  segment_to_ref, parts_to_segments, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, segment_id_to_choordinates, \
            segment_to_gene, gene_to_small_segments, flank_ids, max_intron_chr, \
            exon_ids, chr_to_id, id_to_chr, tiling_structures


def get_sequences_from_choordinates(sequence_choordinates, refs):
    sequence_container = defaultdict(dict)
    for sequence_id in sequence_choordinates:
        chr_id, start, stop = unpack('LLL',sequence_id)
        # start,stop = sequence_choordinates[sequence_id]
        # chr_id = segments_to_ref[sequence_id]
        if chr_id not in refs:
            continue
        else:
            seq = refs[chr_id][start : stop] 
            # key = array("L", [chr_id,start,stop]).tobytes()
            sequence_container[sequence_id] = seq
    # print(segments)
    return sequence_container


# Functions for masking overly abundant kmers 

import itertools
from collections import defaultdict, deque
def kmer_counter(ref_part_sequences, kmer_size):
    count = defaultdict(int)
    position_count = defaultdict(list)
    it_count = 0
    for sequence_id, seq  in ref_part_sequences.items():
        # for part, seq in ref_part_sequences[chr_id].items():
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
    for sequence_id, seq in ref_part_sequences.items():
        # chr_id, start, stop = unpack('LLL',sequence_id)
        # seq = ref_part_sequences[sequence_id]
        # for part, seq in ref_part_sequences[chr_id].items():
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
                ref_part_sequences[sequence_id] = seq_masked
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





