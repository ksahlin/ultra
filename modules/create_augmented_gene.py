
from collections import defaultdict
from array import array

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

def dd_array(): # top level function declaration needed for multiprocessing
    return defaultdict(lambda :array("L"))


def add_tiling(tiling_segment_id, p1, p2, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                 tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                            tiling_parts_to_segments, tiling_gene_to_small_segments):
    k = 0
    while p2 - p1 - k > min_segment_size:
        # tiling_segment_name = "segm_{0}_{1}_{2}".format(chr_id, p1+k, p1 + k + min_segment_size)
        tiling_segment_name = tiling_segment_id
        tiling_segment_id_to_choordinates[tiling_segment_name] = (p1+k, p1 + k + min_segment_size)
        tiling_segment_to_ref[tiling_segment_name] = chr_id
        tiling_parts_to_segments[chr_id][(active_start, active_stop)].append(tiling_segment_name)
        tiling_segment_to_gene[tiling_segment_name] =  active_gene_ids 
        # total_unique_segment_counter += min_segment_size
        k += min_segment_size
        for gene_id in active_gene_ids:
            tiling_gene_to_small_segments[gene_id].append(tiling_segment_name)
        tiling_segment_id += 1

    # tiling_segment_name = "segm_{0}_{1}_{2}".format(chr_id, p1 + k - min_segment_size, p2)
    tiling_segment_name = tiling_segment_id
    tiling_segment_id_to_choordinates[tiling_segment_name] = (p1 + k - min_segment_size, p2)
    tiling_segment_to_ref[tiling_segment_name] = chr_id
    # tiling_parts_to_segments[chr_id][(active_start, active_stop)].add(tiling_segment_name)
    tiling_parts_to_segments[chr_id][(active_start, active_stop)].append(tiling_segment_name)
    tiling_segment_to_gene[tiling_segment_name] =  active_gene_ids 
    # total_unique_segment_counter += p2 - (p1 + k - min_segment_size)
    for gene_id in active_gene_ids:
        tiling_gene_to_small_segments[gene_id].append(tiling_segment_name)
    tiling_segment_id += 1
    return tiling_segment_id # need to return cause immutable



def get_canonical_segments(part_to_canonical_pos, part_count_to_choord, part_to_active_gene, pos_to_exon_ids, exon_id_to_choordinates, small_segment_threshold, min_segment_size):
    # parts_to_exons, exon_to_gene, exon_id_to_choordinates, exons_to_ref,
    segment_id_to_choordinates = {}
    segment_to_gene = {} 
    segment_to_ref = {}
    # parts_to_segments = defaultdict(dd_set)
    # gene_to_small_segments = defaultdict(list)
    parts_to_segments = defaultdict(dd_array)
    gene_to_small_segments = defaultdict(lambda :array("L"))
    temp_already_added_choord = set()
    total_unique_segment_counter = 0
    total_segments_bad = 0
    bad = 0

    # for all the poorly fitting reads to annotated start stop exon sites we employ a tiling segment structure and 
    # align reads in a second pass
    # M2 = defaultdict(lambda: defaultdict(lambda :array("I")))
    tiling_segment_id_to_choordinates = {}
    tiling_segment_to_gene = {} 
    tiling_segment_to_ref = {}
    tiling_parts_to_segments = defaultdict(dd_array)
    tiling_gene_to_small_segments = defaultdict(lambda :array("L"))

    # tiling_segment_id_to_choordinates = {}
    # tiling_segment_to_gene = {} 
    # tiling_segment_to_ref = {}
    # tiling_parts_to_segments = defaultdict(dd_set)
    # tiling_gene_to_small_segments = defaultdict(list)

    tiling_structures = [tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref, tiling_parts_to_segments, tiling_gene_to_small_segments]
    tiling_segment_id = 0
    segment_id = 0
    for (chr_id, part_id) in part_to_canonical_pos:
        active_start, active_stop = part_count_to_choord[(chr_id, part_id)]
        active_gene_ids = part_to_active_gene[(chr_id, part_id)]
        sorted_pos = sorted(part_to_canonical_pos[(chr_id, part_id)])

        open_starts_e_ids = set() #list( pos_to_exon_ids[(chr_id, part_id)][p, is_start] for p, is_start in pos_to_exon_ids[(chr_id, part_id)] if p < ) # exons spanning over the small segment
        # if chr_id == "SIRV3":
        #     print(sorted_pos)
        pos_tuples = [(p1, p2) for p1, p2 in zip(sorted_pos[:-1], sorted_pos[1:])]
        for i, (p1, p2) in enumerate(pos_tuples):
            open_starts_e_ids.update(pos_to_exon_ids[(chr_id, part_id)][p1, True]) # add the exons that start at this point 
            open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p1, False]) # remove the ones that ended here
            tiling_segment_id = add_tiling(tiling_segment_id, p1, p2, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                        tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                        tiling_parts_to_segments, tiling_gene_to_small_segments)
            if p2 - p1 >= min_segment_size:
                # print("here good", p2 - p1)
                # segment_name = "segm_{0}_{1}_{2}".format(chr_id,p1,p2)
                segment_name = segment_id
                temp_already_added_choord.add((chr_id, p1, p2))
                segment_id_to_choordinates[segment_name] = (p1, p2)
                segment_to_ref[segment_name] = chr_id
                parts_to_segments[chr_id][(active_start, active_stop)].append(segment_name)
                segment_to_gene[segment_name] =  active_gene_ids 
                total_unique_segment_counter += p2 - p1
                # print("Std:", segment_name, (chr_id, p1, p2))
                open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p2, False])
                # add small segments
                if p2 - p1 <= small_segment_threshold:
                    for gene_id in active_gene_ids:
                        gene_to_small_segments[gene_id].append(segment_name)
                segment_id += 1

            else:
                # add_earlier = False
                if i > 0: 
                    k = 1
                    while i - k >= 0:
                        if p2 - pos_tuples[i-k][0] >= min_segment_size: # or i - k == 0: # last case is i - k == 0 and we still have not added anything
                            # segment_name = "segm_{0}_{1}_{2}".format(chr_id, pos_tuples[i-k][0],p2)
                            # exon_name = "exon_{0}_{1}_{2}".format(chr_id, pos_tuples[i-k][0], p2) 
                            segment_name = segment_id

                            if (chr_id, pos_tuples[i-k][0],p2) not in temp_already_added_choord: #segment_id_to_choordinates and exon_name not in segment_id_to_choordinates:
                                # print("Added 1:", k, p2 - p1, p2 - pos_tuples[i-k][0])
                                total_segments_bad += p2 - pos_tuples[i-k][0]

                                segment_id_to_choordinates[segment_name] = (pos_tuples[i-k][0], p2)
                                segment_to_ref[segment_name] = chr_id
                                parts_to_segments[chr_id][(active_start, active_stop)].append(segment_name)
                                segment_to_gene[segment_name] =  active_gene_ids 
                                temp_already_added_choord.add( (chr_id, pos_tuples[i-k][0],p2) )
                                # print("while before:", segment_name,(chr_id, pos_tuples[i-k][0],p2))
                                segment_id += 1

                                # add_earlier = True
                                # add small segments
                                if p2 - pos_tuples[i-k][0] <= small_segment_threshold:
                                    for gene_id in active_gene_ids:
                                        gene_to_small_segments[gene_id].append(segment_name)
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
                            segment_name = segment_id

                            if (chr_id, p1, pos_tuples[i+k][1]) not in temp_already_added_choord: # segment_name not in segment_id_to_choordinates and exon_name not in segment_id_to_choordinates: 
                                # print("Added 2:", k, p2 - p1, pos_tuples[i+k][1] - p1)
                                total_segments_bad += pos_tuples[i+k][1] - p1

                                segment_id_to_choordinates[segment_name] = (p1, pos_tuples[i+k][1])
                                segment_to_ref[segment_name] = chr_id
                                parts_to_segments[chr_id][(active_start, active_stop)].append(segment_name)
                                segment_to_gene[segment_name] =  active_gene_ids 
                                temp_already_added_choord.add( (chr_id, p1, pos_tuples[i+k][1]) )
                                # print("while after:", segment_name, (chr_id, p1, pos_tuples[i+k][1]))

                                segment_id += 1

                                # add_later = True
                                # add small segments
                                if pos_tuples[i+k][1] - p1 <= small_segment_threshold:
                                    for gene_id in active_gene_ids:
                                        gene_to_small_segments[gene_id].append(segment_name)
                            break
                        else:
                            pass
                            # print("DID not make it 2:", k, pos_tuples[i+k][1] - p1, p1, p2)
                        k += 1

                if len(pos_tuples) == 1:
                    # segment_name = "segm_{0}_{1}_{2}".format(chr_id,p1,p2)
                    segment_name = segment_id

                    segment_id_to_choordinates[segment_name] = (p1, p2)
                    segment_to_ref[segment_name] = chr_id
                    parts_to_segments[chr_id][(active_start, active_stop)].append(segment_name)
                    segment_to_gene[segment_name] =  active_gene_ids 
                    total_unique_segment_counter += p2 - p1
                    temp_already_added_choord.add( (chr_id, p1, p2) )
                    # print("Third:", segment_name,  (chr_id, p1, p2))
                    segment_id += 1

                    # add small segments
                    if p2 - p1 <= small_segment_threshold:
                        for gene_id in active_gene_ids:
                            gene_to_small_segments[gene_id].append(segment_name)


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
                    tiling_segment_id = add_tiling(tiling_segment_id, e_start, e_stop, active_gene_ids, active_start, active_stop, min_segment_size, chr_id,
                                tiling_segment_id_to_choordinates, tiling_segment_to_gene, tiling_segment_to_ref,
                                tiling_parts_to_segments, tiling_gene_to_small_segments)

                    # exon_name = "exon_{0}_{1}_{2}".format(chr_id,e_start,e_stop) 
                    segment_name = segment_id

                    if (chr_id, e_start, e_stop) not in temp_already_added_choord: #"segm_{0}_{1}_{2}".format(chr_id,e_start,e_stop) not in segment_id_to_choordinates:                
                        segment_id_to_choordinates[segment_name] = (e_start, e_stop) 
                        segment_to_ref[segment_name] = chr_id
                        parts_to_segments[chr_id][(active_start, active_stop)].append(segment_name)
                        segment_to_gene[segment_name] = active_gene_ids
                        total_segments_bad += e_stop - e_start
                        temp_already_added_choord.add( (chr_id, e_start, e_stop) )
                        # print("Fourth:", segment_name,  (chr_id, e_start, e_stop))

                        segment_id += 1

                        if e_stop - e_start <= small_segment_threshold:
                            gene_to_small_segments[gene_id].append(segment_name)
                    else:
                        pass
                        # print("LOOOOOOL", e_stop - e_start)


                open_starts_e_ids.difference_update(pos_to_exon_ids[(chr_id, part_id)][p2, False])



                # relevant_starts = list(pos_to_exon_ids[(chr_id, part_id)][p1, True])
                # relevant_ends = list(pos_to_exon_ids[(chr_id, part_id)][p2, False])
                # all_segm_spanning = set(relevant_starts + relevant_ends)
                # # if chr_id == "SIRV3":
                # #     print("here bad", p1, p2)
                # #     print(relevant_starts)
                # #     print(relevant_ends)
                # if not all_segm_spanning:
                #     # print("BUG", open_starts_e_ids)
                #     # if chr_id == "SIRV5":
                #     # for e_id in open_starts_e_ids:
                #     #     print(exon_id_to_choordinates[e_id])
                #     # sys.exit()
                #     all_segm_spanning = open_starts_e_ids
                # # print(relevant_starts + relevant_ends)
                # bad += p2 - p1
                # for e_id in all_segm_spanning:
                #     e_start, e_stop = exon_id_to_choordinates[e_id]
                #     segment_id_to_choordinates[e_id] = (e_start, e_stop) 
                #     segment_to_ref[e_id] = chr_id
                #     parts_to_segments[chr_id][(active_start, active_stop)].add(e_id)
                #     segment_to_gene[e_id] = active_gene_ids
                #     total_segments_bad += e_stop - e_start
                #     # print("YOO", (e_start, e_stop))

                # # add small segments
                # if p2 - p1 <= small_segment_threshold:
                #     for gene_id in active_gene_ids:
                #         for e_id in all_segm_spanning:
                #             gene_to_small_segments[gene_id].append(e_id)

    print("total_unique_segment_counter", total_unique_segment_counter)
    print("total_segments_bad", total_segments_bad)
    print("bad", bad)
    # for e_id in segment_to_ref:
    #     if segment_to_ref[e_id] == "SIRV3":
    #         print(segment_id_to_choordinates[e_id])
    # sys.exit()
    # for s_id, (start,stop) in sorted(segment_id_to_choordinates.items(), key=lambda x: x[1]):
    #     if start > 34880000 and stop < 34980000:
    #         print(s_id, (start,stop))

    # for p_id, exons_ in sorted(parts_to_segments['chr18'].items(), key=lambda x: x[1]):
    #     if p_id[0] > 34860000 and p_id[1] < 34980000:
    #         print(p_id)
    # sys.exit()
    return parts_to_segments, segment_to_gene, segment_id_to_choordinates, segment_to_ref, gene_to_small_segments, tiling_structures 


def create_graph_from_exon_parts(db, flank_size, small_exon_threshold, min_segment_size): 
    """
        We need to link parts --> exons and exons --> transcripts
    """
    # print(dir(db))
    # genes_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exons_to_ref = {} # gene_id : { (exon_start, exon_stop) : set() }
    exon_to_gene = {} # exon_id : [gene_id ]

    flanks_to_gene2 = defaultdict(dict)  
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
    exon_choordinates_to_id = defaultdict(dd_set)

    for i, exon in enumerate(db.features_of_type('exon', order_by='seqid')):
        chr_id = exon.seqid
        exons_to_ref[exon.id] = chr_id
        # print(dir(exon))
        # print(exon.attributes["gene_id"])
        exon_gene_ids = exon.attributes["gene_id"] # is a list of strings
        exon_to_gene[exon.id] = exon_gene_ids
        exon_id_to_choordinates[exon.id] = (exon.start - 1, exon.stop)
        exon_choordinates_to_id[chr_id][(exon.start - 1, exon.stop)].add(exon.id)
        # creating the augmentation
        if i == 0: # initialization
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            active_exons = set() 
            active_exons.add(exon.id)    
            active_gene_ids = set(exon_gene_ids) 
            # adding very first flank on chromosome
            # print(max(0, exon.start - 1000), exon.start - 1)
            flanks_to_gene2[chr_id][(max(0, exon.start - 2*flank_size), exon.start - 1)] = "flank_{0}".format(i)
            total_flanks2 += 1
            total_flank_size += 2*flank_size
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(set(exon_gene_ids))

            continue
            # print(2*flank_size)

        if chr_id != prev_seq_id: # switching chromosomes
            parts_to_exons[prev_seq_id][(active_start, active_stop)] = active_exons
            part_count_to_choord[(prev_seq_id,part_counter)] = (active_start, active_stop)
            # part_intervals[prev_seq_id].addi(active_start, active_stop, None)
            # get_complementary_exon_seq_per_part(parts_to_exons, exon_to_gene, exon_id_to_choordinates, exons_to_ref, active_exons, active_start, active_stop, prev_seq_id)
            # adding very last flank on chromosome
            flanks_to_gene2[prev_seq_id][(max(0, active_stop), active_stop + 2*flank_size)] = "flank_{0}".format(i)
            total_flanks2 += 1
            total_flank_size += 2*flank_size
            # print(2*flank_size)

            # print(max(0, active_stop), active_stop + 1000)
            prev_seq_id = chr_id
            active_start = exon.start - 1
            active_stop = exon.stop
            part_counter +=1
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.start - 1)
            part_to_canonical_pos[(chr_id, part_counter)].add(exon.stop)
            part_to_active_gene[(chr_id, part_counter)].update(set(exon_gene_ids))
            pos_to_exon_ids[(chr_id, part_counter)][(exon.start - 1, True)].add(exon.id)
            pos_to_exon_ids[(chr_id, part_counter)][ (exon.stop, False)].add(exon.id) 

            active_exons = set() 
            active_exons.add(exon.id)     
            active_gene_ids = set(exon_gene_ids) 


        elif exon.start - 1 > active_stop:
            parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
            part_count_to_choord[(chr_id,part_counter)] = (active_start, active_stop)
            # get_complementary_exon_seq_per_part(parts_to_exons, exon_to_gene, exon_id_to_choordinates, exons_to_ref, active_exons, active_start, active_stop, prev_seq_id)
            # part_intervals[prev_seq_id].addi(active_start, active_stop, None)
            segment_size = 2*flank_size if set(exon_gene_ids).isdisjoint(active_gene_ids) else flank_size
            # part_intervalst(segment_size)
            if exon.start - active_stop > 2*segment_size:
                flanks_to_gene2[chr_id][(max(0, active_stop), active_stop + segment_size)] = "flank_{0}_1".format(i)
                flanks_to_gene2[chr_id][(max(0, exon.start - segment_size), exon.start - 1)] = "flank_{0}_2".format(i)
                total_flanks2 += 2
                total_flank_size += 2*segment_size
                # print(max(0, active_stop), active_stop + 1000)
                # print(max(0, exon.start - 1000), exon.start - 1)

            else: # add the whole intron
                flanks_to_gene2[chr_id][(max(0, active_stop), exon.start - 1)] = "flank_{0}".format(i)
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

    print("NR EXONS + COMPL:", len(exon_to_gene))
    print("total_flanks2:", total_flanks2)
    print("total_flank_size", total_flank_size)
    # print(flanks_to_gene2)
    parts_to_exons[chr_id][(active_start, active_stop)] = active_exons
    part_count_to_choord[(chr_id,part_counter)] = (active_start, active_stop)

    parts_to_segments, segment_to_gene, \
    segment_id_to_choordinates, segment_to_ref, \
    gene_to_small_segments, tiling_structures  = get_canonical_segments(part_to_canonical_pos, part_count_to_choord, part_to_active_gene, pos_to_exon_ids, exon_id_to_choordinates, small_exon_threshold, min_segment_size)

    print("total parts size:", sum( [stop - start for chrrr in parts_to_exons for start,stop in parts_to_exons[chrrr] ]))
    print("total exons size:", sum( [stop - start for start, stop in exon_id_to_choordinates.values() ]))
    # print( 'parts_to_segments', parts_to_segments["SIRV5"])
    # sys.exit()
    # parts_to_exons = parts_to_segments
    # exon_to_gene = segment_to_gene
    # exon_id_to_choordinates = segment_id_to_choordinates
    # exons_to_ref = segment_to_ref

    # print("parts:", [(start, stop) for chrrr in parts_to_exons for start,stop in parts_to_exons[chrrr] ])
    # sys.exit()

    # print(parts_to_exons)
    # part_intervals[prev_seq_id].addi(active_start, active_stop, None)
    min_intron = 2**32

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

            if i2 - i1 < min_intron:
                min_intron = i2 - i1
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
    print("min_intron:", min_intron)
    # sys.exit()
    # transcripts_to_splices = reverse_mapping(splices_to_transcripts)
    transcripts_to_splices = defaultdict(dict)
    for chr_id, chr_sp_sites_dict in splices_to_transcripts.items():
        for unique_sp_sites, tr_ids in chr_sp_sites_dict.items():
            for tr_id in tr_ids:
                transcripts_to_splices[chr_id][tr_id] = unique_sp_sites

    # gene_to_small_exons = {} # gene_id : [exon_id ]
    # # flanks_to_gene = defaultdict(dict)   
    # # flanks_not_overlapping = 0
    # # total_flanks = 0
    # for gene in db.features_of_type('gene', order_by='seqid'):
    #     gene_to_small_exons[gene.id] = []
    #     exons_list = [exon for exon in  db.children(gene, featuretype='exon', order_by='start')]
    #     # chr_id = gene.seqid
    #     if exons_list:
    #         # ovl = part_intervals[chr_id].overlaps(max(0, exons_list[0].start - flank_size), exons_list[0].start - 1)
    #         # if not ovl:
    #         #     flanks_to_gene[chr_id][(max(0, exons_list[0].start - flank_size), exons_list[0].start - 1)] = gene.id
    #         #     flanks_not_overlapping +=1
    #         # total_flanks +=1            

    #         # ovl = part_intervals[chr_id].overlaps(exons_list[-1].stop, exons_list[-1].stop + flank_size)
    #         # if not ovl:
    #         #     flanks_to_gene[chr_id][(exons_list[-1].stop, exons_list[-1].stop + flank_size)] = gene.id
    #         #     flanks_not_overlapping +=1
    #         # total_flanks +=1            

    #         for exon in exons_list:
    #             if exon.stop - exon.start < small_exon_threshold:
    #                 gene_to_small_exons[gene.id].append(exon.id)

    # # print(gene_to_small_segments["SIRV6"])
    # # sys.exit()
    # gene_to_small_exons = gene_to_small_segments

    return  segment_to_ref, parts_to_segments, splices_to_transcripts, \
            transcripts_to_splices, all_splice_pairs_annotations, \
            all_splice_sites_annotations, segment_id_to_choordinates, \
            segment_to_gene, gene_to_small_segments, flanks_to_gene2, max_intron_chr, \
            exon_choordinates_to_id, tiling_structures



def get_part_sequences_from_choordinates(parts_to_segments, flanks_to_gene, refs):
    part_segments = {}
    flank_segments = {}
    tot_flanks = 0
    tot_parts = 0

    for chr_id in parts_to_segments:
        if chr_id not in refs:
            continue
        else:
            parts_instance = parts_to_segments[chr_id]
            # chromosome = genes_to_ref[chr_id]
            part_segments[chr_id] = {}
            for part in parts_instance:
                start,stop = part[0], part[1]
                seq = refs[chr_id][start : stop] 
                part_segments[chr_id][part] = seq
                tot_parts += stop - start

            flank_instances = flanks_to_gene[chr_id]
            flank_segments[chr_id] = {}

            for flank in flank_instances:
                start,stop = flank[0], flank[1]
                seq = refs[chr_id][start : stop] 
                flank_segments[chr_id][flank] = seq
                tot_flanks += stop - start

    print("Total parts size:", tot_parts)
    print("Total flanks size:", tot_flanks)

    return part_segments, flank_segments


def get_segment_sequences_from_segment_id(segment_id_to_choordinates, segments_to_ref, refs):
    segment_sequences = defaultdict(dict)
    for segment_id in segment_id_to_choordinates:
        start,stop = segment_id_to_choordinates[segment_id]
        chr_id = segments_to_ref[segment_id]
        if chr_id not in refs:
            continue
        else:
            seq = refs[chr_id][start : stop] 
            segment_sequences[chr_id][(start,stop)] = seq
    # print(segments)
    return segment_sequences


def get_sequences_from_choordinates(flanks_to_gene, refs):
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





