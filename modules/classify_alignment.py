
def contains(sub, pri):
    M, N = len(pri), len(sub)
    i, LAST = 0, M-N+1
    while True:
        try:
            found = pri.index(sub[0], i, LAST) # find first elem in sub
        except ValueError:
            return False
        if pri[found:found+N] == sub:
            return True
        else:
            i = found+1

def main(chr_id, predicted_exons, predicted_transcript, exons_to_transcripts, parts_to_transcript_annotations, transcripts_to_parts_annotations, all_parts_pairs_annotations, all_part_sites_annotations):

    # FSM
    transcripts = ''
    if tuple(predicted_transcript) in parts_to_transcript_annotations[chr_id]:
        transcript = ",".join( tr for tr in parts_to_transcript_annotations[chr_id][tuple(predicted_transcript)])
        print()
        print('Found, FSM to:', transcript)
        print()
        return "FSM", transcript
    elif tuple(predicted_exons) in exons_to_transcripts[chr_id]:
        transcript = ",".join( tr for tr in exons_to_transcripts[chr_id][tuple(predicted_exons)])  
        print()
        print('Found, FSM but not classified by parts to:', transcript)
        print()
        return "FSM", transcript

    # else:

    #     print('Did not find FSM', predicted_transcript)
    #     for ann_tr in parts_to_transcript_annotations[chr_id]:
    #         print(parts_to_transcript_annotations[chr_id][ann_tr] ,ann_tr)

    # ISM
    hits = [all_parts_pairs_annotations[chr_id][part_pair] for part_pair in all_parts_pairs_annotations[chr_id]]
    # print(hits)
    in_all_pairs = set.intersection(*hits)
    # print(in_all_pairs)
    # print(predicted_transcript)
    for transcript_id in in_all_pairs:
        transcript_parts = transcripts_to_parts_annotations[chr_id][transcript_id]
        if contains(predicted_transcript, transcript_parts):
            # print("Found, ISM to", transcript_id )
            transcript = transcript_id
            return "ISM", transcript
        else:
            print(predicted_transcript, transcript)

    # if set(predicted_transcript).issubset(all_parts_pairs_annotations[chr_id]): #change code here for true ISM they have to be in consecutive order!!!
    #     print()
    #     print('Found, ISM to:', tuple(predicted_transcript)  )
    #     print(parts_to_transcript_annotations[chr_id])
    # # else:
    # #     print('Did not find ISM', predicted_transcript)
    # print(predicted_transcript)
    # print(contains([(1,2),(3,4)], [1,(1,2),(3,4),5,6]))

    # else:
    all_sites_annotations_chr  = all_part_sites_annotations[chr_id] 
    is_nic = True
    for start, stop in predicted_transcript:
        if start not in all_sites_annotations_chr or stop not in all_sites_annotations_chr:
            is_nic = False
    if is_nic:
        all_pairs_annotations_chr = all_parts_pairs_annotations[chr_id]
        is_nic_comb = True
        for start, stop in predicted_transcript:
            if (start, stop) not in all_pairs_annotations_chr:
                is_nic_comb = False


        if is_nic_comb:
            print()
            print('Found, NIC (new combination of exons):', tuple(predicted_transcript) )
            print()             
            for ann_tr in parts_to_transcript_annotations[chr_id]:
                print(parts_to_transcript_annotations[chr_id][ann_tr] ,ann_tr)
            return  "NIC_comb", transcript

        else:
            print()
            print('Found, NIC (new donor-acceptor pair):', tuple(predicted_transcript) )
            print()   
            return   "NIC_novel", transcript          
    # else:
    #     print('Did not find NIC', predicted_transcript)

    # else:
    print()
    print('NNC:', tuple(predicted_transcript) )
    print()       
    return "NNC", transcript
