
import sys
import re
import os
import errno
import itertools

import parasail
import edlib
import dill as pickle 
import collections.abc

# def check_reference_headers(refs):
#     modified = False
#     for header in list(refs.keys()):
#         if header.isdigit() or header == 'X' or header == 'Y':
#             chr_id = 'chr'+ header
#         elif header == 'MT':
#             chr_id = 'chrM'
#         else:
#             chr_id = header

#         # we have modified 
#         if chr_id != header:
#             modified = True
#             seq = refs[header]
#             del refs[header]
#             refs[chr_id] = seq
#     return modified


def update_nested(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_nested(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def remove_read_polyA_ends(seq, threshold_len, to_len):
    end_length_window = min(len(seq)//2, 100)
    seq_list = [ seq[:-end_length_window] ]

    for ch, g in itertools.groupby(seq[-end_length_window:]):
        h_len = sum(1 for x in g)
        # print(ch, h_len, g )
        if h_len > threshold_len and (ch == "A" or ch == "T"):
            seq_list.append(ch*to_len)
        else:
            seq_list.append(ch*h_len)

    seq_mod = "".join([s for s in seq_list])
    return seq_mod

def pickle_dump(args, data, filename):
    with open(os.path.join(args.outfolder,filename), 'wb') as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break



def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    cig_pos = 0
    for length in result[:-1]:
        cig_pos += len(length)
        type_ = cigar[cig_pos]
        cig_pos += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln]), cigar_tuples


def edlib_alignment(read_seq, ref_seq, aln_mode="NW"):
    result = edlib.align(read_seq, ref_seq, task="path", mode=aln_mode)
    cigar_string = result["cigar"]
    start, stop = result['locations'][0]
    read_alignment, ref_alignment, cigar_tuples = cigar_to_seq(cigar_string, read_seq, ref_seq[start: stop])
    read_alignment = "-"*start + read_alignment + "-"* (len(ref_seq)-stop - 1)
    ref_alignment = ref_seq[:start] + ref_alignment + ref_seq[stop:]
    # cigar_tuples.insert(0, (len(ref_seq[:start]),"D"))
    # cigar_tuples.append((len(ref_seq[stop:]),"D"))
    # print( len(read_alignment), len(ref_alignment))
    assert len(read_alignment) == len(ref_alignment)
    return read_alignment, ref_alignment, result['editDistance']

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def parasail_alignment(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 3, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!",len(s1), len(s2))
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        print("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1, s2)
    # print(result.score, len(s1), len(s2))
    # print(s1_alignment)
    # print(s2_alignment)
    # print(cigar_string)
    # sys.exit()
    # print(dir(result))
    # print(result.end_query, result.end_ref, result.len_query, result.len_ref, result.length, result.matches)
    # print()
    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score

    # # Rolling window of matching blocks
    # match_vector = [ 1 if n1 == n2 else 0 for n1, n2 in zip(s1_alignment, s2_alignment) ]    
    # match_window = deque(match_vector[:k]) # initialization
    # current_match_count = sum(match_window)
    # aligned_region = []
    # if current_match_count >= match_id:
    #     aligned_region.append(1)
    # else:
    #     aligned_region.append(0)


    # for new_m_state in match_vector[k:]:
    #     prev_m_state = match_window.popleft()
    #     current_match_count = current_match_count - prev_m_state + new_m_state 
    #     match_window.append(new_m_state)
        
    #     if current_match_count >= match_id:
    #         aligned_region.append(1)
    #     else:        
    #         aligned_region.append(0)

    # # print("".join([str(m) for m in aligned_region]))
    # # print("Aligned ratio (tot aligned/len(seq1):", sum(aligned_region)/float(len(s1)))
    # alignment_ratio = sum(aligned_region)/float(len(s1))
    # return (s1, s2, (s1_alignment, s2_alignment, alignment_ratio))

def ssw_alignment(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 3, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.ssw(s1, s2, opening_penalty, gap_ext, user_matrix)
    print(result, type(result), dir(result))
    print(dir(result))
    for attr, value in result.__dict__.items():
        print(attr, value)
    # print(result.ref_begin1, result.ref_end1, result.read_begin1, result.read_end1)
    # print()
    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score

def parasail_local(s1, s2, match_score = 2, mismatch_penalty = -2, opening_penalty = 3, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sw_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!",len(s1), len(s2))
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        print("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    s1_alignment, s2_alignment, cigar_tuples = cigar_to_seq(cigar_string, s1[result.cigar.beg_query:result.end_query], s2[result.cigar.beg_ref: result.end_ref])
    # print(result.traceback.ref)
    # print(result.traceback.comp)
    # print(result.traceback.query)
    # print(result.score, len(s1), len(s2))
    print("read",s1_alignment)
    print("Rref",s2_alignment)
    print(result.cigar.beg_query,result.end_query)
    print(result.cigar.beg_ref, result.end_ref)
    print(cigar_string)
    # print(result.cigar.seq)

    # sys.exit()
    # print(dir(result))  
    # for attr, value in result.__dict__.items():
    #     print(attr, value)
    # print(result.end_query, result.end_ref, result.len_query, result.len_ref, result.length, result.matches)
    # print()
    return s1_alignment, s2_alignment, cigar_string, cigar_tuples, result.score


def find_all_paths(graph, start, end):
    path  = []
    paths = []
    queue = [(start, end, path)]
    while queue:
        start, end, path = queue.pop()
        # print( 'PATH', path)

        path = path + [start]
        if start == end:
            paths.append(path)
        for node in set(graph[start]).difference(path):
            queue.append((node, end, path))
    return paths
