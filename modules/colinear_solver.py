

#! /usr/bin/env python

import os
import sys
import os
import itertools
import argparse
import errno
import math
import copy


from collections import namedtuple


from modules import range_query_max_search_tree as RMaxQST


def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def max_both(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])

def traceback(index, C):
    return 

def all_solutions_c_max_indicies(C, C_max):
    return [i for i, c in enumerate(C) if c == C_max] 

# def all_solutions_c_max_indicies_mam(C, C_max):
#     return [i for i, c in enumerate(C) if c >= C_max - 1 ] 

# def reconstruct_solution(mems, C, trace_vector):
#     solution_index = argmax(C)
#     value = C[solution_index]
#     # print()
#     solution = []
#     while solution_index > 0:
#         solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
#         solution_index = trace_vector[solution_index]
#         # print(solution_index)
#         # if solution_index is None:
#         #     break
#     return value, solution[::-1]

def reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector, C, mam_mode = False):
    # solution_index = argmax(C)
    solutions = []
    for solution_index in all_C_max_indicies:
        value = C[solution_index]
        solution = []
        while solution_index > 0:
            if mam_mode:
                solution.append(mems[solution_index])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
            else:
                solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
            solution_index = trace_vector[solution_index]
        solutions.append( solution[::-1] )
    return value, solutions

def make_leafs_power_of_2(mems):
    nodes = []
    nodes.append( RMaxQST.Node(-1, -1, -2**32, -1) ) # add an start node in case a search is smaller than any d coord in tree
    for i, mem in enumerate(mems):
        # if i > 3: break
        m = RMaxQST.Node(mem.d, mem.j, -2**32, mem.j)
        nodes.append(m)

    for i in range(20):
        if len(nodes) == 2**i or len(nodes) == 2**(i+1):
            break
        elif 2**i < len(nodes) < 2**(i+1):
            remainder = 2**(i+1) - len(nodes) 
            for i in range(remainder):
                nodes.append( RMaxQST.Node(-1, -i - 2, -2**32, -i - 2) ) # fill up nodes to have leaves a power of 2
            break

    leafs = sorted(nodes, key= lambda x: x.d)
    # n = len(leafs)
    return leafs

def n_logn_read_coverage(mems):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using Binary search trees for range max queries,
        so n log n time complexity. Each mem is an Namedtuple. python object

    """
    # assert mems == sorted(mems, key=lambda x: x.y)


    T_leafs = make_leafs_power_of_2(mems)
    I_leafs = make_leafs_power_of_2(mems)
    n = len(T_leafs)
    T = [0 for i in range(2 * n) ]  
    I = [0 for i in range(2 * n) ]  
    # T_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(T, T_leafs, n)
    # I_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(I, I_leafs, n)

    mem_to_leaf_index = {l.j : i for i,l in enumerate(T_leafs)}

    C = [0]* (len(mems) + 1) #(len(leafs))
    trace_vector = [None]*(len(mems) + 1)

    RMaxQST.update(T, 0, 0, n) # point update 
    RMaxQST.update(I, 0, 0, n) # point update 

    for j, mem in enumerate(mems):
        leaf_to_update = mem_to_leaf_index[j]
        c = mem.c
        T_max, j_prime_a, node_pos  = RMaxQST.range_query(T, -1, c-1, len(T_leafs)) 
        # print("C_a:",  T_max +  mem.d - mem.c + 1, j_prime_a, node_pos, leaf_to_update )
        # print("T TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
        C_a =  T_max +  mem.d - mem.c + 1  # add the mem_length to T since disjoint

        if T_max < 0:
            print("BUG", T_max)
            sys.exit()

        
        d = mem.d
        I_max, j_prime_b, node_pos  = RMaxQST.range_query(I, c, d, len(I_leafs))         
        # print("C_b:", I_max +  mem.d, I_max, j_prime_b, node_pos, leaf_to_update )
        # print( I_max, mem.d, mems[j_prime_b].d, mems[j_prime_b])
        # print("I TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(I) if type(zz) != int])
        C_b =  I_max +  mem.d #- mems[j_prime_b].d   # add the part of the mem that is not overlapping

        # if C_b < 0:
        #     print("BUG")
        #     sys.exit()

        index, value = max_both([C_a, C_b])
        C[j+1] = value
        if index == 0: # Updating with C_a
            j_prime = j_prime_a
        else: # Updating with C_b
            j_prime = j_prime_b


        if j_prime < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
            trace_vector[j+1] = 0
        elif value == 0: # first j (i.e. j=0) 
            trace_vector[j+1]= 0
        else:
            trace_vector[j+1] = j_prime +1

        RMaxQST.update(T, leaf_to_update, value, n) # point update 
        RMaxQST.update(I, leaf_to_update, value - mem.d, n) # point update 


    # C_max, solution = reconstruct_solution(mems, C, trace_vector)
    # print("C", C)
    # print(trace_vector)

    solution_index = argmax(C)
    C_max = C[solution_index]
    all_C_max_indicies = all_solutions_c_max_indicies(C, C_max)
    # print("number solutions with the same score:", all_solutions_c_max_indicies(C, C_max))
    C_max, solutions = reconstruct_all_solutions(mems, all_C_max_indicies, trace_vector, C)

    return solutions, C_max #, is_unique_solution(C)

def read_coverage(mems, max_intron):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using lists instead of Binary search trees for range max queries,
        so n^2 time complexity instead of O(n log n).

        each mem is an Namedtuple. python object

    """
    # assert mems == sorted(mems, key = lambda x: x.y )

    if len(mems) > 1000:
        print('MEM',len(mems))

    # print("Going in to mem chaining:", mems)
    T = [ (v.d, v.val)  for v in mems]
    I = [ (v.d, v.val)  for v in mems]
    
    # T_dict = {mems[i][3] : -10000 for i in range(len(mems))}
    # T_dict[0] = -10000
    # I_dict = {mems[i][3] : -10000 for i in range(len(mems))}
    # I_dict[0] = -10000

    # C_a = [0]*(len(T))
    # C_b = [0]*(len(T))
    C = [0]*(len(T)+1)
    traceback_vector = [None]*(len(T)+1)

    for j in range(len(T)):
        v =  mems[j]

        # linear scan -- replace with range max Q tree
        T_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C[1:]) if  mems[j_prime].d < v.c and j_prime < j and v.y - mems[j_prime].y < max_intron ]
        if T_values:
            # print(T_values)
            T_traceback_index, max_c_value_case_a = max(reversed(T_values), key=lambda x: x[1])
        else:
            max_c_value_case_a = 0
            T_traceback_index = -1

        I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C[1:]) if v.c <= mems[j_prime].d  <= v.d and j_prime < j and v.y - mems[j_prime].y < max_intron ]
        if I_values:
            # print(I_values)
            I_values_plus_chord_diff = [ (j_prime, c_val + (v.d - mems[j_prime].d)) for j_prime, c_val in I_values]
            I_traceback_index, max_c_value_case_b = max(reversed(I_values_plus_chord_diff), key=lambda x: x[1])
            # I_v_prev_coord = mems[I_traceback_index].d
            # C_b[j] = (v.d - I_v_prev_coord) + max_c_value_case_b # shouldnt it be v.d - v_tmp.d
            C_b = max_c_value_case_b # shouldnt it be v.d - v_tmp.d

        else:
            I_v_prev_coord = v.c - 1
            I_traceback_index = -1
            max_c_value_case_b = 0
            C_b = 0


        C_a = (v.d - v.c + 1) +  max_c_value_case_a

        index, value = max_both([C_a, C_b])
        C[j+1] = value
        if index == 0: # Updating with C_a
            j_prime = T_traceback_index
        else: # Updating with C_b
            j_prime = I_traceback_index

        # C[j+1] = max(C_a, C_b)

        if j_prime < 0: # first j (i.e. j=0) 
            traceback_vector[j+1]= 0
        else:
            traceback_vector[j+1]= j_prime + 1

        # if j_prime == 0: # first j (i.e. j=0) 
        #     traceback_vector[j+1]= 0
        # elif C_a >= C_b:
        #     traceback_vector[j+1] = T_traceback_index + 1
        # else:
        #     traceback_vector[j+1] = I_traceback_index + 1

        # print(v.c, v.d, v.d -v.c, C_a[j], C_b[j], v.d, I_values, T_values)


    # solution_index = argmax(C)
    # value = C[solution_index]
    # solution = []
    # while True:
    #     solution.append(mems[solution_index])
    #     solution_index = traceback_vector[solution_index]
    #     if solution_index is None:
    #         break


    solution_index = argmax(C)
    # print(C)
    # print(traceback_vector)
    C_max = C[solution_index]
    all_C_max_indicies = all_solutions_c_max_indicies(C, C_max)
    # print(all_C_max_indicies)
    # print("number solutions with the same score:", all_solutions_c_max_indicies(C, C_max))
    C_max, solutions = reconstruct_all_solutions(mems, all_C_max_indicies, traceback_vector, C)
    # solutions = []
    # print("MEM Solution:", solution[::-1])
    return solutions, C_max #, is_unique_solution(C)
    # traceback(C, best_solution_index)


def read_coverage_mam_score(mams, overlap_threshold = 20):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using lists instead of Binary search trees for range max queries,
        so n^2 time complexity instead of O(n log n).

        each mem is an Namedtuple. python object

    """

    # assert mams == sorted(mams, key=lambda x: x.y)
    # mams = sorted(mams, key = lambda x: x.y )

    # print("MAM INSTANCE", mams)
    # for mam in mams:
    #     print(mam.mam_id, mam.x, mam.y, mam.c, mam.d, '\t', mam.val, mam.min_segment_length)
    if len(mams) > 1000:
        print('MAM',len(mams))
    T = [ (v.d, v.val)  for v in mams]
    I = [ (v.d, v.val)  for v in mams]
    
    # T_dict = {mams[i][3] : -10000 for i in range(len(mams))}
    # T_dict[0] = -10000
    # I_dict = {mams[i][3] : -10000 for i in range(len(mams))}
    # I_dict[0] = -10000

    C_a = [0]*(len(T))
    C_b = [0]*(len(T))
    C = [0]*(len(T))
    traceback_pointers = [None]*(len(T))

    for j in range(len(T)):
        v =  mams[j]

        # linear scan -- replace with range max Q tree
        T_values = [(j_prime, c_val - 0.1* (v.c - mams[j_prime].d - 1) ) for j_prime, c_val in enumerate(C) if  mams[j_prime].d < v.c and j_prime < j]
        # T_values2 = [(j_prime, c_val - max(0, mams[j_prime].y - v.x)) for j_prime, c_val in enumerate(C) if  mams[j_prime].d < v.c and j_prime < j] # Is this proper symmetric variant subtracting overlapping genomic positions?
        if T_values:
            # print(j, T_values)
            T_traceback_index, max_c_value_case_a = max(reversed(T_values), key=lambda x: x[1])
        else:
            max_c_value_case_a = 0
            T_traceback_index = None

        # I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C) if v.c <= mams[j_prime].d <= v.d and j_prime < j]
        I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C) if v.c <= mams[j_prime].d <= v.c - 1 + overlap_threshold and j_prime < j]
        # I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C) if v.c <= mams[j_prime].d  <= v.c + (v.d - v.c)*(1.0 - v.val/(v.d - v.c)) ]  #Maybe bug here, check the objective!! need to be careful
        if I_values:
            # print("here", j, I_values)
            # I_values_plus_chord_diff = [ (j_prime, c_val + (v.val - (mams[j_prime].d - v.c + 1 ) - 0.01 if v.x != mams[j_prime].y else v.val - 0.01)) for j_prime, c_val in I_values]  # Penalize an extra 0.1 to prefer non-overlapping solutions in case of tie breakers
            I_values_plus_chord_diff = [ (j_prime, c_val + (v.val - (mams[j_prime].d - v.c + 1 ) - 0.1*(mams[j_prime].d - v.c + 1 )) ) for j_prime, c_val in I_values]  # Penalize an extra 0.1 to prefer non-overlapping solutions in case of tie breakers

            I_traceback_index, max_c_value_case_b = max(reversed(I_values_plus_chord_diff), key=lambda x: x[1])
            C_b[j] = max_c_value_case_b

            # I_traceback_index, max_c_value_case_b = max(I_values, key=lambda x: x[1])
            # I_v_prev_coord = mams[I_traceback_index].d
            # C_b[j] = max_c_value_case_b + (v.val - v.val*(v.d - I_v_prev_coord)) # (v.d - I_v_prev_coord) +

        else:
            I_v_prev_coord = None
            I_traceback_index = None
            C_b[j] = 0

        C_a[j] = max_c_value_case_a + v.val # + (v.d - v.c + 1)

        if C_a[j] >= C_b[j]:
            traceback_pointers[j] = T_traceback_index
        else:
            # print("here")
            traceback_pointers[j] = I_traceback_index

        C[j] =  max(C_a[j], C_b[j]) #C_a[j] 

    # print('Value vector Max approx matches:', C)
    # print(traceback_pointers)

    # print(argmax(C))
    # print(traceback_pointers)
    solution_index = argmax(C)
    # solution_index = len(C) - argmax(C[::-1]) -1
    # print()
    # help_functions.eprint(C)
    value = C[solution_index]
    # print("index best sol:", solution_index, argmax(C), len(C), C)
    # print(traceback_pointers)

    # all_C_max_indicies = all_solutions_c_max_indicies_mam(C, value)
    # print(C)
    # print([m.x for m in mams])
    # print(traceback_pointers)
    # print(all_C_max_indicies)
    # print("number solutions with the same score:", all_solutions_c_max_indicies_mam(C, value))
    # C_max, solutions = reconstruct_all_solutions(mams, all_C_max_indicies, traceback_pointers, C, mam_mode = True)
    # for s in solutions:
    #     print(C_max, s)

    solution = []
    while True:
        solution.append(mams[solution_index])
        solution_index = traceback_pointers[solution_index]
        # print(solution_index)
        if solution_index is None:
            break

    non_unique = [x for x in set(C) if C.count(x) > 1]
    unique = True
    if non_unique:
        unique = False
        # print(traceback_pointers)

    # print("Solution:", solution[::-1])
    return tuple(solution[::-1]), value, unique
    # traceback(C, best_solution_index)



def n_logn_read_coverage_mams(mams, overlap_threshold = 5):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using Binary search trees for range max queries,
        so n log n time complexity. Each mem is an Namedtuple. python object

    """
    # assert mams == sorted(mams, key=lambda x: x.y)
    # for mam in mams:
    #     print(mam.mam_id, mam.x, mam.y, mam.c, mam.d, '\t', mam.val, mam.min_segment_length)
    # overlap_threshold = 20
    T_leafs = make_leafs_power_of_2(mams)
    I_leafs = make_leafs_power_of_2(mams)
    n = len(T_leafs)
    T = [0 for i in range(2 * n) ]  
    I = [0 for i in range(2 * n) ]  
    # T_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(T, T_leafs, n)
    # I_leafs = copy.deepcopy(leafs)
    RMaxQST.construct_tree(I, I_leafs, n)

    mem_to_leaf_index = {l.j : i for i,l in enumerate(T_leafs)}

    C = [0]* (len(mams) + 1) #(len(leafs))
    trace_vector = [None]*(len(mams) + 1)

    # aux_score = [0]* (len(mams) + 1) #(len(leafs))

    RMaxQST.update(T, 0, 0, n) # point update 
    RMaxQST.update(I, 0, 0, n) # point update 
    for j, mam in enumerate(mams):
        # print()
        # print( 'starting:', mam.c, mam.d)
        leaf_to_update = mem_to_leaf_index[j]
        c = mam.c
        T_max, j_prime_a, node_pos  = RMaxQST.range_query(T, -1, c-1, len(T_leafs)) 
        # print("C_a:", T_max, T_max + mam.val , j_prime_a, node_pos, leaf_to_update )
        # print("T TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
        C_a =  T_max + mam.val #mam.d - mam.c + 1   # add the mam_score to T since disjoint

        if T_max < 0:
            print("BUG", T_max)
            sys.exit()

        
        d = mam.d
        I_max, j_prime_b, node_pos  = RMaxQST.range_query(I, c, c - 1 + overlap_threshold, len(I_leafs))         
        # print( I_max, mam.d, mams[j_prime_b].d, mams[j_prime_b])
        # print("I TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(I) if type(zz) != int])
        prev_end = max(mam.c - 1, mams[j_prime_b].d)
        prev_end_diff = mams[j_prime_b].d - (mam.c - 1)
        # (mams[j_prime].d - v.c + 1 ) - 0.0001 if v.x != mams[j_prime].y else v.val - 0.0001)
        ovl_penalty = mams[j_prime_b].d - (mam.c - 1) + 0.0001 if mam.x != mams[j_prime_b].y else 0.0001
        # prev_score = aux_score[j_prime_b]
        # print("prev_score", prev_score)
        # print("C_b:", I_max, I_max  + mam.val - ovl_penalty , j_prime_b, node_pos, leaf_to_update )
        C_b =  I_max  + mam.val - ovl_penalty # (mam.d - prev_end)*mam.val/(mam.d-mam.c + 1) - 0.001 #mam.d #- mams[j_prime_b].d   # add the part of the mam that is not overlapping


        index, value = max_both([C_a, C_b])
        # print("taking:", index, value)
        C[j+1] = value
        # aux_score[j+1] = value
        if index == 0: # Updating with C_a
            j_prime = j_prime_a
        else: # Updating with C_b
            j_prime = j_prime_b


        if j_prime < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
            trace_vector[j+1] = 0
        elif value == 0: # first j (i.e. j=0) 
            trace_vector[j+1]= 0
        else:
            trace_vector[j+1] = j_prime +1

        RMaxQST.update(T, leaf_to_update, value, n) # point update 
        RMaxQST.update(I, leaf_to_update, value, n) # point update 


    solution_index = argmax(C)

    value = C[solution_index]
    solution = []
    while solution_index > 0:
        solution.append(mams[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
        solution_index = trace_vector[solution_index]
    # solutions.append( solution[::-1] )

    non_unique = [x for x in set(C) if C.count(x) > 1]
    unique = True
    if non_unique:
        unique = False
    # sys.exit()
    return tuple(solution[::-1]), value, unique





if __name__ == '__main__':
    mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', 'exon_part_id'])

    mems = [mem(0, 21,1000, 1021, 22, 'id1'), mem(100, 130, 1030,1060, 31, 'id2'), mem(65, 90, 6050, 6075, 26, 'id3'), mem(75, 101, 6060, 6086, 27, 'id4')]
    main(mems)

    import random
    # random.randint(0,100)
    mem_lengths = [ random.randint(1, 10)  for i in range(100)]
    print(mem_lengths)
    positions = [random.randint(0,100) for i in range(100)]
    read_intervals = [(pos, pos + mem_lengths[i] - 1) for i, pos in enumerate(positions)]
    print(read_intervals)

    positions = [random.randint(0,100) for i in range(100)]
    ref_intervals = [(pos, pos + mem_lengths[i] - 1) for i, pos in enumerate(positions)]
    print(ref_intervals)

    mems = [mem(read_iv[0],read_iv[1], ref_intervals[i][0], ref_intervals[i][1], mem_lengths[i], i ) for i, read_iv in enumerate(read_intervals)]
    print(mems)
    # sys.exit()
    main(mems)

