

#! /usr/bin/env python

import os
import sys
import os

import itertools
import argparse
import errno
import math


from collections import namedtuple

def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def traceback(index, C):

    return 

def main(mems):
    """
        Algorithm 15.1 in Genome scale algorithmic design, Makinen et al.

        Using lists instead of Binary search trees for range max queries,
        so n^2 time complexity instead of O(n log n).

        each mem is an Namedtuple. python object

    """
    mems = sorted(mems, key = lambda x: x.y )
    print(mems)
    T = [ (v.d, v.val)  for v in mems]
    I = [ (v.d, v.val)  for v in mems]
    
    # T_dict = {mems[i][3] : -10000 for i in range(len(mems))}
    # T_dict[0] = -10000
    # I_dict = {mems[i][3] : -10000 for i in range(len(mems))}
    # I_dict[0] = -10000

    C_a = [0]*(len(T))
    C_b = [0]*(len(T))
    C = [0]*(len(T))
    traceback_pointers = [None]*(len(T))

    for j in range(len(T)):
        v =  mems[j]

        # linear scan -- replace with range max Q tree
        T_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C) if  mems[j_prime].d < v.c]
        if T_values:
            # print(T_values)
            T_traceback_index, max_c_value_case_a = max(T_values, key=lambda x: x[1])
        else:
            max_c_value_case_a = 0
            T_traceback_index = None

        I_values = [(j_prime, c_val) for j_prime, c_val in enumerate(C) if v.c <= mems[j_prime].d  <= v.d]
        if I_values:
            # print(I_values)
            I_traceback_index, max_c_value_case_b = max(I_values, key=lambda x: x[1])
            I_v_prev_coord = mems[I_traceback_index].d
        else:
            I_value = 0
            I_v_prev_coord = None
            I_traceback_index = None


        C_a[j] = (v.d - v.c + 1) +  max_c_value_case_a
        C_b[j] = (v.d - I_v_prev_coord) + max_c_value_case_b # shouldnt it be v.d - v_tmp.d

        if C_a[j] >= C_b[j]:
            traceback_pointers[j] = T_traceback_index
        else:
            traceback_pointers[j] = I_traceback_index

        C[j] = max(C_a[j], C_b[j])

    # print('Value vector:', C)
    # print(argmax(C))
    # print(traceback_pointers)
    solution_index = argmax(C)
    # print()
    # print(C[solution_index], C)
    value = C[solution_index]
    # print()
    solution = []
    while True:
        solution.append(mems[solution_index])
        solution_index = traceback_pointers[solution_index]
        print(solution_index)
        if solution_index is None:
            break


    # print("Solution:", solution)
    return solution[::-1], value
    # traceback(C, best_solution_index)


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

