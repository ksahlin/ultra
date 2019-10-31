
from time import time
import queue 
import random
import copy

from collections import namedtuple

import colinear_solver 


def construct_tree(tree, leafs, n):  
      
    # assign values to leaves of  
    # the segment tree  
    for i in range(n):  
        tree[n + i] = leafs[i] 
      
    # assign values to remaining nodes  
    # to compute maximum in a given range  
    for i in range(n - 1, 0, -1):  
        # print(i)
        # print(tree)
        min_coord_node = min([tree[2 * i ], tree[2 * i + 1]], key = lambda x: x.d ) # should always be left one though
        tree_node = Node(min_coord_node.d, min_coord_node.j, min_coord_node.Cj, min_coord_node.j_max)
        tree[i] = tree_node
                          
def range_query(tree, l, r, n): 
    assert l <= r
    # for right search coord
    pos = 1 # root position
    left_subtree_root_pos = set()
    while True: # not yet reached a leaf
        left_child_position = 2*pos
        right_child_position = 2*pos + 1
#         left_child_index = tree[pos].d # always holds the left most value in subtree
        if r >= tree[right_child_position].d : # 
            if pos != 1 and tree[left_child_position].d >= l: # We are not at root
                left_subtree_root_pos.add(left_child_position)        
            pos = right_child_position 

        else: # tree[left_child_index].d <= c:
            pos = left_child_position         

        if pos >= n:
            break

    R = tree[pos].d
    if R <= r:
        left_subtree_root_pos.add(pos)

    # print("search coord: {0}, node pos:{1}, leaf values:{2}.".format(r, pos, (tree[pos].d, tree[pos].Cj,tree[pos].j ) ))

    # for left search coord
    pos = 1 # root position
    right_subtree_root_pos = set()
    while True: # not yet reached a leaf
        left_child_position = 2*pos
        right_child_position = 2*pos + 1
#         left_child_index = tree[pos].d # always holds the left most value in subtree
        if l > tree[right_child_position].d : # 
            pos = right_child_position         
        else: # tree[left_child_index].d <= c:
            if pos != 1 and tree[right_child_position].d <= r: # We are not at root    
                right_subtree_root_pos.add(right_child_position)        
            pos = left_child_position     

        if pos >= n:
            break
   
    L = tree[pos].d
    if L >= l:
        right_subtree_root_pos.add(pos)

    # print("search coord: {0}, node pos:{1}, leaf values:{2}.".format(l, pos, (tree[pos].d, tree[pos].Cj,tree[pos].j ) ))

    # print("right subtrees:", right_subtree_root_pos, "left subtrees:", left_subtree_root_pos)

    # We now have the rightmost leaf node
    # now trace back the path back to the root to find maximum value
    # From chapter 3 in GSAD book.
    V_prime = left_subtree_root_pos #-  right_subtree_root_pos
    # V_prime.add(L)
    # if l == L:
    #     V_prime.add(L)
    #     # vl = max([vl,L_pos], key = lambda x: tree[x].Cj)

    if V_prime:
        vl = max(V_prime, key = lambda x: tree[x].Cj)
    else:
        vl = None
    #     if l == L:
    #         vl = max([vl,L_pos], key = lambda x: tree[x].Cj)
    # elif l == L:
    #     vl = L_pos
    # else:
    #     print("BUGGGG")


    V_biss = right_subtree_root_pos #- left_subtree_root_pos 
    # V_biss.add(R)
    # if r == R:
    #     V_biss.add(R)
    
    if V_biss:
        vr = max(V_biss, key = lambda x: tree[x].Cj)
    else:
        vr = None
       
    #     if r == R:
    #         vr = max([vr,R_pos], key = lambda x: tree[x].Cj)
    # elif r == R:
    #     vr = R_pos
    # else:
    #     print("BUGGGG")

    # if l == L:
    #     vl = max([vl,L_pos], key = lambda x: tree[x].Cj)
    # if r == R:
    #     vr = max([vr,R_pos], key = lambda x: tree[x].Cj)
    if vl is not None and vr is not None:
        v_max_pos = max([vl,vr], key = lambda x: tree[x].Cj)
        C_max = tree[v_max_pos].Cj
    elif vl is not None:
        v_max_pos = vl
        C_max = tree[vl].Cj
    elif vr is not None:
        v_max_pos = vr
        C_max = tree[vr].Cj
    else:
        print("BUG", l,r)
        print(tree)
        sys.exit()

    # trace down the j index of the leaf that produced the max here
    # while (pos > 1): 
          
    #     # move up one level at a time in the tree  
    #     pos >>= 1;  
    #     print('tracing max to: ', pos)
    #     # update the values in the nodes  
    #     # in the next higher level 
    #     C_max = max(tree[pos].Cj, C_max)  
    # print(C_max)
    return C_max, tree[v_max_pos].j_max, v_max_pos

  

def update(tree, leaf_pos, value, n): 
    # change the index to leaf node first  
    pos = leaf_pos + n
    # tree[pos].Cj = value
    # print('updating: ', pos)


      
    # update the value at the leaf node  
    # at the exact index 
    tree[pos].Cj = value  
    while (pos > 1): 
          
        # move up one level at a time in the tree  
        pos >>= 1  
        # print('updating: ', pos)
        # update the values in the nodes  
        # in the next higher level 
        cur_best = max([tree[2 * pos], tree[2 * pos + 1]], key=lambda x: x.Cj) 
        # tree[pos].Cj = max(tree[2 * pos].Cj,  
        #                    tree[2 * pos + 1].Cj)  
        tree[pos].Cj = cur_best.Cj
        tree[pos].j_max = cur_best.j_max


class Node:
    __slots__ = ['d', 'j', 'Cj', "j_max"]
    def __init__(self, d, j, Cj, j_max):
        self.d = d
        self.j = j
        self.Cj = Cj
        self.j_max = j_max


def argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

def max_both(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])


def reconstruct_solution(mems, C, trace_vector):
    solution_index = argmax(C)
    value = C[solution_index]
    # print()
    solution = []
    while solution_index > 0:
        solution.append(mems[solution_index - 1])  # trace vector is shifted on to the right so need to remove 1 from vectore to get j_index 
        solution_index = trace_vector[solution_index]
        # print(solution_index)
        # if solution_index is None:
        #     break
    return value, solution[::-1]



# construct sorted leafs

order_in_ref = [j  for j in range(10)] #[5,1,3,2,4,7,6,8]
choord_range = 10*max(order_in_ref)

mem_lengths = [10]*len(order_in_ref)

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val','j'])
mems = []
for (ref_index, mem_length) in zip(order_in_ref, mem_lengths):
    pos = 10*ref_index
    read_pos = random.randint(1,100)
    m = mem(pos, pos+mem_length,  read_pos, read_pos + mem_length, mem_length, ref_index)
    mems.append(m)


nodes = []
nodes.append( Node(0, -1, -2**32, -1) ) # add an start node in case a search is smaller than any d coord in tree
for i, mem in enumerate(mems):
    # if i > 3: break
    m = Node(mem.d, mem.j, -2**32, mem.j)
    nodes.append(m)

for i in range(20):
    if len(nodes) == 2**i or len(nodes) == 2**(i+1):
        break
    elif 2**i < len(nodes) < 2**(i+1):
        remainder = 2**(i+1) - len(nodes) 
        for i in range(remainder):
            nodes.append( Node(0, -i - 2, -2**32, -i - 2) ) # fill up nodes to have leaves a power of 2
        break

leafs = sorted(copy.deepcopy(nodes), key= lambda x: x.d)
n = len(leafs)
# print(len(leafs))
# print([l.j for l in  leafs] )
print([(m.y, m.c, m.d) for m in  mems])



# print(len(T),[ (t.j, t.d) for t in T if type(t) != int]) 
# construct_segment_tree(I, leafs_I, n); 

mem_to_leaf_index = {l.j : i for i,l in enumerate(leafs)}


############  T only stable #############
#########################################

st = time()
T = [0 for i in range(2 * n) ]  
construct_tree(T, leafs, n)

C = [0]* (len(mems) + 1) #(len(leafs))
trace_vector = [None]*(len(mems) + 1)

# print(mem_to_leaf_index)
update(T, 0, 0, n) # point update 
# sys.exit()
for j, mem in enumerate(mems):
    # print(mem)
    # print("vals:", [l.Cj for l in leafs])
    
    c = mem.c
    C_a_max, j_prime_a, node_pos  = range_query(T, 0, c - 1, len(leafs)) 
    leaf_to_update = mem_to_leaf_index[j]
    # print("TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
    # print("C_a:", C_a_max, j_prime_a, node_pos, leaf_to_update )
    if C_a_max < 0:
        print("BUG")
        sys.exit()
    C_a =  C_a_max +  mem.d - mem.c + 1  # add the mem_length to T since disjoint
    update(T, leaf_to_update, C_a, n) # point update 

    C[j+1] = C_a
    if j_prime_a < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
        trace_vector[j+1] = 0
    elif C_a_max == 0: # first j (i.e. j=0) 
        trace_vector[j+1]= 0
    else:
        trace_vector[j+1] = j_prime_a +1

C_max, solution = reconstruct_solution(mems, C, trace_vector)

print(C)
print(trace_vector)
print(C_max, [mem.j for mem in solution]) #, solution)
print("Total time RQmax T:", time()- st)
print()
print()
########################################
########################################


st = time()

solution, mem_solution_value, unique = colinear_solver.read_coverage(mems)
print(mem_solution_value, [mem.j for mem in solution])
print("Total time quadratic method:", time()- st)
print()
print()

############  T and I ###################
#########################################
st = time()
T = [0 for i in range(2 * n) ]  
I = [0 for i in range(2 * n) ]  
leafs = sorted(copy.deepcopy(nodes), key= lambda x: x.d)
construct_tree(T, leafs, n)
I_leafs = copy.deepcopy(leafs)
# for leaf in I_leafs:
#     leaf.Cj -= leaf.d 
# print([l.Cj for l in I_leafs])
construct_tree(I, I_leafs, n)

C = [0]* (len(mems) + 1) #(len(leafs))
trace_vector = [None]*(len(mems) + 1)

update(T, 0, 0, n) # point update 
update(I, 0, 0, n) # point update 

for j, mem in enumerate(mems):
    # print(mem)
    # print("vals T:", [l.Cj for l in leafs])
    # print("vals I:", [l.Cj for l in I_leafs])
    leaf_to_update = mem_to_leaf_index[j]

    c = mem.c
    T_max, j_prime_a, node_pos  = range_query(T, 0, c-1, len(leafs)) 
    # print("C_a:",  T_max +  mem.d - mem.c + 1, j_prime_a, node_pos, leaf_to_update )
    # print("T TREE:", [(s, zz.j, zz.d, zz.Cj, zz.j_max) for s, zz in enumerate(T) if type(zz) != int])
    C_a =  T_max +  mem.d - mem.c + 1  # add the mem_length to T since disjoint

    if T_max < 0:
        print("BUG", T_max)
        sys.exit()

    
    d = mem.d
    I_max, j_prime_b, node_pos  = range_query(I, c, d, len(I_leafs))         
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

    update(T, leaf_to_update, value, n) # point update 
    update(I, leaf_to_update, value - mem.d, n) # point update 

# print(trace_vector)

C_max, solution = reconstruct_solution(mems, C, trace_vector)
print(C)
print(trace_vector)
print(C_max , [mem.j for mem in solution])
print("Total time RQmax I and T:", time()- st)


# print("time querying RQ method 2:", time()- st)  




