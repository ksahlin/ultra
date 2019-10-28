
from time import time

import queue 
import random

from collections import namedtuple



from collections import namedtuple
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
        tree_node = Node(min_coord_node.d, min_coord_node.j, 0, min_coord_node.j_max)
        tree[i] = tree_node
                          
def range_query(tree, l, r, n): 
    # for right search coord
    pos = 1 # root position
    left_subtree_root_pos = set()
    while True: # not yet reached a leaf
        left_child_position = 2*pos
        right_child_position = 2*pos + 1
#         left_child_index = tree[pos].d # always holds the left most value in subtree
        if r > tree[right_child_position].d : # 
            if pos != 1: # We are not at root
                left_subtree_root_pos.add(left_child_position)        
            pos = right_child_position 

        else: # tree[left_child_index].d <= c:
            pos = left_child_position         

        if pos >= n:
            break

    left_subtree_root_pos.add(pos)

    print("search coord: {0}, node pos:{1}, leaf values:{2}.".format(r, pos, (tree[pos].d, tree[pos].Cj,tree[pos].j ) ))

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
            if pos != 1: # We are not at root    
                right_subtree_root_pos.add(right_child_position)        
            pos = left_child_position     

        if pos >= n:
            break
   
    right_subtree_root_pos.add(pos)

    # left_leaf_node_pos = pos
    print("right subtrees:", right_subtree_root_pos, "left subtrees:", left_subtree_root_pos)
    print("search coord: {0}, node pos:{1}, leaf values:{2}.".format(l, pos, (tree[pos].d, tree[pos].Cj,tree[pos].j ) ))

    # We now have the rightmost leaf node
    # now trace back the path back to the root to find maximum value
    # From chapter 3 in GSAD book.
    V_prime = left_subtree_root_pos -  right_subtree_root_pos
    # if V_prime:
    vl = max(V_prime, key = lambda x: tree[x].Cj)
    #     if l == L:
    #         vl = max([vl,L_pos], key = lambda x: tree[x].Cj)
    # elif l == L:
    #     vl = L_pos
    # else:
    #     print("BUGGGG")


    V_biss = right_subtree_root_pos - left_subtree_root_pos 
    # if V_biss:
    vr = max(V_biss, key = lambda x: tree[x].Cj)
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

    v_max_pos = max([vl,vr], key = lambda x: tree[x].Cj)
    C_max = tree[v_max_pos].Cj
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


st = time()

# construct sorted leafs

order_in_ref = [j  for j in range(4)] #[5,1,3,2,4,7,6,8]
choord_range = 10*max(order_in_ref)

mem_lengths = [10]*len(order_in_ref)

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val','j'])
mems = []
for (ref_index, mem_length) in zip(order_in_ref, mem_lengths):
    pos = 10*ref_index
    read_pos = random.randint(1,100)
    m = mem(pos, pos+mem_length,  read_pos, read_pos + mem_length, mem_length, ref_index)
    mems.append(m)


# node = namedtuple('node', ['d', 'j', 'Cj'])
nodes = []
nodes.append( Node(0, -1, 0, -1) ) # add an start node in case a search is smaller than any d coord in tree
for i, mem in enumerate(mems):
    # if i > 3: break
    m = Node(mem.d, mem.j, 0, mem.j)
    nodes.append(m)

for i in range(20):
    if len(nodes) == 2**i or len(nodes) == 2**(i+1):
        break
    elif 2**i < len(nodes) < 2**(i+1):
        remainder = 2**(i+1) - len(nodes) 
        for i in range(remainder):
            nodes.append( Node(0, -i - 2, 0, -i - 2) ) # fill up nodes to have leaves a power of 2
        break

leafs = sorted(nodes, key= lambda x: x.d)
n = len(leafs)
print(len(leafs))
print([l.j for l in  leafs] )
print([(m.y, m.c, m.d) for m in  mems])
# T = {} 
T = [0 for i in range(2 * n) ]  
print(len(T))
# I = [0 for i in range(2 * n)]  
construct_tree(T, leafs, n)
print(len(T),[ (t.j, t.d) for t in T if type(t) != int]) 
# construct_segment_tree(I, leafs_I, n); 
print("time init RQmax:", time()- st)
st = time()
# print(nodes)
print("remainder", remainder)
C = [0]* (len(mems) + 1) #(len(leafs))
trace_vector = [None]*(len(mems) + 1)

mem_to_leaf_index = {l.j : i for i,l in enumerate(leafs)}

for j, mem in enumerate(mems):
    print(mem)
    c = mem.c
    print("vals:", [l.Cj for l in leafs])
    C_a_max, traceback_index, node_pos  = range_query(T, 0, c, len(leafs)) 
    leaf_to_update = mem_to_leaf_index[j]
    print(C_a_max, traceback_index, node_pos, leaf_to_update )
    # assert node_pos == leaf_to_update
    # j = node_pos - remainder - n
    # print(j, C_a_max, traceback_index, node_pos, leaf_to_update)

    # print(C_a_max, traceback_index )
    C_a =  C_a_max +  mem.d - mem.c   # add the mem_length to T since disjoint
    update(T, leaf_to_update, C_a, n) # point update 
    C[j+1] = C_a
    if traceback_index < 0: # any of the additional leaf nodes (with negative index) we add to make number of leafs 2^n
        trace_vector[j+1] = 0
    elif C_a_max == 0: # first j (i.e. j=0) 
        trace_vector[j+1]= 0
    else:
        trace_vector[j+1] = traceback_index +1

    # trace_vector[]
    # # print(c, d, j_T,  val_T, C_a, l, r)
    # C[j_T] = C_a
    # trace_vector[j_T] = l

print(C)
print(trace_vector)
print("time querying RQ method 2:", time()- st)  
