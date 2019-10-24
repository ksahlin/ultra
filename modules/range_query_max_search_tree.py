
from time import time

import queue 
import random

from collections import namedtuple
# def construct_tree(tree, mems_index, leafs):  
    
#     # assign values to leaves of  
#     # the range query tree  
#     n = len(leafs)
#     q = queue.Queue(2*n)
#     lvl = {}
#     for i in range(n-1,-1,-1):  
#         print(i)
#         # tree[n + i - 1] = leafs[i]
#         mems_index[n + i - 1] = leafs[i]
#         tree[n + i - 1] = (None, None) #leaf
#         lvl[n + i - 1] = 0
#         q.put(n + i - 1)
#         print("n", n + i - 1)
#     print()
#     node_index = n - 2
#     print("starting index:", node_index)
#     while not q.empty():
#         n1 = q.get()
#         n2 = q.get()
#         lvl_n1, lvl_n2 = lvl[n1], lvl[n2]
#         if lvl_n1 == lvl_n2:
#             new_internal_node = node_index
#             lvl[new_internal_node] = lvl_n1 +1  # n1 and n2 are the same level by def here
#             q.put(new_internal_node)
#         else:
#             n1
#             lvl[new_internal_node] = max(lvl_n1,lvl_n2) + 1
#         # print(q.get())
#     sys.exit()
#     # Build up graph 
#     # assign values to remaining nodes  
#     # to compute maximum in a given range  
#     for i in range(n - 2, -1, -1):
#         print(i, lvl)
#         mems_index[i] = min([mems_index[2 * i + 1], mems_index[2 * i + 2]], key = lambda x: x.d )  

#         if lvl[2 * i + 1] == lvl[2 * i + 2]:
#             tree[i] = (2 * i + 1, 2 * i + 2 )
#         elif lvl[2 * i + 1] > lvl[2 * i + 2]:
#             tree[i] = (2 * i + 2, 2 * i + 1 )
#         else:
#             print('BUG, should not happend because we fill nodes from right to left!')
#             sys.exit()

#         lvl[i] = max(lvl[2 * i + 1], lvl[2 * i + 2]) + 1
                          
# def range_query(tree, c, n): 
#     """ Basically the left and right indices  
#         will move towards right and left respectively  
#         and with every each next higher level and  
#         compute the minimum at each height change  
#         the index to leaf node first """
#     ma = 0 # initialize maximum value to 0 
#     curr_node = 0
#     while curr_node < n: # not yet reached a leaf
#         right_child_index = curr_node +1
#         left_child_index = tree[curr_node].d # always hods the left most value in subtree
#         if tree[right_child_index].d > c: # 

#         elif tree[left_child_index].d <= c:

#     while (left < right): 
#         if (left & 1): # if left index in odd  
#                 ma = max(ma, tree[left]) 
#                 left = left + 1
#         if (right & 1): # if right index in odd  
#                 right -= 1
#                 ma = max(ma, tree[right]) 
                  
#         # move to the next higher level 
#         left = left // 2
#         right = right // 2
#     return ma, left, right 
 


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
        tree[i] = min([tree[2 * i ], tree[2 * i + 1]], key = lambda x: x.d ) 
                          
def range_query(tree, search_coord, n): 
    ma = 0 # initialize minimum to 0 
    curr_node = 1
    while True: # not yet reached a leaf
        left_child_position = 2*curr_node
        right_child_position = 2*curr_node + 1
#         left_child_index = tree[curr_node].d # always holds the left most value in subtree
        if c > tree[right_child_position].d : # 
            curr_node = right_child_position         
        else: # tree[left_child_index].d <= c:
            curr_node = left_child_position         

        if curr_node >= n:
            break

    print("search coord: {0}, curr_node:{1}, leaf values:{2}.".format(search_coord, curr_node, tree[curr_node]))
    return tree[curr_node].Cj, tree[curr_node].j, curr_node

  

def update(tree, pos, value, n): 
      
    # change the index to leaf node first  
    pos += n  
      
    # update the value at the leaf node  
    # at the exact index 
    tree[pos] = value  
    while (pos > 1): 
          
        # move up one level at a time in the tree  
        pos >>= 1;  
          
        # update the values in the nodes  
        # in the next higher level 
        tree[pos] = max(tree[2 * pos],  
                           tree[2 * pos + 1])  
  





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


node = namedtuple('node', ['d', 'j', 'Cj'])
nodes = []
nodes.append( node(0, -1, 0) ) # add an start node in case a search is smaller than any d coord in tree
for i, mem in enumerate(mems):
    # if i > 3: break
    m = node(mem.d, mem.j, 0)
    nodes.append(m)

for i in range(20):
    if len(nodes) == 2**i or len(nodes) == 2**(i+1):
        break
    elif 2**i < len(nodes) < 2**(i+1):
        remainder = 2**(i+1) - len(nodes) 
        for i in range(remainder):
            nodes.append( node(0, -i - 2, 0) ) # fill up nodes to have leaves a power of 2
        break

leafs = sorted(nodes, key= lambda x: x.d)
n = len(leafs)
print(len(leafs))
print(leafs)
# T = {} 
T = [0 for i in range(2 * n) ]  
print(len(T))
# I = [0 for i in range(2 * n)]  
C = [0]*(len(leafs))
construct_tree(T, leafs, n)
print(len(T),T) 
# construct_segment_tree(I, leafs_I, n); 
print("time init RQmax:", time()- st)
st = time()
# print(nodes)
print("remainder", remainder)
trace_vector = [None]*(len(mems) + 1)
for j, mem in enumerate(mems):
    print(mem)
    c = mem.c
    # d = mem.d
    C_a_max, traceback_index, node_pos  = range_query(T, c, len(leafs)) 
    # j = node_pos - remainder - n
    print(j, C_a_max, traceback_index, node_pos)
    if traceback_index < 0:
        trace_vector[j+1] = 0
    else:
        trace_vector[j+1] = traceback_index

    # print(C_a_max, traceback_index )
    C_a =  C_a_max +  mem.d - mem.c   # add the mem_length to T since disjoint
    # update(T, node_pos, C_a, n); # point update 

    # trace_vector[]
    # # print(c, d, j_T,  val_T, C_a, l, r)
    # C[j_T] = C_a
    # trace_vector[j_T] = l

print(C)
print(trace_vector)
print("time querying RQ method 2:", time()- st)  
