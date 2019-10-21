import os
import subprocess
from sys import stdout

from collections import defaultdict
from collections import namedtuple

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', "exon_part_id"])
globals()[mem.__name__] = mem # Global needed for multiprocessing

def find_mems(outfolder, refs_sequences, read_path, refs_path, mummer_out_path, min_mem):
    refs_path = open(refs_path, 'w') #open(os.path.join(outfolder, "refs_sequences_tmp.fa"), "w")
    for chr_id  in refs_sequences:
        for (start,stop), seq  in refs_sequences[chr_id].items():
            refs_path.write(">{0}\n{1}\n".format(chr_id + str("_") + str(start) + "_" + str(stop), seq))
    refs_path.close()
    # mummer_out_path = os.path.join( outfolder, "mummer_mems.txt" )
    with open(mummer_out_path, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open(os.path.join(outfolder, "mummer_errors.1") , "w")
        subprocess.check_call([ 'mummer',   '-maxmatch', '-l' , str(min_mem),  refs_path.name, read_path], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()
    # mummer_file = open(mummer_out_path, "r").readlines()
    # for line in mummer_file:
    #     print(line)
    # consensus = l[1].strip()
    # msa = [s.strip() for s in l[3:]]
    # print("regular spoa:", consensus)
    # print(len(consensus))
    # print(msa)
    # r = open(ref_out_file, "w")
    # r.write(">{0}\n{1}".format("reference", consensus))
    # r.close()
    # return consensus


def parse_results(mems_path):
    # file = open(os.path.join(mems_folder, "mummer_mems.txt"), 'r')
    mems_db = {}
    # read_mems_tmp = {}

    for i, line in enumerate(open(mems_path, 'r')):
        if line[0] == '>':
            if i == 0:
                read_acc = line.split()[1].strip()  # mems_db[line.split()[1].strip)()] = [] 
            else:
                mems_db[read_acc] = read_mems_tmp 
                read_acc = line.split()[1].strip() 
            
            read_mems_tmp = defaultdict(list)

        else:
            vals =  line.split() #11404_11606           1     11405       202
            exon_part_id = vals[0]
            chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('_')
            mem_len = int(vals[3])
            mem_ref_exon_part_start = int(vals[1])
            mem_read_start = int(vals[2])
            # convert to 0-indexed as python, however last coordinate is inclusive of the hit, not as in python end-indexing
            mem_tuple = mem(int(ref_coord_start) - 1 + mem_ref_exon_part_start - 1, int(ref_coord_start) - 1 + mem_ref_exon_part_start -1 + mem_len - 1,
                            mem_read_start-1, mem_read_start-1 + mem_len - 1, 
                            mem_len, exon_part_id)
            
            read_mems_tmp[chr_id].append( mem_tuple )
        # print(line)
    # add last read
    mems_db[read_acc] = read_mems_tmp 


    return mems_db
