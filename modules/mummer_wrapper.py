import os
import subprocess
from sys import stdout

from collections import defaultdict
from collections import namedtuple


def find_mems(refs_sequences, reads, outfolder):
    refs_path = open(os.path.join(outfolder, "refs_sequences_tmp.fa"), "w")
    for chr_id  in refs_sequences:
        for (start,stop), seq  in refs_sequences[chr_id].items():
            refs_path.write(">{0}\n{1}\n".format(chr_id + str("_") + str(start) + "_" + str(stop), seq))
    refs_path.close()
    mummer_out_file = os.path.join( outfolder, "mummer_mems.txt" )
    with open(mummer_out_file, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open(os.path.join(outfolder, "mummer_errors.1") , "w")
        subprocess.check_call([ 'mummer',   '-maxmatch', '-l' , '20',  refs_path.name, reads], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()
    mummer_file = open(mummer_out_file, "r").readlines()
    for line in mummer_file:
        print(line)
    # consensus = l[1].strip()
    # msa = [s.strip() for s in l[3:]]
    # print("regular spoa:", consensus)
    # print(len(consensus))
    # print(msa)
    # r = open(ref_out_file, "w")
    # r.write(">{0}\n{1}".format("reference", consensus))
    # r.close()
    # return consensus

def parse_results(mems_folder):
    mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val'])

    file = open(os.path.join(mems_folder, "mummer_mems.txt"), 'r')
    mems_db = {}
    # read_mems_tmp = {}
    for i, line in enumerate(file):
        if line[0] == '>':
            if i == 0:
                acc = line.split()[1].strip()  # mems_db[line.split()[1].strip)()] = [] 
            else:
                mems_db[acc] = read_mems_tmp 
                acc = line.split()[1].strip() 
            
            read_mems_tmp = defaultdict(list)

        else:
            vals =  line.split() #11404_11606           1     11405       202
            mem_len = int(vals[3])
            read_start = int(vals[1])
            ref_start = int(vals[2])
            chr_id = vals[0].split('_')[0]
            mem_tuple = mem(read_start, read_start + mem_len, ref_start, ref_start + mem_len,  mem_len)
            
            read_mems_tmp[chr_id].append( mem_tuple )
        # print(line)
    return mems_db
