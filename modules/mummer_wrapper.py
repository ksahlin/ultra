import os
import subprocess
from sys import stdout

def find_mems(refs_sequences, reads, outfolder):
    refs_path = open(os.path.join(outfolder, "refs_sequences_tmp.fa"), "w")
    for (start,stop), seq  in refs_sequences.items():
        refs_path.write(">{0}\n{1}\n".format(str(start) + "_" + str(stop), seq))
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