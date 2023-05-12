import os
import subprocess
from sys import stdout, exit

from collections import defaultdict
from collections import namedtuple

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', 'j', "exon_part_id"])
globals()[mem.__name__] = mem # Global needed for multiprocessing

def find_mems_mummer(outfolder, read_path, refs_path, mummer_out_path, min_mem):
    with open(mummer_out_path, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open(os.path.join(outfolder, "mummer_errors.1") , "w")
        subprocess.check_call([ 'mummer',   '-maxmatch', '-l' , str(min_mem),  refs_path, read_path], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()


def find_mems_slamem(outfolder, read_path, refs_path, out_path, min_mem):
    # time slaMEM -l 14 /Users/kxs624/tmp/ULTRA/human_test/refs_sequences.fa /Users/kxs624/tmp/ULTRA/human_test_new_flanking_strat/reads_tmp.fq -o /Users/kxs624/tmp/ULTRA/human_test/slamem_test.tx
    # with open(out_path, "w") as output_file:
    tmp = out_path.split("seeds_batch_")[1]
    batch_id =  tmp.split('.')[0]
    stdout.flush()
    stderr_file = open(os.path.join(outfolder, "slamem_stderr_{0}.1".format(batch_id)) , "w")
    stdout_file = open(os.path.join(outfolder, "slamem_stdout_{0}.1".format(batch_id)) , "w")
    try: # slaMEM throws error if no MEMs are found in any of the sequences
        subprocess.check_call([ 'slaMEM', '-l' , str(min_mem),  refs_path, read_path, '-o', out_path ], stdout=stdout_file, stderr=stderr_file)
        print("Using SLAMEM")
    except:
        find_mems_mummer(outfolder, read_path, refs_path, out_path, min_mem)
        print("Using MUMMER")

    # print('Done.')
    stdout.flush()
    # output_file.close()


def find_nams_strobemap(outfolder, read_path, refs_path, out_path, nr_cores, min_mem, strobe_size, thinning_level):
    # /usr/bin/time -l ./StrobeMap -n 2 -k 9 -w 30 -t 3 -s  -o /Users/kxs624/tmp/ULTRA/human_test/refs_sequences.fa /Users/kxs624/tmp/ULTRA/human_test_new_flanking_strat/reads_tmp.fq
    # StrobeMap -n 2 -k 9 -w 30 -t 3 -s  -o ~/tmp/STROBEMERS/multithreading/ /Users/kxs624/tmp/ULTRA/dros_tmp/refs_sequences.fa /Users/kxs624/tmp/ULTRA/dros_test/reads_16xrep.fa
    # with open(out_path, "w") as output_file:
    stdout.flush()
    stderr_file = open(os.path.join(outfolder, "namfinder_stderr.1") , "w")
    stdout_file = open(os.path.join(outfolder, "namfinder_stdout.1") , "w") # should be empty
    outfile = os.path.join(outfolder, "seeds.txt")
    if thinning_level == 0:
        s = strobe_size
        l = strobe_size + 1
        u = 35 # seems to be ok value based on some tests
    elif thinning_level == 1: # expected seed distance: 3
        s = strobe_size - 2
        l = (strobe_size + 1)//3 # Need to scale down offsets since seeds subsampled
        u = 35//3 # seems to be ok value based on some tests
    elif thinning_level == 2: # expected seed distance: 5
        s = strobe_size - 4 
        l = (strobe_size + 1)//5 # Need to scale down offsets since seeds subsampled
        u = 35//5 # seems to be ok value based on some tests

    try: # slaMEM throws error if no MEMs are found in any of the sequences
        print("Using namfinder")
        print([ 'namfinder', '-k' , str(strobe_size), '-s' , str(s), '-l' , str(l), '-u' , str(u), '-C' , '500', '-L' , '1000', '-t', str(nr_cores), '-S', '-o', outfile, refs_path, read_path ])
        subprocess.check_call([ 'namfinder', '-k' , str(strobe_size), '-s' , str(s), '-l' , str(l), '-u' , str(u), '-C' , '500', '-L' , '1000', '-t', str(nr_cores), '-S', '-o', outfile, refs_path, read_path ], stdout=stdout_file, stderr=stderr_file)
    except:
        print("An unexpected error happend in namfinder, check error log at:", stderr_file)
        print("If you beileive this is a bug in namfinder, report an issue at: https://github.com/ksahlin/namfinder")
        sys.exit()
    stdout.flush()
    return outfile

def get_mem_records(mems_path, reads):
    '''
        Reads contains all the relevant reads in the batch to read mems from 
    '''
    relevant = False
    relevant_read_cnt = 0
    for i, line in enumerate(open(mems_path, 'r')):
        if line[0] == '>':
            acc = line[1:].strip()
            if acc not in reads:
                relevant = False
                continue
            else:
                relevant = True
                relevant_read_cnt +=1

            if relevant_read_cnt == 1:
                read_acc = acc  
            else:

                for chr_id in list(read_mems_tmp.keys()):
                    coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
                    sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
                    read_mems_tmp[chr_id] = sorted_mems

                yield read_acc, read_mems_tmp
                read_acc = acc 
            
            read_mems_tmp = defaultdict(list)

        elif relevant:
                vals =  line.split() #11404_11606           1     11405       202
                exon_part_id = vals[0]
                chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('^')
                chr_id = int(chr_id)
                mem_len = int(vals[3])
                # convert to 0-indexed reference as in python
                # however, for MEM length last coordinate is inclusive of the hit in MEM solvers, not as in python end-indexing
                mem_ref_exon_part_start = int(vals[1]) - 1
                mem_read_start = int(vals[2]) - 1
                ref_coord_start = int(ref_coord_start) # has already been 0-indexed when constructing parts
                mem_genome_start = ref_coord_start + mem_ref_exon_part_start
                                
                info_tuple = ( mem_genome_start, mem_genome_start + mem_len - 1,
                                mem_read_start, mem_read_start + mem_len - 1, 
                                mem_len, exon_part_id)
                read_mems_tmp[chr_id].append( info_tuple )


    for chr_id in list(read_mems_tmp.keys()):
        coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
        sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
        read_mems_tmp[chr_id] = sorted_mems
    print("READ {0} RECORDS.".format(relevant_read_cnt))
    yield read_acc, read_mems_tmp



def read_seeds(seeds):
    curr_acc = ''
    curr_acc_rev = ''
    read_mems = defaultdict(list)
    read_mems_rev = defaultdict(list)
    is_rc = False
    nr_reads = 0
    for i, line in enumerate(open(seeds, 'r')):
        if line[0] == '>':
            if curr_acc and curr_acc_rev:
                for chr_id in list(read_mems.keys()):
                    coordinate_sorted_tuples = sorted(read_mems[chr_id], key = lambda x: x[1])
                    sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
                    read_mems[chr_id] = sorted_mems

                for chr_id in list(read_mems_rev.keys()):
                    coordinate_sorted_tuples = sorted(read_mems_rev[chr_id], key = lambda x: x[1])
                    sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
                    read_mems_rev[chr_id] = sorted_mems

                yield curr_acc, read_mems, curr_acc_rev, read_mems_rev

                # Reset
                curr_acc = ''
                curr_acc_rev = ''
                read_mems = defaultdict(list)
                read_mems_rev = defaultdict(list)   
                is_rc = False

            if 'Reverse' in line:
                curr_acc_rev = line[1:].strip()
                is_rc = True
            else:
                curr_acc = line[1:].strip()
                is_rc = False

            nr_reads += 1

        else:
            vals =  line.split() #11404_11606           1     11405       202
            exon_part_id = vals[0]
            chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('^')
            # chr_id, ref_coord_start, ref_coord_end = ['1', '1', '1'] # CURRENT DUMMY LINE FOR TESTING OUTSIDE ULTRA'S FORMAT
            chr_id = int(chr_id)
            mem_len = int(vals[3])
            
            
            mem_ref_exon_part_start = int(vals[1]) - 1 # convert to 0-indexed reference as in python
            mem_read_start = int(vals[2]) - 1
            ref_coord_start = int(ref_coord_start) # has already been 0-indexed when constructing parts
            mem_genome_start = ref_coord_start + mem_ref_exon_part_start
            
            info_tuple = ( mem_genome_start, mem_genome_start + mem_len - 1,
                            mem_read_start, mem_read_start + mem_len - 1, 
                            mem_len, exon_part_id) # however, for MEM length last coordinate is inclusive of the hit in MEM solvers, not as in python end-indexing
            if is_rc:
                read_mems_rev[chr_id].append( info_tuple )
            else:
                read_mems[chr_id].append( info_tuple )


    # Last record
    if curr_acc and curr_acc_rev:
        for chr_id in list(read_mems.keys()):
            coordinate_sorted_tuples = sorted(read_mems[chr_id], key = lambda x: x[1])
            sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
            read_mems[chr_id] = sorted_mems

        for chr_id in list(read_mems_rev.keys()):
            coordinate_sorted_tuples = sorted(read_mems_rev[chr_id], key = lambda x: x[1])
            sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
            read_mems_rev[chr_id] = sorted_mems
        
        yield curr_acc, read_mems, curr_acc_rev, read_mems_rev

    print("READ {0} RECORDS (FW and RC counted as 2 records).".format(nr_reads))