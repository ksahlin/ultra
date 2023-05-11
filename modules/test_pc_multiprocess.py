import argparse
import multiprocessing as mp
from queue import Empty
import dill as pickle

from collections import namedtuple, defaultdict

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', 'j', "exon_part_id"])
globals()[mem.__name__] = mem # Global needed for multiprocessing

'''
    readfq function from https://github.com/lh3/readfq/blob/master/readfq.py
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


def write(outfile, output_sam_buffer, tot_written):
    rec_cnt = 0
    # print('HEREEE')
    while True:
        try:
            record = output_sam_buffer.get( False )
            outfile.write(''.join([r for r in record]))
            rec_cnt += 1
            tot_written += len(record)
        except Empty:
            print("Wrote {0} batches of reads.".format(rec_cnt))
            break
    return tot_written

def file_IO(input_queue, reads, seeds, output_sam_buffer, outfile_name):
    outfile = open(outfile_name, 'w')
    tot_written = 0
    batch = []
    read_cnt = 1
    batch_id = 1
    # generate reads and their seeds
    for (acc, (seq, _)), (r_acc, read_mems, r_acc_rev, r_mems_rev) in zip(readfq(open(reads,"r")), read_seeds(seeds)):
        assert acc == r_acc
        batch.append((acc, seq, read_mems, r_mems_rev))

        if read_cnt % 10 == 0:
            input_queue.put((batch_id, batch))
            batch = []
            batch_id += 1
        read_cnt += 1

        if output_sam_buffer.qsize() >= 100:
            tot_written = write(outfile, output_sam_buffer, tot_written)

    # last batch
    input_queue.put((batch_id, batch))
    # print(batch)
    # print()
    # print(input_queue.qsize())
    input_queue.put(None)
    # print(input_queue.qsize())
    tot_written = write(outfile, output_sam_buffer, tot_written)
    print('file_IO: Reading records done. Tot read:', read_cnt - 1)
    print('file_IO: Written records in producer process:', tot_written)
    outfile.close()
    return tot_written


def consumer(c_id, input_queue, output_sam_buffer):
    while True:
        batch = input_queue.get()
        # check for stop
        if batch is None:
            # add the signal back for other consumers
            input_queue.put(batch)
            # stop running
            break

        output = []
        # Simulate alignment effort
        p = 'AACGAGCTAGGTCAGGCATCACTGCGTA'
        r = 'AACGAGCTAGTTCAGGAATCACTGCGTA'
        for i in range(100):
            h = sum([1 for n1, n2 in zip(p,r) if n1==n2])
        for b in batch[1]:
            output.append(b[0] + '\t' + p+'\t'+r+'\t'+ str(h) + '\n')
        output_sam_buffer.put(output)
        print('Consumer {0} got batch id {1}, size now {2}'.format(c_id, batch[0], input_queue.qsize()))




class Managers:
    def __init__(self, reads, seeds, outfile_name, n_proc):
        self.reads = reads
        self.seeds = seeds    
        self.outfile_name = outfile_name    
        self.m = mp.Manager()
        self.input_queue = self.m.Queue(200)
        self.output_sam_buffer = self.m.Queue()
        self.n_proc = n_proc

    def start(self):
        self.p = mp.Process(target=file_IO, args=(self.input_queue, self.reads, self.seeds, self.output_sam_buffer, self.outfile_name))
        self.p.start()
        self.workers = [mp.Process(target=consumer, args=(i, self.input_queue, self.output_sam_buffer)) for i in range(self.n_proc)]
        for w in self.workers:
            w.start()

    def join(self):
        for w in self.workers:
            w.join()

        self.p.join()

        f = open(self.outfile_name,'a')
        tot_written = write(f, self.output_sam_buffer, 0)
        f.close()
        print('file_IO: Remainig written records after consumer join:', tot_written)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Parses seeds")
    parser.add_argument('nams', type=str, help='seeds')
    parser.add_argument('reads', type=str, help='reads fast(a/q)')
    parser.add_argument('outfile', type=str, help='SAM file output')
    parser.add_argument('--t', type=int, default= 1, help='Nr processes')
    args = parser.parse_args()

    m = Managers(args.reads, args.nams, args.outfile, args.t)
    m.start()
    m.join()

