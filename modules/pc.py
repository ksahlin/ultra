import sys
import multiprocessing as mp
from queue import Empty

from collections import defaultdict, namedtuple

from modules import align
from modules import seed_wrapper
from modules import help_functions



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
    outfile = open(outfile_name, 'a')
    tot_written = 0
    batch = []
    read_cnt = 1
    batch_id = 1
    reads_aln = []
    buffer_write_cnt = 0
    # generate reads and their seeds
    for (acc, (seq, _)), (r_acc, read_mems, r_acc_rev, r_mems_rev) in zip(help_functions.readfq(open(reads,"r")), seed_wrapper.read_seeds(seeds)):
        assert acc == r_acc
        batch.append([acc, seq, read_mems, r_mems_rev])

        if read_cnt % 1000 == 0:
            input_queue.put([batch_id, batch])
            batch = []
            batch_id += 1
        read_cnt += 1

        # reads_aln.append('\t'.join([s for s in [acc,seq, r_acc, r_acc_rev]])) 
        # if len(reads_aln) > 1000:
        #     output_sam_buffer.put(reads_aln)
        #     reads_aln = []
 
        if buffer_write_cnt >= 50000:
           tot_written = write(outfile, output_sam_buffer, tot_written)
           buffer_write_cnt = 0

    # last batch
    input_queue.put((batch_id, batch))
    input_queue.put(None)
    tot_written = write(outfile, output_sam_buffer, tot_written)
    print('file_IO: Reading records done. Tot read:', read_cnt - 1)
    print('file_IO: Written records in producer process:', tot_written)
    outfile.close()
    return tot_written



class Managers:
    def __init__(self, reads, seeds, outfile_name, n_proc, args):
        self.reads = reads
        self.seeds = seeds    
        self.outfile_name = outfile_name    
        self.m = mp.Manager()
        self.input_queue = self.m.Queue(200)
        self.output_sam_buffer = self.m.Queue()
        self.classification_and_aln_cov = self.m.Queue()
        self.n_proc = n_proc
        self.args = args

    def start(self):
        self.p = mp.Process(target=file_IO, args=(self.input_queue, self.reads, self.seeds, self.output_sam_buffer, self.outfile_name))
        self.p.start()
        self.workers = [mp.Process(target=align.align_single, args=(i, self.input_queue, self.output_sam_buffer, self.classification_and_aln_cov, self.args)) for i in range(self.n_proc)]
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

        tot_counts = [0, 0, 0, 0, 0, 0, 0, 0] # entries: [aln_cov, 'FSM', 'unaligned', 'NO_SPLICE', 'Insufficient_junction_coverage_unclassified', 'ISM/NIC_known', 'NIC_novel', 'NNC']
        while True:
            try:
                res = self.classification_and_aln_cov.get( False )
                for i in range(len(res)):
                    tot_counts[i] += res[i]
            except Empty:
                break
        print('Done joining processes.')
        return tot_counts


def main(reads, seeds, outfile, args):

    # # for profiling
    # m = mp.Manager()
    # input_queue = m.Queue(1000)
    # output_sam_buffer = m.Queue()
    # file_IO(input_queue, reads, seeds, output_sam_buffer, outfile)
    # sys.exit()
    # ############

    m = Managers(reads, seeds, outfile, args.nr_cores, args)
    m.start()
    tot_counts = m.join()
    return tot_counts
