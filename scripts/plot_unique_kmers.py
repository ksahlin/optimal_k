import sys,os
import subprocess

import argparse

from Bio import SeqIO

from collections import Counter
import matplotlib.pyplot as plt


def get_kmer_counter(fast_file,k,file_type='fasta'):
    kmer_counter = Counter()
    infile = open(fast_file)

    iter_thresh = 0
    for seq_record in SeqIO.parse(infile, file_type):
        iter_thresh += 1
        ref =  str(seq_record.seq)
        ref_rc =  str(seq_record.seq.reverse_complement())
        for i in range(0,len(ref)-k+1):
            kmer = ref[i:i+k].upper()
            kmer_rc = ref_rc[i:i+k].upper()
            if kmer.count('N')  == 0 and  kmer.count('R') == 0 and  kmer.count('Y') == 0:
                kmer_counter[ kmer ] += 1
            if kmer_rc.count('N')  == 0 and  kmer_rc.count('R') == 0 and  kmer_rc.count('Y') == 0:
                kmer_counter[ kmer_rc ] += 1

            # elif kmer not in ['A','C', 'G', 'T']:
            #     print kmer
            # else:
            #     or kmer_rc.count('N')
                
            #     kmer_counter[ kmer_rc ] += 1

        #if iter_thresh >= 1000:
        #   break
        
    #print kmer_counter
    infile.close() #seek(0)
    print 'finished counting k-mers'
    return kmer_counter

def plot_abundances(kmer_dict, args, k):
    abundance_counts = map(lambda x: kmer_dict[x], kmer_dict)
    plt.hist(abundance_counts,bins=20, log=True)
    plt.savefig(os.path.join(args.outfolder, 'k_{0}'.format(k)))

def main(args):
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    stats_file = open(os.path.join(args.outfolder, 'stats.txt'), 'w')
    number_ref_kmers = {}
    number_read_kmers = {}
    
    print >> stats_file, 'k\ton_reference\tin_reads_with_abundance>{0}'.format(args.a)

    for k in range(10, 40, 3):
        print 'k={0}'.format(k)
        kmer_dict_ref = get_kmer_counter(args.reference_file,k)
        #kmer_dict_reads = get_kmer_counter(args.read_file,k,file_type='fastq')
        print kmer_dict_reads, kmer_dict_ref
        number_ref_kmers[k] = len(kmer_dict_ref)
        print filter(lambda x: kmer_dict_reads[x] >= args.a, kmer_dict_reads)
        #number_read_kmers[k] = len(filter(lambda x: kmer_dict_reads[x] >= args.a, kmer_dict_reads)) 

        #plot_abundances(kmer_dict, args, k)
        print >> stats_file, '{0}\t{1}\t{2}'.format(k, number_ref_kmers[k], number_read_kmers)

    ref_x,ref_y = zip(*number_ref_kmers.iteritems())
    #read_x,read_y = zip(*number_read_kmers.iteritems())
    plt.plot(ref_x,ref_y,'go-')
    #plt.plot(read_x,read_y,'b^-')
    plt.savefig(os.path.join(args.outfolder, 'unique_kmers_{0}'.format(k)))





if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Calculate and plot kmers on reference ")
    parser.add_argument('reference_file', type=str, help='Fasta file. ')
    parser.add_argument('read_file', type=str, help='Fasta file. ')
    parser.add_argument('a', type=int, help='Minimum abundance. ')

    #parser.add_argument('fm_index', type=str, help='FM index of reads. ')

    parser.add_argument('outfolder', type=str, help='Outfolder. ')
    args = parser.parse_args()
    main(args)















