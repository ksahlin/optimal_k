import sys,os
import subprocess

import argparse

from Bio import SeqIO


import random
from collections import Counter
import matplotlib.pyplot as plt


def get_kmer_counter(args,k):
    kmer_counter = Counter()
    reference_file = open(args.reference_file)
    for seq_record in SeqIO.parse(reference_file, "fasta"):
        ref =  str(seq_record.seq)
        ref_rc =  str(seq_record.seq.reverse_complement())
        for i in range(0,len(ref)-k+1):
            kmer_counter[ ref[i:i+k] ] += 1
            kmer_counter[ ref_rc[i:i+k] ] += 1
            #if i >= 1000:
            #    break
        
    #print kmer_counter
    reference_file.close() #seek(0)
    print 'finished counting k-mers'
    return kmer_counter

def plot_abundances(kmer_dict, args, k):
    abundance_counts = map(lambda x: kmer_dict[x], kmer_dict)
    plt.hist(abundance_counts, log=True)
    plt.savefig(os.path.join(args.outfolder, 'random_genome_k_{0}'.format(k)))

def main(args):
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    stats_file = open(os.path.join(args.outfolder, 'stats.txt'), 'w')

    for k in range(15, 31, 5):
        print 'k={0}'.format(k)
        kmer_dict = get_kmer_counter(args,k)
        #print kmer_dict
        plot_abundances(kmer_dict, args, k)
        print >> stats_file, 'k={0}: {1}'.format(k,len(kmer_dict))


    # x,y = zip(*min_function_vals.iteritems())
    # plt.plot(x, y,'go')
    # x,y = zip(*e_size_vals.iteritems())
    # plt.plot(x, y,'b^')




if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Calculate and plot kmers on reference ")
    parser.add_argument('reference_file', type=str, help='Fastq file. ')
    parser.add_argument('outfolder', type=str, help='Outfolder. ')
    args = parser.parse_args()
    main(args)
