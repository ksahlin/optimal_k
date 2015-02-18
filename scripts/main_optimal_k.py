import sys,os
import subprocess

import argparse

import shellinford
from Bio import SeqIO


import random
from collections import Counter
import matplotlib.pyplot as plt

def is_internal(sample,fm,a):
    # Abundance is over a here already good kmers reach here
    in_nbrs = 0
    for nucl in 'ACGT':
        hits = sum([ hit.count[0] for hit in fm.search(nucl+sample[:-1]) ])
        if hits >= a:
            in_nbrs += 1

    if in_nbrs != 1:
        return False


    # unitig starts
    else:
        out_nbrs = 0
        for nucl in 'ACGT':
            hits = sum([ hit.count[0] for hit in fm.search(sample[1:]+nucl) ])
            if hits >= a:
                out_nbrs += 1
        if out_nbrs != 1:

            return False

    # normal "inner kmers"
    return True


def build_fm_index(fm, args):
    index = []
    #iter_thresh= 100
    #i = 0
    readfile = open(args.readfile)
    if args.fasta:
        read_file = SeqIO.parse( readfile , "fasta")
    else:
        read_file = SeqIO.parse( readfile , "fastq")
    for seq_record in read_file:

        #i+=1
        #if i >= iter_thresh:
        #    break

        read =  str(seq_record.seq)
        index.append(read) 
        read_rc =  str(seq_record.seq.reverse_complement())
        index.append(read_rc) 
        
    fm.build(index, args.readfile +'.fm')
    readfile.close() #seek(0)
    print 'finished creating fm index'

def load_fm_index(fm,args):
    #print os.path.join(args.readfile,'.fm')
    fm.read(args.readfile+'.fm')

def calc_sample_weight(fm, sample, a):
    abundance = sum([ hit.count[0] for hit in fm.search(sample) ])
    weight = 1 / float( abundance ) # / float(total_kmers)
    return  weight
    
def number_edges(sample,fm,a):
    edges = 0
    for nucl in 'ACGT':
        hits = sum([ hit.count[0] for hit in fm.search(nucl+sample[:-1]) ])
        # out_edges = sum([ hit.count[0] for hit in fm.search(sample[1:]+nucl) ])
        if hits >= a:
            edges += 1 #3e+ out_edges
    return edges

def get_sample_size(p_est_k, delta_pi):
    delta_max = delta_pi/float(2 + delta_pi)
    delta_p_external_k_plus_one = float(p_est_k * delta_max)
    #delta_p_external_k_plus_one = delta_p_external*p_est_k
    

    #p_est_k_plus_one = (1.96/
    sample_size =  (1.96/delta_p_external_k_plus_one)**2 * (1 - p_est_k)*p_est_k
    return sample_size

def sample_internal_nodes(args, fm, k, a , r, n, sample_size):
    sample_proportion = 0.2
    sample_nr = 0
    max_weight = 1/ float(a) 

    readfile = open(args.readfile)

    #avg_kmer_abundance = n*(r-k+1)/args.genome_length
    total_kmers = n*(r-k+1)
    #required_avg_tightness = ((r - k +1)/ float(a) )
    total_internal = 0
    total_non_internal = 0
    total_weight = 0

    kmers_above_abundance = 0
    kmers_below_abundance = 0
    edges_above_abundance = 0
    edges_below_abundance = 0
    if args.fasta:
        read_file = SeqIO.parse( readfile , "fasta")
    else:
        read_file = SeqIO.parse( readfile , "fastq")

    for seq_record in read_file:
        
        if sample_nr >= sample_size:
            break

        if random.random() < sample_proportion:
            #samples.append(read)
            if random.random() < 0.5:
                read = str(seq_record.seq)
            else:
                read =  str(seq_record.seq.reverse_complement())

            pos = random.randint(0, len(seq_record.seq)-k+1)
            sample = read[pos:pos+k]

            if sample.count('N') > 0:
                continue
            
            sample_weight = calc_sample_weight(fm, sample, a)
            edge_number = number_edges(sample,fm,a)
            #internal = is_internal(sample,fm,a)
            #print sample_weight
            if sample_weight <= max_weight :
                total_weight += sample_weight
                #tot_kmers_sampled += 1 
                edges_above_abundance += edge_number * sample_weight
                kmers_above_abundance += 1 * sample_weight
                internal = is_internal(sample,fm,a)
                #print internal

                if internal:
                    total_internal += 1 * sample_weight
                else:
                    total_non_internal += 1 * sample_weight

                sample_nr += 1   
            
            else:
                edges_below_abundance += edge_number * sample_weight
                kmers_below_abundance += 1 * sample_weight


            
    estimated_kmers_in_graph = (total_kmers/sample_size)*kmers_above_abundance
    estimated_kmers_erroneos = (total_kmers/sample_size)*kmers_below_abundance

    estimated_edges_in_graph = (total_kmers/sample_size)*edges_above_abundance
    estimated_edges_erroneos = (total_kmers/sample_size)*edges_below_abundance

    readfile.close()
    try:
        objective = total_internal/(float(total_non_internal)/2) #+ k + 1
    except ZeroDivisionError:
        # We did not sample any unitig breaks at all so the
        # best prediction is that the whole genome is in one unitig
        objective = args.genome_length 

    print 'OBJ:',objective
    print 'Estimated nodes in graph:', estimated_kmers_in_graph
    print 'Estimated kmers not in graph:', estimated_kmers_erroneos
    print 'Estimated edges in graph:', estimated_edges_in_graph
    print 'Estimated edges not in graph:', estimated_edges_erroneos

    return objective

def main(args):
    fm = shellinford.FMIndex()
    if args.load_index:
        load_fm_index(fm, args)
    #fm.read(args.fm)
    else:
        build_fm_index(fm,args)
  
    min_function_vals = {}
    p_est_k = 0.5 # initial for proportion of external nodes in sample
    delta_pi = 0.1 # maximum error 10% of our estimator of average nr of nodes in unitig
    for k in range(args.min_k, 50,2): #[50,60,70,75,80]: #range(args.min_k,80,5): #[50,60,70,75,80]: #
        print 'k={0}'.format(k)
        sample_size = get_sample_size(p_est_k, delta_pi)
        print 'Determined smaple size: {0}'.format(sample_size)
        objective_value = sample_internal_nodes(args, fm, k, args.a, args.read_length , args.n, sample_size)
        if not os.path.exists('/tmp/minia_supporting_kmers_{0}/'.format(k)):
            os.makedirs('/tmp/minia_supporting_kmers_{0}/'.format(k))

        min_function_vals[k] = objective_value 
        p_est_k = 1/float(objective_value)


    print min_function_vals
    x,y = zip(*min_function_vals.iteritems())
    plt.plot(x, y,'go')
    plt.savefig(os.path.join(args.outfolder,'estimated_nr_nodes.png'))
    out_file = open(os.path.join(args.outfolder,'estimated_values.txt'),'w')
    for k,v in sorted(min_function_vals.items()):
        print >> out_file, '{0}\t{1}'.format(k,v)


if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Finds optimal k ")
    parser.add_argument('readfile', type=str, help='Fast(a/q) file. ')
    parser.add_argument('a', type=int, help='Abundance. ')
    parser.add_argument('min_k', type=int, help='Minimum k. ')
    #parser.add_argument('fm', type=str, help='FM-index file.')
    read_file_args = parser.add_mutually_exclusive_group(required=True)
    read_file_args.add_argument('--fasta', dest='fasta', action='store_true', help='Fasta file. ')
    read_file_args.add_argument('--fastq', dest='fastq', action='store_true', help='Fastq file. ')

    fm_index_args = parser.add_mutually_exclusive_group(required=True)
    fm_index_args.add_argument('--build_index', dest='build_index', action='store_true',  help='Buid FM index. ')
    fm_index_args.add_argument('--load_index', dest='load_index',action='store_true', help='Buid FM index. ')

    parser.add_argument('outfolder', type=str, help='outfolder. ')
    parser.add_argument('genome_length', type=int, help='Estimated genome size. ')

    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    parser.add_argument('n', type=int, help='Total number of reads. ')

    args = parser.parse_args()
    main(args)
