import sys,os
import subprocess

import argparse

import matplotlib.pyplot as plt



def parse_unitiger_stdout(unitiger_file):
    avg_nr_nodes = 'error occured'
    for line in unitiger_file:
        if line[:3] == 'avg':
            avg_nr_nodes = line.strip().split()[-1]
    return avg_nr_nodes


def calc_E_size(contigfile, genome_length):
    cont_dict = {}
    k = 0
    temp = ''
    accession = ''
    for line in open(contigfile,'r'):
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            cont_dict[accession] = ''
            k += 1
        elif line[0] == '>':
            cont_dict[accession] = temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    cont_dict[accession] = temp

    # if contig_filter_length:
    #     singled_out = 0
    #     for contig in cont_dict.keys():
    #         if len(cont_dict[contig]) < contig_filter_length:
    #             del cont_dict[contig]
    #             singled_out += 1
    #     print >> Information, 'Number of contigs discarded from further analysis (with -filter_contigs set to {0}): {1}'.format(contig_filter_length,singled_out)
    E_tot = 0
    assm_length = 0
    for acc,seq in cont_dict.iteritems():
        E_tot += len(seq)**2
        assm_length += len(seq)
    E_size = E_tot/ assm_length #genome_length
    return(E_size)


def main(args):

    unitig_stats = {}
    e_size_vals = {}
    for k in range(args.min_k, args.max_k + 1,1): #[50,60,70,75,80]: #range(args.min_k,80,5): #[50,60,70,75,80]: #
        print 'k={0}'.format(k)
        #objective_value = sample_internal_nodes(args, fm, k, args.a, args.read_length , args.n)
        if not os.path.exists('/tmp/minia_supporting_kmers_{0}/'.format(k)):
            os.makedirs('/tmp/minia_supporting_kmers_{0}/'.format(k))

        subprocess.check_call([ "unitiger", args.readfile, str(k), str(args.a)], stdout=open('/tmp/unitiger_k_{0}.stdout'.format(k),'w'), stderr=open('/tmp/unitiger_k_{0}.stdout'.format(k),'w'))
        unitiger_result = parse_unitiger_stdout(open('/tmp/unitiger_k_{0}.stdout'.format(k),'r'))

        if os.path.isfile( args.readfile +'.unitigs' ): #'/tmp/minia_out_{0}/minia.contigs.fa'.format(k)): 
            e_size = calc_E_size( args.readfile+'.unitigs', args.genome_length )
        else:
            e_size = 0
        e_size_vals[k] = e_size 

        unitig_stats[k] = unitiger_result 


    print  unitig_stats
    x,y = zip(*unitig_stats.iteritems())
    node_nr , = plt.plot(x, y,'xr-',label="Avg nr of nodes")
    x,y = zip(*e_size_vals.iteritems())
    e_size , = plt.plot(x, y,'.b-', label="E size")
    plt.legend()
    plt.savefig(os.path.join(args.outfolder,'true_nr_nodes.png'))
    out_file = open(os.path.join(args.outfolder,'true_values.txt'),'w')
    for k,v in sorted(unitig_stats.items()):
        print >> out_file, '{0}\t{1}\t{2}'.format(k,v, e_size_vals[k])


if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Finds optimal k ")
    parser.add_argument('readfile', type=str, help='Fastq file. ')
    parser.add_argument('a', type=int, help='Abundance. ')
    parser.add_argument('min_k', type=int, help='Minimum k. ')
    parser.add_argument('max_k', type=int, help='Minimum k. ')
    parser.add_argument('genome_length', type=int, help='Estimated genome length. ')


    parser.add_argument('outfolder', type=str, help='outfolder. ')

    args = parser.parse_args()
    main(args)
