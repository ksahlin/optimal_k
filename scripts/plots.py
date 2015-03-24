
import argparse
import sys
import os
import csv
import pandas as pd
import itertools

try:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_palette("husl", desat=.6)
except ImportError:
    print "You don't have matplotlib or seaborn installed, or no access to at least one of them"
    sys.exit()

class ResultContainer(object):
    """docstring for ResultContainer"""
    def __init__(self, name):
        super(ResultContainer, self).__init__()
        self.name = name
        
    def read_in_result_file(self,result_file):
        self.results = pd.read_csv(result_file)


def scatterplot(x_axis,y_axis,methods,outfolder):
    for method in methods:
        # print type(method.results)
        # print x_axis, y_axis
        # try:
        #     ax = method.results.plot( x=x_axis, y=y_axis)
        # except AttributeError:
        #     print 'Skipping to plot x_axis:{0} to y_axis:{1} due to ValueError. \
        #     probably because csv column contains ".". This is expected for e.g. "estiamator". '.format(x_axis,y_axis)
        #     return
        try:
            plt.plot(method.results[x_axis], method.results[y_axis], '-', label=method.name )
        except ValueError:
            print 'Skipping to plot x_axis:{0} to y_axis:{1} due to ValueError. \
            probably because csv column contains ".". This is expected for e.g. "estiamator". '.format(x_axis,y_axis)

    # fig = ax.get_figure()
    # fig.savefig(os.path.join(outfolder,'x_axis='+x_axis+',y_axis='+y_axis))

    plt.ylabel(y_axis)
    plt.xlabel(x_axis)
    title = ""
    plt.title(title)
    plt.legend( )
    plt.grid()
    plt.savefig(os.path.join(outfolder,'x_axis='+x_axis+',y_axis='+y_axis))
    plt.close()
    plt.clf()

def histogram():
    pass

def main(args):
    methods= []
    for i,result_file_path in enumerate(args.result_files):
        method_results = ResultContainer(args.names[i])
        method_results.read_in_result_file(result_file_path)
        methods.append(method_results)

    axis_iter = itertools.combinations(['k','a','nr_nodes','avg_internal_nodes','avg_length_unitigs','e_size'], 2)
    for x,y in axis_iter:
        scatterplot(x,y,methods,args.outfolder)


if __name__ == '__main__':

    # create the top-level parser
    parser = argparse.ArgumentParser('Plot script for optimal k')
    subparsers = parser.add_subparsers(help='help for subcommand')

    # create the parser for the "all_plots" command
    all_plots = subparsers.add_parser('all_plots', help='Draw all plots')
    all_plots.add_argument('--result_files', dest='result_files', type=str, nargs='+', help='Paths to the result files. ')
    all_plots.add_argument('--names', dest='names',type=str, nargs='+', help='One name (label) for each result file. This will be used in the figure legend. e.g. "estimates" "validator" "unitiger" ')
    all_plots.add_argument('--outfolder', dest='outfolder', type=str, help='Outfolder. ')
    all_plots.set_defaults(which='all_plots')

    # group = all_plots.add_mutually_exclusive_group(required=True)
    # group.add_argument('--fastq', type=str, help='Fastq file with reads ')
    # group.add_argument('--fasta', type=str, help='Fasta file with reads ')
    # # create the parser for the "nr_nodes" command    
    # nr_nodes_parser = subparsers.add_parser('nr_nodes', help='Filters bam file for better uniform coverage.')
    # nr_nodes_parser.add_argument('bampath', type=str, help='bam file with mapped reads. ')
    # nr_nodes_parser.add_argument('outfolder', type=str, help='Outfolder. ')
    # nr_nodes_parser.set_defaults(which='nr_nodes')


    
    args = parser.parse_args()

    assert len(args.result_files) == len(args.names)

    # if args.which == 'pipeline' or args.which == 'get_bp_stats' or args.which == 'lib_est':
    #     try:
    #         open(args.bampath)
    #     except IOError as e:
    #         sys.exit("couldn't find BAM file: " + args.bampath + " check that the path is correct and that the file exists")
    #     try:
    #         open(args.bampath + '.bai')
    #     except IOError as e:
    #         print "couldn't find index file: ", args.bampath + '.bai', " check that the path is correct and that the bam file is sorted and indexed"
    #         sys.exit(0)
    
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    if args.which == 'all_plots':
        main(args)
    # elif args.which == 'nr_nodes_parser':
    #     nr_nodes(args)
    else:
        print 'invalid call'


