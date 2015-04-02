
#include <gatb/gatb_core.hpp>
#include <string>
#include <fstream>
#include <unordered_set>
#include <stdlib.h>
#include <iostream>

#include "OptionParser.h"
#include "utils.h"

uint64_t UNITIG_BUFFER = 1048576;
// uint64_t UNITIG_BUFFER = 32768;

using namespace std;

// We define a functor that will be cloned by the dispatcher
struct ExploreBranchingNodeFunctor { 

    ISynchronizer* synchro;
    Graph &graph;
    unordered_set<string> &fingerprints_of_unitigs;
    ofstream &unitigsFile;
    uint64_t &unitigsCounter;
    uint64_t &sum_unitig_lengths;
    uint64_t &sum2_unitig_lengths;
    bool print_output_unitigs;
    uint64_t &print_times;
    vector<string> &assembled_unitigs;

    ExploreBranchingNodeFunctor (ISynchronizer* synchro, 
        Graph &graph,
        unordered_set<string> &fingerprints_of_unitigs, 
        ofstream &unitigsFile, 
        uint64_t &unitigsCounter,
        uint64_t &sum_unitig_lengths,
        uint64_t &sum2_unitig_lengths,
        bool print_output_unitigs,
        uint64_t &print_times,
        vector<string> &assembled_unitigs)  : 
            synchro(synchro), 
            graph(graph),
            fingerprints_of_unitigs(fingerprints_of_unitigs),
            unitigsFile(unitigsFile),
            unitigsCounter(unitigsCounter),
            sum_unitig_lengths(sum_unitig_lengths),
            sum2_unitig_lengths(sum2_unitig_lengths),
            print_output_unitigs(print_output_unitigs),
            print_times(print_times),
            assembled_unitigs(assembled_unitigs)
        {}

    void operator() (BranchingNode current_node2)
    {
        BranchingNode current_node = current_node2;
        string current_unitig;
        current_unitig.reserve(10000);

        for (uint64_t strand = 0; strand <= 1; strand++)
        {
            if (strand == 1)
            {
                current_node = graph.reverse(current_node);
            }
        
            // for each out-neighbor, we traverse as long as we see unary nodes
            Graph::Vector<Node> out_neighbors = graph.successors<Node>(current_node);
            bool self_loop_isolated = false;
            for (size_t i = 0; i < out_neighbors.size(); i++)
            {
                current_unitig = graph.toString(current_node);
                string last_node;
                string last_node_previous = graph.toString(current_node);
                bool start_of_unitig = true;
                Node current_successor = out_neighbors[i];
                if ( is_unary_node(graph,current_node) and 
                     (current_node == current_successor) )
                {
                    self_loop_isolated = true;
                    /*LOCK*/ synchro->lock();
                    cout << "found isolated self-loop" << endl;
                    /*UNLOCK*/ synchro->unlock();
                }

                // if current extension has not been reported before
                string fingerprint = last_node_previous;
                merge_kmer_at_end(fingerprint, graph.toString(current_successor));
                if (fingerprints_of_unitigs.count(reverse_complement(fingerprint)) != 0)
                {
                    continue;
                }

                // as long as we see unary nodes, we go on
                while (true)
                {
                    last_node = graph.toString(current_successor);
                    merge_kmer_at_end(current_unitig, graph.toString(current_successor));    
                    if (is_unary_node(graph,current_successor) and (not self_loop_isolated))
                    {
                        current_successor = graph.successors<Node> (current_successor)[0];
                        last_node_previous = last_node;
                        last_node = graph.toString(current_successor); 
                    }
                    else
                    {
                        // this is the end of the unitig
                        merge_kmer_at_end(last_node_previous, last_node);
                        
                        /*LOCK*/ synchro->lock();
                        // make sure that no other thread has processed this unitig from the other side in the mean time
                        if (fingerprints_of_unitigs.count(reverse_complement(fingerprint)) != 0)
                        {
                            /*LOCK*/ synchro->unlock();
                            break;
                        }

                        // we are here if we actually assembled a new unitig
                        // /*LOCK*/ synchro->lock();
                        fingerprints_of_unitigs.insert(last_node_previous);
                        unitigsCounter++;
                        if (print_output_unitigs)
                        {
                            assembled_unitigs.push_back(current_unitig);
                            if (unitigsCounter % UNITIG_BUFFER == 0)
                            {
                                uint64_t i = 0;
                                for (auto &unitig : assembled_unitigs)
                                {
                                    unitigsFile << ">UNITIG_" << (print_times * UNITIG_BUFFER + i) << endl;
                                    unitigsFile << unitig << endl;        
                                    i++;
                                }
                                print_times++;
                                assembled_unitigs.clear();
                            }                              
                        }
                        sum_unitig_lengths += current_unitig.length();
                        sum2_unitig_lengths += current_unitig.length() * current_unitig.length();
                        /*UNLOCK*/ synchro->unlock();
                        break;    
                    }
                }
            }

            // if it was isolated node
            if ((graph.outdegree(current_node) == 0) and (graph.indegree(current_node) == 0) and (strand == 0))
            {
                string current_unitig = graph.toString(current_node);
                /*LOCK*/ synchro->lock();
                unitigsCounter++;    
                if (print_output_unitigs)
                {
                    assembled_unitigs.push_back(current_unitig);
                    if (unitigsCounter % UNITIG_BUFFER == 0)
                    {
                        uint64_t i = 0;
                        for (auto &unitig : assembled_unitigs)
                        {
                            unitigsFile << ">UNITIG_" << (print_times * UNITIG_BUFFER + i) << endl;
                            unitigsFile << unitig << endl;        
                            i++;
                        }
                        print_times++;
                        assembled_unitigs.clear();
                    }
                }
                
                sum_unitig_lengths += current_unitig.length();
                sum2_unitig_lengths += current_unitig.length() * current_unitig.length();    
                /*UNLOCK*/ synchro->unlock();
            }
        }
    }
};


void compute_and_print_unitigs(Graph& graph,
    int nbCores,
    ofstream &unitigsFile,
    ofstream &metricsFile,
    bool print_output_unitigs,
    const int &k,
    const int &abundance
    )
{
    uint64_t unitigsCounter = 0;
    uint64_t sum_unitig_lengths = 0;
    uint64_t sum2_unitig_lengths = 0;

    unordered_set<string> fingerprints_of_unitigs;
    vector<string> assembled_unitigs;
    assembled_unitigs.reserve(UNITIG_BUFFER);
    uint64_t print_times = 0;

    // Graph::Iterator<Node> it = graph.iterator<Node> ();
    Graph::Iterator<BranchingNode> iter = graph.iterator<BranchingNode> ();

    ISynchronizer* synchro = System::thread().newSynchronizer();

    // We create a dispatcher configured for 'nbCores' cores.
    Dispatcher dispatcher (nbCores, 1);

    IDispatcher::Status status = dispatcher.iterate (iter, ExploreBranchingNodeFunctor(synchro,
        graph,
        fingerprints_of_unitigs, 
        unitigsFile, 
        unitigsCounter,
        sum_unitig_lengths,
        sum2_unitig_lengths,
        print_output_unitigs,
        print_times,
        assembled_unitigs) );

    uint64_t i = 0;
    for (auto &unitig : assembled_unitigs)
    {
        unitigsFile << ">UNITIG_" << (print_times * UNITIG_BUFFER + i) << endl;
        unitigsFile << unitig << endl;        
        i++;
    }

    std::cout << "we used " << status.nbCores << " cores, traversal time " << (double)status.time/1000 << " sec" << std::endl;

    double average_length = (double)sum_unitig_lengths / unitigsCounter;
    double average_internal_nodes = average_length - k - 1;
    double e_size = (double) sum2_unitig_lengths / sum_unitig_lengths;

    metricsFile << k << ",";
    metricsFile << abundance << ",";
    metricsFile << count_nodes(graph,nbCores) << ","; // number of nodes
    metricsFile << ".,"; // count_arcs(graph) << ","; // number of edges
    metricsFile << average_internal_nodes << ","; // average number of internal nodes in unitigs
    metricsFile << average_length << ","; // average length of unitigs
    metricsFile << ".,"; // estimated sample size
    metricsFile << unitigsCounter << ",";
    metricsFile << e_size; // e-size
    metricsFile << endl; 

}

void compute_and_print_unitigs_OLD(const Graph& graph,
    int nb_cores,
    ofstream &unitigsFile,
    ofstream &metricsFile,
    bool print_output_unitigs,
    const int &k,
    const int &abundance
    )
{
    uint64_t unitigsCounter = 0;
    uint64_t sum_unitig_lengths = 0;
    uint64_t sum2_unitig_lengths = 0;

    unordered_set<string> fingerprints_of_unitigs;

    // if (nb_cores == 0)
    // {
    //     nb_cores = omp_get_num_procs();
    // }

    // we get an iterator over the branching nodes
    Graph::Iterator<BranchingNode> iter = graph.iterator<BranchingNode> ();

    // we gather all iterators in a vector
    vector< Node > branching_nodes;
    // we iterate over the branching nodes
    for (iter.first(); !iter.isDone(); iter.next())
    {
        branching_nodes.push_back(iter.item());
    }
    uint64_t n_nodes = branching_nodes.size();

    // #pragma omp parallel for num_threads(nb_cores)
    for (uint64_t j = 0; j < n_nodes; j++)
    {
        // The currently iterated branching node is available with it.item()
        // We dump an ascii representation of the current node.
        // std::cout << "[" << it.rank() << "] " << graph.toString (it.item()) << std::endl;

        //Graph::Iterator<BranchingNode> it = branching_nodes[j];

        Node current_node;
        for (uint64_t strand = 0; strand <= 1; strand++)
        {
            if (strand == 0)
            {
                current_node = branching_nodes[j];       
            }
            else
            {
                current_node = graph.reverse(branching_nodes[j]);
            }
        
            // for each out-neighbor, we traverse as long as we see unary nodes
            Graph::Vector<Node> out_neighbors = graph.successors<Node> (current_node);
            for (size_t i = 0; i < out_neighbors.size(); i++)
            {
                string current_unitig = graph.toString(current_node);
                string last_node;
                string last_node_previous = graph.toString(current_node);
                bool start_of_unitig = true;
                Node current_successor = out_neighbors[i];

                // if current extension has not been reported before
                string fingerprint = last_node_previous;
                merge_kmer_at_end(fingerprint, graph.toString(current_successor));
                if (fingerprints_of_unitigs.count(reverse_complement(fingerprint)) != 0)
                {
                    continue;
                }

                // as long as we see unary nodes, we go on
                while (true)
                {
                    // std::cout << "visited by successors " << graph.toString(current_successor) << std::endl;
                    last_node = graph.toString(current_successor);
                     
                    merge_kmer_at_end(current_unitig, graph.toString(current_successor));    
                    if (is_unary_node(graph,current_successor))
                    {
                        current_successor = graph.successors<Node> (current_successor)[0];
                        last_node_previous = last_node;
                        last_node = graph.toString(current_successor); 
                    }
                    else
                    {
                        // this is the end of the unitig
                        // we are here if we actually assembled a new unitig
                        merge_kmer_at_end(last_node_previous, last_node);
                        // #pragma omp critical
                        {
                            fingerprints_of_unitigs.insert(last_node_previous);
                            if (print_output_unitigs)
                            {
                                //unitigs.insert(current_unitig);
                                unitigsFile << "UNITIG" << unitigsCounter << endl;
                                unitigsFile << current_unitig << endl;
                                unitigsCounter++;    
                            }
                            sum_unitig_lengths += current_unitig.length();
                            sum2_unitig_lengths += current_unitig.length() * current_unitig.length();
                        }    
                        break;    
                    }
                }
            }

            // if it was isolated node
            if ((graph.outdegree(current_node) == 0) and (graph.indegree(current_node) == 0) and (strand == 0))
            {
                // #pragma omp critical
                {
                    string current_unitig = graph.toString(current_node);
                    if (print_output_unitigs)
                    {
                        //unitigs.insert(graph.toString(current_node));
                        unitigsFile << "UNITIG" << unitigsCounter << endl;
                        unitigsFile << current_unitig << endl;
                    }
                    unitigsCounter++;    
                    sum_unitig_lengths += current_unitig.length();
                    sum2_unitig_lengths += current_unitig.length() * current_unitig.length();    
                }
            }
        }
    }

    double average_length = (double)sum_unitig_lengths / unitigsCounter;
    double average_internal_nodes = average_length - k - 1;
    double e_size = (double) sum2_unitig_lengths / sum_unitig_lengths;

    metricsFile << k << ",";
    metricsFile << abundance << ",";
    metricsFile << count_nodes(graph, nb_cores) << ","; // number of nodes
    metricsFile << ".,"; // count_arcs(graph) << ","; // number of edges
    metricsFile << average_internal_nodes << ","; // average number of internal nodes in unitigs
    metricsFile << average_length << ","; // average length of unitigs
    metricsFile << ".,"; // estimated sample size
    metricsFile << unitigsCounter << ",";
    metricsFile << e_size; // e-size
    metricsFile << endl; 

}


int main (int argc, char* argv[])
{
    size_t k, abundance, nb_cores;
    string readFileName, outputFileName;
    bool print_output_unitigs;

    string usage = "\n  %prog OPTIONS";
    const string version = "%prog 0.1\nCopyright (C) 2014-2015\n"
        "License GPLv3+: GNU GPL version 3 or later "
        "<http://gnu.org/licenses/gpl.html>.\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.";
    const string desc = "";
    const string epilog = "";
    
    optparse::OptionParser parser = optparse::OptionParser()
        .usage(usage)
        .version(version)
        .description(desc)
        .epilog(epilog);

    parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("a file containing a list of FASTA/Q(.gz) file names, one per line");
    parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output file");
    parser.add_option("-k", "--kmersize") .type("int") .dest("k") .action("store") .set_default(31) .help("kmer size (default: %default)");
    parser.add_option("-a", "--abundance") .type("int") .dest("a") .action("store") .set_default(3) .help("minimum abundance (default: %default)");
    parser.add_option("-s", "--silentoutput") .action("store_true") .dest("not_print_output_unitigs") .set_default(false) .help("this option suppresses writing the unitigs to file");
    parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(0) .help("number of threads (0 for using cores; default: %default)");

    optparse::Values& options = parser.parse_args(argc, argv);
    readFileName = (string) options.get("r");
    outputFileName = (string) options.get("o");
    k = (size_t) options.get("k");
    abundance = (int) options.get("a");
    nb_cores = (size_t) options.get("t");
    print_output_unitigs = (options.get("not_print_output_unitigs") ? false : true);

    ofstream metricsFile;
    ofstream unitigsFile;
    string filePrefix = outputFileName + ".k" + int_to_string(k) + ".a" + int_to_string(abundance);

    Graph graph;
    try
    {
        initialize_de_bruijn_graph(graph, readFileName, k, abundance, nb_cores);
        metricsFile.open((filePrefix + ".csv").c_str());
        metricsFile << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
        if (print_output_unitigs)
        {
            unitigsFile.open((filePrefix + ".unitigs").c_str());
        }
    }
    catch (gatb::core::system::Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }
    compute_and_print_unitigs(graph, 
        nb_cores, 
        unitigsFile, 
        metricsFile, 
        print_output_unitigs, 
        k, 
        abundance);
    
    metricsFile.close();
    if (print_output_unitigs)
    {
        unitigsFile.close();
    }

    return EXIT_SUCCESS;
}