
#include <gatb/gatb_core.hpp>
#include <string>
#include <fstream>
#include <unordered_set>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "OptionParser.h"

using namespace std;

string merge_kmer_at_end (string s, string kmer)
{
    s = s + kmer[kmer.size()-1];
    return s;
}

string merge_kmer_at_beginning (string s, string kmer)
{
    s = kmer[0] + s;
    return s;
}

inline char reverse_complement_char(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

inline string reverse_complement(string s)
{
    string reverse;

    for (int i = s.length()-1; i >= 0; i--)
    {
        reverse += reverse_complement_char(s[i]);
    }

    return reverse;
}

string int_to_string(size_t x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}

// Check if a file is readable
bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 

void make_upper_case(string& s)
{
    for (size_t i = 0; i < s.length(); i++)
    {
        if (s[i] == 'a') s[i] = 'A';
        else if (s[i] == 'c') s[i] = 'C';
        else if (s[i] == 'g') s[i] = 'G';
        else if (s[i] == 't') s[i] = 'T';
        else if (s[i] == 'n') s[i] = 'N';
    }
}


int initialize_de_bruijn_graph(Graph& graph, string reads, size_t k, size_t abundance, size_t nb_cores)
{

    string readsString = "";
    string line;

    ifstream readFile;
    readFile.open(reads);

    if (readFile.is_open())
    {
        getline(readFile , line);
        readsString = line;

        while (getline(readFile , line)) // this is the comment line
        {
            readsString += "," + line;
        }
        readFile.close();            
    }
    else
    {
        cout << "Error: couldn't open file " << reads << endl;
        return EXIT_FAILURE;
    }

    // Tokenize reads file (list of files separated by ,)
    int filecount=1;
    char *readcstr = (char *)readsString.c_str();
    for(int i = 0; i < strlen(readcstr); i++) 
    {
        if (readcstr[i] == ',')
            filecount++;
    }

    if (filecount > 1) 
    {
        char **files = new char*[filecount];
        files[0] = readcstr;
        int j = 1;
        int l = strlen(readcstr);
        for (int i = 0; i < l; i++) {
            if (readcstr[i] == ',')
            {
                readcstr[i] = '\0';
                files[j] = &readcstr[i+1];
                j++;
            }
        }

        BankFasta *b = new BankFasta(filecount, files);
        graph = Graph::create(b, (char const *)"-kmer-size %d -abundance %d -verbose 0 -nb-cores %d", k, abundance, nb_cores);
    } 
    else 
    {
        graph = Graph::create ((char const *)"-in %s -kmer-size %d -abundance %d -verbose 0 -nb-cores %d", readsString.c_str(), k, abundance, nb_cores);
    }
  
	std::cout << graph.getInfo() << std::endl;


	return EXIT_SUCCESS;
}


size_t count_nodes(const Graph& graph)
{
    size_t nb_kmers = 0;

    // We get an iterator for all nodes of the graph.
    Graph::Iterator<Node> it = graph.iterator<Node> ();
    // We loop each node. 
    for (it.first(); !it.isDone(); it.next())
    {
        nb_kmers++;
    }

    return nb_kmers;
}

size_t count_arcs(const Graph& graph)
{
    // we iterate over every node and cout the size of its out-neighborghood

    size_t nb_arcs = 0;

    // We get an iterator for all nodes of the graph.
    Graph::Iterator<Node> it = graph.iterator<Node> ();
    // We loop each node. 
    for (it.first(); !it.isDone(); it.next())
    {
        nb_arcs += graph.outdegree(it.item());
        // cout << "current node " << graph.toString(it.item()) << " has " << graph.outdegree(it.item()) << " out-neighbors ";
        Node other = graph.reverse (it.item());
        // cout << "its reverse complement is " << graph.toString(other) << " which has " << graph.outdegree(other) << " out-neighbors " << std::endl;
        nb_arcs += graph.outdegree(other);
    }

    return nb_arcs/2;
}

bool is_unary_node(Graph& graph, Node node)
{
    return ((graph.successors<Node>(node).size() == 1) and (graph.predecessors<Node>(node).size() == 1));
}

void print_nodes_and_neighbors(Graph& graph)
{
    std::cout << "All nodes: " << std::endl;
    Graph::Iterator<Node> allit = graph.iterator<Node> ();
    // We loop each node. Note the structure of the for loop.
    for (allit.first(); !allit.isDone(); allit.next())
    {
        // The currently iterated node is available with it.item()
        // We dump an ascii representation of the current node.
        std::cout << "node " << graph.toString (allit.item()) << " has predecessors: ";
        Graph::Vector<Node> in_n = graph.predecessors<Node> (allit.item());
        for (size_t i=0; i < in_n.size(); i++)
        {
            std::cout << graph.toString(in_n[i]) << " ";
        }
        std::cout << " and successors: ";
        Graph::Vector<Node> out_n = graph.successors<Node> (allit.item());
        for (size_t i=0; i < out_n.size(); i++)
        {
            std::cout << graph.toString(out_n[i]) << " ";
        }
        std::cout << endl;
    }
}


void initialize_dummy_graph(Graph& graph, size_t& k)
{
    // k = 4;
    // IBank* bank = new BankStrings (
    //     "AGGTTCA",
    //     "GTTA",
    //     "AACC",
    //     NULL
    // );

    // k = 4;
    // IBank* bank = new BankStrings (
    //     "AGGTTCATG",
    //     "GGTACAT",
    //     "TTACG",
    //     NULL
    // );

    // k = 5;
    // IBank* bank = new BankStrings (
    //     "AGTCATA",
    //     NULL
    // );

    k = 5;
    IBank* bank = new BankStrings (
        "AGTCATC",
        "TCATA",
        NULL
    );

    // We create the graph from a given sequence, and for a given kmer size
    graph = Graph::create (bank,  "-kmer-size %d  -abundance 1 -verbose 0", k);
}


inline unordered_set<string> compute_unitigs(Graph& graph, int nb_cores)
{
    unordered_set<string> unitigs;

    if (nb_cores == 0)
    {
        nb_cores = omp_get_num_procs();
    }

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

    #pragma omp parallel for num_threads(nb_cores)
    for (uint64_t j = 0; j < n_nodes; j++)
    {
        // The currently iterated branching node is available with it.item()
        // We dump an ascii representation of the current node.
        // std::cout << "[" << it.rank() << "] " << graph.toString (it.item()) << std::endl;

        //Graph::Iterator<BranchingNode> it = branching_nodes[j];

        Node current_node = branching_nodes[j];

        // for each out-neighbor, we traverse as long as we see unary nodes
        Graph::Vector<Node> out_neighbors = graph.successors<Node> (current_node);
        for (size_t i = 0; i < out_neighbors.size(); i++)
        {
            string current_contig = graph.toString(current_node);
            // as long as we see unary nodes, we go on
            for (Node current_successor = out_neighbors[i]; ; )
            {
                // std::cout << "visited by successors " << graph.toString(current_successor) << std::endl;
                current_contig = merge_kmer_at_end(current_contig, graph.toString(current_successor));    

                if (is_unary_node(graph,current_successor))
                {
                    current_successor = graph.successors<Node> (current_successor)[0];
                }
                else
                    break;
            }

            // std::cout << ">> Constructed contig: " << current_contig << std::endl;
            if ((unitigs.count(current_contig) == 0) and (unitigs.count(reverse_complement(current_contig)) == 0))
            {
                #pragma omp critical
                {
                    unitigs.insert(current_contig);    
                }
            }
        }

        // for each in-neighbor, we traverse as long as we see unary nodes
        Graph::Vector<Node> in_neighbors = graph.predecessors<Node> (current_node);
        for (size_t i = 0; i < in_neighbors.size(); i++)
        {
            string current_contig = graph.toString(current_node);
            // as long as we see unary nodes, we go on
            for (Node current_predecessor = in_neighbors[i]; ; )
            {
                // std::cout << "visited by predecessors " << graph.toString(current_predecessor) << std::endl;
                current_contig = merge_kmer_at_beginning(current_contig, graph.toString(current_predecessor));    

                if (is_unary_node(graph,current_predecessor))
                {
                    current_predecessor = graph.predecessors<Node> (current_predecessor)[0];
                }
                else
                    break;
            }

            // std::cout << ">> Constructed contig: " << current_contig << std::endl;
            if ((unitigs.count(current_contig) == 0) and (unitigs.count(reverse_complement(current_contig)) == 0))
            {
                #pragma omp critical
                {
                    unitigs.insert(current_contig);    
                }
            }
        }

    }

    return unitigs;
}

double compute_average_length(unordered_set<string>& unitigs)
{
    size_t sum_lengths = 0;
    for (unordered_set<string>::iterator itr = unitigs.begin(); itr != unitigs.end(); ++itr) 
        {
            sum_lengths += (*itr).length();
        }

    return sum_lengths / (double)unitigs.size();
}

double compute_e_size(unordered_set<string>& unitigs)
{
    double e_size;
    size_t sum_lengths = 0;
    size_t sum_lengths_squared = 0;
    for (unordered_set<string>::iterator itr = unitigs.begin(); itr != unitigs.end(); ++itr) 
        {
            int ctg_len = (*itr).length();
            sum_lengths += ctg_len;  //(*itr).length();

            sum_lengths_squared += pow( static_cast<double>(ctg_len) ,2);
        }
    e_size = sum_lengths_squared / (double)sum_lengths;
    return e_size;
}

int print_unitigs(Graph& graph, unordered_set<string>& unitigs, string reads)
{
    try 
    {
        ofstream output_file;

        output_file.open((reads + ".unitigs").c_str());

        size_t i = 1;
        for (unordered_set<string>::iterator itr = unitigs.begin(); itr != unitigs.end(); ++itr) 
        {
            output_file << ">UNITIG" << i << endl;
            output_file << *itr << endl;
            i++;
        }

        output_file.close();    
    }
    catch (gatb::core::system::Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;

}

void print_metrics(const Graph& graph,
    size_t k,
    size_t abundance,
    unordered_set<string>& unitigs,
    ofstream& metricsFile
    )
{
    double average_length = compute_average_length(unitigs);
    double average_internal_nodes = average_length - k - 1;

    metricsFile << k << ",";
    metricsFile << abundance << ",";
    metricsFile << count_nodes(graph) << ","; // number of nodes
    metricsFile << count_arcs(graph) << ","; // number of edges
    metricsFile << average_internal_nodes << ","; // average number of internal nodes in unitigs
    metricsFile << average_length << ","; // average length of unitigs
    metricsFile << ".,"; // estimated sample size
    metricsFile << unitigs.size() << ",";
    metricsFile << compute_e_size(unitigs); // e-size
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
    parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(1) .help("number of threads for graph construction (0 for using cores; default: %default)");


    optparse::Values& options = parser.parse_args(argc, argv);
    readFileName = (string) options.get("r");
    outputFileName = (string) options.get("o");
    k = (size_t) options.get("k");
    abundance = (int) options.get("a");
    nb_cores = (size_t) options.get("t");
    print_output_unitigs = (options.get("not_print_output_unitigs") ? false : true);

    ofstream metricsFile;

    metricsFile.open((outputFileName + ".k" + int_to_string(k) + ".a" + int_to_string(abundance) + ".metrics.csv").c_str());
    metricsFile << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    
    unordered_set<string> unitigs;
    Graph graph;
    try
    {
        initialize_de_bruijn_graph(graph, readFileName, k, abundance, nb_cores);
    }
    catch (gatb::core::system::Exception& e)
    {
        std::cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }
    unitigs = compute_unitigs(graph, nb_cores);
    if (print_output_unitigs)
    {
        print_unitigs(graph, unitigs, outputFileName);
    }
    print_metrics(graph, k, abundance, unitigs, metricsFile);

    metricsFile.close();

    return EXIT_SUCCESS;
}