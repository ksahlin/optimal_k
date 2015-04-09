#ifndef __Unitiger_utils_h
#define __Unitiger_utils_h

inline void merge_kmer_at_end (string &s, const string &kmer)
{
    s = s + kmer[kmer.size()-1];
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

inline string int_to_string(size_t x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}

inline uint64_t count_nodes(const Graph& graph, uint64_t nbCores)
{
    uint64_t nb_kmers = 0;

    // We get an iterator for all nodes of the graph.
    Graph::Iterator<Node> it = graph.iterator<Node>();

    Dispatcher dispatcher (nbCores, 1);
    dispatcher.iterate (it, [&] (Node n)  {  __sync_fetch_and_add (&nb_kmers, 1);  });

    // // We loop each node. 
    // for (it.first(); !it.isDone(); it.next())
    // {
    //     nb_kmers++;
    // }

    return nb_kmers;

    // ofstream ofs;
    // ofs.open("temp.txt");
    // ofs << graph.getInfo();
    // ofs.close();
    
    // int line_number = 1;
    // string s;

    // ifstream ifs;
    // ifs.open("temp.txt");

    // while (getline(ifs, s))
    // {
    //     cout << line_number << endl;
    //     cout << line_number << " " << s << endl;
    //     if (line_number == 37)
    //     {
    //         cout << "xx " << s << endl;
    //         stringstream ss2(s); 
    //         getline(ss2, s, ':');
    //         cout << "xx " << s << endl;
    //         getline(ss2, s);
    //         cout << "xx " << s << endl;
    //         return stoi(s);
    //     }
    //     line_number++;
    // }
    // ifs.close();
    //remove("temp.txt");

    return 0;
}

inline size_t count_arcs(const Graph& graph)
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

inline bool is_unary_node(const Graph& graph, const Node &node)
{
    return (graph.outdegree(node) == 1) and (graph.indegree(node) == 1);
}


inline void initialize_dummy_graph(Graph& graph, size_t& k)
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

// Check if a file is readable
bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 

inline int initialize_de_bruijn_graph(Graph& graph, 
    string reads, 
    size_t k, 
    size_t abundance, 
    size_t nb_cores,
    bool load_graph)
{

    // //string reads_tmp = reads + ".h5";
    // string reads_tmp = "frag_1.h5";
    // if (load_graph and is_readable(reads_tmp))
    // {
    //     std::cout << "Loading from " << reads_tmp << std::endl;
    //     graph = Graph::load("frag_1");
    //     std::cout << graph.getInfo() << std::endl;
    //     return EXIT_SUCCESS;
    // }

    graph = Graph::create ((char const *)"-in %s -kmer-size %d -abundance-min %d -verbose 0 -nb-cores %d", reads.c_str(), k, abundance, nb_cores);
  
	std::cout << graph.getInfo() << std::endl;

	return EXIT_SUCCESS;
}

#endif
