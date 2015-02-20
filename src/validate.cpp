
#include "utils.h"

using namespace CSA;
using namespace std;

#define twosided95p_quantile 1.96

int N_THREADS;

char reverse_complement_char(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

inline int calc_abundance(const RLCSA* rlcsa, 
	const string& sample
	)
{
	pair_type result = rlcsa->count(sample);
	return length(result);
}

void get_in_out_degrees(const string& node, 
	const RLCSA* rlcsa, 
	const int abundance,
	int &in_degree,
	int &out_degree
	)
{
	// at this point, the abundance of 'node' is more than 'abundance'
	in_degree = 0;
	out_degree = 0;
	string neighbor;
	pair_type result;
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = nucl + node.substr(0,node.length()-1);
	 	if (calc_abundance(rlcsa,neighbor) >= abundance)
	 	{
	 		in_degree++;
	 	}
	}

	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = node.substr(1,node.length()-1) + nucl;
		if (calc_abundance(rlcsa,neighbor) >= abundance)
		{
			out_degree++;
		}
	}
  
}


void sample_all_nodes(const RLCSA* rlcsa, 
	const int k,
	const int abundance,
	vector<string>& reads,
	double &n_internal,
	double &n_starts,
	uint64_t &n_nodes
	)
{
	double kmers_above_abundance = 0;
	uint64_t total_kmers = reads.size() * (reads[0].length() - k + 1);
	uint64_t sample_size = total_kmers;
	double n_internal_local = 0;
	double n_starts_local = 0;

	int n_reads = reads.size();
	
	int progress_reads = -1;

	int n_palindrome_kmers = 0;

	// omp_set_dynamic(0);
	#pragma omp parallel for shared(n_internal,n_starts) num_threads(N_THREADS)
	for (int i = 0; i < n_reads; i = i + 1)
	{
		//#pragma omp critical
		{
			#pragma omp atomic
			progress_reads++;
		}
		if (progress_reads % 100000 == 0)
		{
			cout << "processing read " << progress_reads << "/" << n_reads << endl;
		}

		string read = reads[i];
		for (int pos = 0; pos <= read.length() - k; pos++)
		{
			string sample = read.substr(pos,k);

			// if (sample == reverse_complement(sample))
			// {
			// 	n_palindrome_kmers++;
			// }

	        size_t foundN = sample.find('N');
	        if (foundN == std::string::npos)
	        {
	        	int sample_abundance = calc_abundance(rlcsa, sample);
	        	assert(sample_abundance > 0);

	        	if (sample_abundance >= abundance)
	        	{
	        		double sample_weight = 1 / (double)sample_abundance;
	        		int in_degree, out_degree;
	        		get_in_out_degrees(sample,rlcsa,abundance,in_degree,out_degree);

	        		//#pragma omp critical
	        		{
	        			kmers_above_abundance += 1 * sample_weight;

	        			if ((in_degree == 1) and (out_degree == 1)) // is internal
		            	{
		            		#pragma omp atomic
		            		n_internal_local += 1 * sample_weight;
		            	}

		            	if ((out_degree > 1) or ((out_degree == 1) and (in_degree != 1))) // is start of some unitigs
		            	{
		            		#pragma omp atomic
		            		n_starts_local += 1 * sample_weight * out_degree;
		            	}

		            	if ((in_degree == 0) and (out_degree == 0))
		            	{
		            		#pragma omp atomic
		            		n_starts_local += 1 * sample_weight;
		            	}
	        		}
	        	}
	        }
		}
	}

	n_internal = n_internal_local;
	n_starts = n_starts_local;
	n_nodes = (total_kmers / sample_size) * kmers_above_abundance;

	//cout << "n_palindrome_kmers=" << n_palindrome_kmers << ". " << endl;
}


int main(int argc, char** argv)
{
	int opt_index = 0;
	N_THREADS = 8;
    int opt;
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName;
    string outputFileName;
    bool buildindex = false;
	vector<string> reads;
	int abundance;
	int k;

 	double startTime = readTimer();

	// initializing random number generator
	srand(time(NULL));

	// command line argument parser
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

	parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("input fastq file");
	parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output file");
	parser.add_option("-b", "--buildindex") .action("store_true") .dest("buildindex") .help("if the index on the fastq file is not built");
	parser.add_option("-k", "--kmersize") .type("int") .dest("k") .action("store") .set_default(0) .help("kmer size (default: %default)");
	parser.add_option("-a", "--abundance") .type("int") .dest("a") .action("store") .set_default(3) .help("minimum abundance (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(8) .help("number of threads, in [1..16] (default: %default)");


	optparse::Values& options = parser.parse_args(argc, argv);
	buildindex = (options.get("buildindex") ? true : false);
	readFileName = (string) options.get("r");
	outputFileName = (string) options.get("o");
	k = (size_t) options.get("k");
	abundance = (int) options.get("a");
	N_THREADS = (int) options.get("t");

	// if (k == 0)
	// {
	// 	cerr << "Parameter -k|--kmersize is needed" << endl;
	// 	return EXIT_FAILURE;
	// }

	if (outputFileName == "")
	{
		cerr << "Parameter -o|--outputfile is needed" << endl;
		return EXIT_FAILURE;
	}

	// we need to build the index and exit
	if (buildindex) 
	{
		if (EXIT_FAILURE == get_data_for_rlcsa(readFileName, data, char_count))
		{
			return EXIT_FAILURE;
		}
		// Build RLCSA and report some information.
		RLCSA rlcsa(data, char_count, 32, 0, 1, true); // parameter with value '1' is number of threads; available is compiles with muti-thread support
 		data = 0; // The constructor deleted the data.

		if ( !(rlcsa.isOk()) ) 
 		{
 			return EXIT_FAILURE;
 		}
 		rlcsa.printInfo();
 		rlcsa.reportSize(true);
 		rlcsa.writeTo(readFileName);
 		return EXIT_SUCCESS;
	}

	// we need to load the index
 	const RLCSA* rlcsa = new RLCSA(readFileName, false);
 	if (!(rlcsa->isOk())) 
 	{
 		return EXIT_FAILURE;
 	}
 	rlcsa->printInfo();
 	rlcsa->reportSize(true);

 	// we load the reads
 	if (EXIT_FAILURE == get_reads(readFileName, reads))
 	{
		return EXIT_FAILURE;
	}
	
	cout << "Time for loading the index and the reads: " << readTimer() - startTime << "sec" << endl;
 	
 	// we sample internal nodes and check their abundances
 	startTime = readTimer();
 	
 	uint64_t sample_size;
    double prop_external_k = 0.5; // initial for proportion of external nodes in sample
    double delta_avg_unitig_length = 0.1; // maximum error 10% of our estimator of average nr of nodes in unitig

    ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;

    for (int kk = 15; kk <= 85; kk++)
    {
		double n_internal, n_starts;
 		uint64_t n_nodes;
 		sample_all_nodes(rlcsa, kk, abundance, reads, n_internal, n_starts, n_nodes);

 		cout << "********* k = " << kk << "****************";

 		outputFile << kk << ",";
 		outputFile << abundance << ",";
 		outputFile << n_nodes << ","; // number of nodes
 		outputFile << ".,"; // number of edges
 		outputFile << (int)(n_internal / (double)(n_starts)) << ","; // average number of internal nodes in unitigs
		outputFile << (int)(n_internal / (double)(n_starts)) + kk + 1 << ","; // average length of unitigs
		outputFile << ".,"; // estimated sample size
		outputFile << (int)n_starts << ","; // number of unitigs
		outputFile << "."; // e-size
		outputFile << endl; 
		outputFile.flush();	
    }

    outputFile.close();
 	
 	cout << "CPU time sampling all kmers for all k: " << readTimer() - startTime << "sec" << endl;

	return EXIT_SUCCESS;
}