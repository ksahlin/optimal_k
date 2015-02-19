#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>
#include <omp.h>
#include <algorithm>

#include "../bin/rlcsa/rlcsa.h"
#include "OptionParser.h"

using namespace CSA;
using namespace std;

#define twosided95p_quantile 1.96

#define MIN(a,b) (a <= b) ? a : b

int N_THREADS;

char reverse_complement_char(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

string reverse_complement(string& s)
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

inline int calc_abundance(const RLCSA* rlcsa, 
	const string& sample
	)
{
	pair_type result = rlcsa->count(sample);
	return length(result);
}

void get_in_out_degrees(const string& node, 
	const RLCSA* rlcsa, 
	const int &min_abundance,
	const int &max_abundance,
	vector<int> &in_degree,
	vector<int> &out_degree
	)
{
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		in_degree[a] = 0;
		out_degree[a] = 0;	
	}
	
	string neighbor;
	pair_type result;
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = nucl + node.substr(0,node.length()-1);
		int neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (int a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		in_degree[a]++;
	 	}
	}

	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = node.substr(1,node.length()-1) + nucl;
		int neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (int a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		out_degree[a]++;
	 	}
	}
}

int get_reads(const string readFileName, 
	vector<string>& reads
	)
{
	try 
	{
		ifstream readFile;
		readFile.open(readFileName);
	   	string line;

	   	while (getline(readFile , line)) // this is the comment line
	   	{
	    	getline(readFile , line); // the actual read
	    	reads.push_back(line);
	    	//reads.push_back(reverse_complement(line));
	    	assert(getline(readFile , line)); // the +/- sign
	    	assert(getline(readFile , line)); // the quality values
	   	}
	   	readFile.close();
	} catch (exception& error) 
	{ // check if there was any error
		std::cerr << "Error: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

	cout << "Input file contains " << reads.size() << " reads." << endl;

   	return EXIT_SUCCESS;
}

int get_data_for_rlcsa(const string& readFileName, 
	uchar*& data, 
	uint64_t& char_count
	)
{
	ifstream readFile;
	readFile.open(readFileName);
   	char_count = 0;
   	string line;
   	vector<string> reads;

   	// counting total number of reads and their length
   	while (getline(readFile , line)) // this is the comment line
   	{
    	getline(readFile , line); // the actual read
    	reads.push_back(line);
    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
 		reads.push_back(reverse_complement(line));
    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'   	
    	assert(getline(readFile , line)); // the +/- sign
    	assert(getline(readFile , line)); // the quality values
   	}
   	readFile.close();

   	cout << "Input file " << readFileName << " contains " << reads.size() << " reads." << endl;
   	cout << "The temporary data array will have size " << (double)char_count/1000000000 << "GB." << endl;

   	uint64_t i = 0;
	data = new uchar[char_count];
	for (auto read : reads)
	{
    	for (uint64_t j = 0; j < read.length(); j++)
    	{
    		data[i] = (uchar)read[j];
    		i++;
    		
    	}
    	data[i] = '\0';
    	i++;		
	}
   	cout << "Created the data array" << endl;

   	return EXIT_SUCCESS;
}

void sample_nodes(const RLCSA* rlcsa, 
	const int k,
	const int min_abundance,
	const int max_abundance,
	const vector<string>& reads,
	const vector<uint64_t>& sample_size,
	vector<double> &n_internal,
	vector<double> &n_starts,
	vector<uint64_t> &n_nodes
	)
{
	uint64_t total_kmers = reads.size() * (reads[0].length() - k + 1);
	vector<double> kmers_above_abundance(max_abundance + 1, 0);
	vector<double> n_internal_local(max_abundance + 1, 0);
	vector<double> n_starts_local(max_abundance + 1, 0);

	double n_internal_local_a = 0;
	double n_starts_local_a = 0;

	uint64_t kmers_tried = 0;
	uint64_t n_reads = reads.size();
	bool sampled_enough = false;
	vector<uint64_t> sampled_so_far(max_abundance + 1, 0);


	vector<int> shuffle_vector;
	// create a random permutation of [0..n_reads-1]
	for (int i = 0; i < n_reads; i++)
	{
		shuffle_vector.push_back(i);
	}
	std::random_shuffle(shuffle_vector.begin(), shuffle_vector.end());

	int a_w_max_ss = -1;
	int max_sample_size = 0;
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		if (sample_size[a] >= max_sample_size)
		{
			max_sample_size = sample_size[a];
			a_w_max_ss = a;
		}
	}

	// omp_set_dynamic(0);
	#pragma omp parallel for shared(n_internal_local,n_starts_local,sampled_enough,sampled_so_far) num_threads(N_THREADS)
	for (int i = 0; i < n_reads; i++)
	{
		//#pragma omp flush (sampled_enough)
		if (!sampled_enough)
		{
			if (sampled_so_far[a_w_max_ss] >= sample_size[a_w_max_ss])
			{
				sampled_enough = true;
				continue;
			}
			
			string read;
			string sample;

			read = reads[shuffle_vector[i]];

			// MAKE SURE THIS IS OK!
			if (rand() / (double)RAND_MAX < 0.5)
			{
				read = reverse_complement(read);
			}

			int pos = rand() % (read.length() - k + 1);
	        sample = read.substr(pos,k);

	        // fix kmers tried
	        #pragma omp atomic
	        kmers_tried++;

	        // FIX THIS:
        	if (sample.find('N') != string::npos)
        	{
        		continue;
        	}

        	int sample_abundance = calc_abundance(rlcsa, sample);
        	assert(sample_abundance > 0);
        	double sample_weight = 1 / (double)sample_abundance;
        	vector<int> in_degree(max_abundance + 1), out_degree(max_abundance + 1);

        	get_in_out_degrees(sample,rlcsa,min_abundance,max_abundance,in_degree,out_degree);

			int for_limit = MIN(max_abundance,sample_abundance);	

        	for (int a = min_abundance; a <= for_limit; a++)
        	{
        		// cout << k << "," << a << " in_degree: " << in_degree[a] << " out_degree: " << out_degree[a] << endl;
        	}
        	
        	for (int a = min_abundance; a <= for_limit; a++)
        	{
        		#pragma omp critical
        		{
        			kmers_above_abundance[a] += 1 * sample_weight;	
        		}

    			if ((out_degree[a] == 1) and (in_degree[a] == 1)) // is internal
            	{
            		#pragma omp critical
            		{
            			n_internal_local[a] += 1 * sample_weight;
	            	   	sampled_so_far[a]++;
	            	   	// cout << "internal" << endl;
            		}
            	}

				if ((out_degree[a] > 1) or ((out_degree[a] == 1) and (in_degree[a] != 1))) // is start of some unitigs
            	{
            		#pragma omp critical
            		{
            			n_starts_local[a] += 1 * sample_weight * out_degree[a];	
            			//cout << "start" << endl;
            			sampled_so_far[a]++;	
            		}
            	}

            	if ((out_degree[a] == 0) and (in_degree[a] == 0)) // is isolated node
            	{
            		#pragma omp critical
            		{
            			n_starts_local[a] += 1 * sample_weight;	
            			//cout << "isolated" << endl;	
            			sampled_so_far[a]++;
            		}
            	}
        	}

        	// FIX THIS
        	// else
        	// {
        	// 	#pragma omp atomic
        	// 	kmers_below_abundance += 1 / (double)sample_abundance;
        	// }
    	}
	}
	
	if (not sampled_enough)
	{
		//cout << "I didn't sample enough" << endl;
		cout << "I sampled only " << sampled_so_far[a_w_max_ss] << " out of " <<  sample_size[a_w_max_ss] << endl;
	}

	//cout << "sampled_so_far" << sampled_so_far << endl;

	for (int a = min_abundance; a <= max_abundance; a++)
	{
		n_internal[a] = n_internal_local[a];
		n_starts[a] = n_starts_local[a];
		// this should be correct, but it doesn't give the right estimates
		// n_nodes = total_kmers * kmers_above_abundance / (double)(kmers_above_abundance + kmers_below_abundance);
		// fix kmers_tried
		assert(kmers_tried > 0);
		n_nodes[a] = total_kmers / kmers_tried * kmers_above_abundance[a];
	}
}

uint64_t get_sample_size(const double &prop_external_k, 
	const double &delta_avg_unitig_length)
{
    double delta_max = delta_avg_unitig_length / (double)(2 + delta_avg_unitig_length);
    double delta_p_external_k_plus_one = (double)(prop_external_k * delta_max);

    assert(delta_p_external_k_plus_one > 0);
    return pow(twosided95p_quantile / delta_p_external_k_plus_one, 2) * (1 - prop_external_k) * prop_external_k;
}

int main(int argc, char** argv)
{
	int opt_index = 0;
	N_THREADS;
    int opt;
    int mink;
    int maxk;
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName, outputFileName;
    bool buildindex = false;
	vector<string> reads;
	int min_abundance,max_abundance;

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
	parser.add_option("-a", "--minabundance") .type("int") .dest("a") .action("store") .set_default(3) .help("runs for all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("int") .dest("A") .action("store") .set_default(3) .help("runs for all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(8) .help("number of threads, in [1..16] (default: %default)");
	parser.add_option("-k", "--mink") .type("int") .dest("k") .action("store") .set_default(5) .help("minimum kmer size to try (default: %default)");
	parser.add_option("-K", "--maxk") .type("int") .dest("K") .action("store") .set_default(85) .help("maximum kmer sizeto try (default: %default)");

	
	optparse::Values& options = parser.parse_args(argc, argv);
	buildindex = (options.get("buildindex") ? true : false);
	readFileName = (string) options.get("r");
	outputFileName = (string) options.get("o");
	// kmersize = (size_t) options.get("k");
	min_abundance = (int) options.get("a");
	max_abundance = (int) options.get("A");
	mink = (int) options.get("k");
	maxk = (int) options.get("K");
	N_THREADS = (int) options.get("t");

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
 		cout << "Constructed the index successfully. Now run again the program without -b|--buildindex option." << endl;
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
 	get_reads(readFileName, reads);
	
	cout << "Time for loading the index and the reads: " << readTimer() - startTime << "sec" << endl;
 	
 	// we sample internal nodes and check their abundances
 	startTime = readTimer();
 	
 	vector<uint64_t> sample_size(max_abundance + 1, 0);
    vector<double> prop_external_k(max_abundance + 1, 0.5); // initial for proportion of external nodes in sample
    double delta_avg_unitig_length = 0.1; // maximum error 10% of our estimator of average nr of nodes in unitig


    // This is kind of deprecated: update the code to avoid its use
 	vector<double> average_unitig_length(101,0);

    vector<ofstream> outputFile(max_abundance + 1);
    for (int a = 1; a <= max_abundance; a++)
    {
    	outputFile[a].open((outputFileName + ".a" + int_to_string(a)).c_str());
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    } 

 	for (int k = mink; k <= maxk; k++)
 	{
 		vector<double> n_internal(max_abundance + 1,0), n_starts(max_abundance + 1,0);
 		vector<uint64_t> n_nodes(max_abundance + 1,0);

 		for (int a = min_abundance; a <= max_abundance; a++)
 		{
 			sample_size[a] = get_sample_size(prop_external_k[a], delta_avg_unitig_length);	
 		}

 		sample_nodes(rlcsa, k, min_abundance, max_abundance, reads, sample_size, n_internal, n_starts, n_nodes);	

 		for (int a = min_abundance; a <= max_abundance; a++)
 		{
 			assert(n_starts[a] > 0);
			average_unitig_length[k] = n_internal[a] / (double)n_starts[a];

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << ".,"; // number of edges
	 		outputFile[a] << (int)average_unitig_length[k] << ","; // average number of internal nodes in unitigs
			outputFile[a] << (int)average_unitig_length[k] + k + 1 << ","; // average length of unitigs
			outputFile[a] << sample_size[a] << ","; // estimated sample size
			outputFile[a] << (int)n_starts[a] << ","; // number of unitigs
			outputFile[a] << "."; // e-size
			outputFile[a] << endl; 

	 		cout << k << " " << a << " avg internal nodes=" << (int)average_unitig_length[k] << " avg length=" << (int)average_unitig_length[k] + k + 1 << " n_nodes=" << n_nodes[a] << " ess=" << sample_size[a] << endl;
	 		cout.flush();
	 		// update estimator
	 		prop_external_k[a] = (2 * n_starts[a]) / (n_internal[a] + 2 * n_starts[a]);
	 		// prop_external_k[a] = (2 * n_starts[a]) / (n_internal[a] + 2 * n_starts[a]);
 		}

 	}
 	
 	for (int a = 1; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}
 	
 	
 	cout << "Time for sampling: " << readTimer() - startTime << "sec" << endl;

	return EXIT_SUCCESS;
}