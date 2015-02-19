// g++ -Wall -O3 -std=c++0x -DMASSIVE_DATA_RLCSA -o main main.cpp rlcsa/rlcsa.a

#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>
#include <omp.h>
#include <algorithm>

#include "rlcsa/rlcsa.h"
#include "OptionParser.h"

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

string reverse_complement(string& s)
{
    string reverse;

    for (int i = s.length()-1; i >= 0; i--)
    {
        reverse += reverse_complement_char(s[i]);
    }

    return reverse;
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
	const int abundance,
	vector<string>& reads,
	const int sample_size,
	double &n_internal,
	double &n_starts,
	uint64_t &n_nodes
	)
{
	uint64_t sampled_so_far = 0;

	double kmers_above_abundance = 0;
	uint64_t kmers_tried = 0;
	uint64_t total_kmers = reads.size() * (reads[0].length() - k + 1);
	double n_internal_local = 0;
	double n_starts_local = 0;

	int n_reads = reads.size();
	vector<int> shuffle_vector;
	bool sampled_enough = false;

	int n_rejected_N = 0;

	// create a random permutation of [0..n_reads-1]
	for (int i = 0; i < n_reads; i++)
	{
		shuffle_vector.push_back(i);
	}
	std::random_shuffle(shuffle_vector.begin(), shuffle_vector.end());

	
	// omp_set_dynamic(0);
	#pragma omp parallel for shared(n_internal_local,n_starts_local,sampled_enough,sampled_so_far) num_threads(N_THREADS)
	for (int i = 0; i < n_reads; i++)
	{
		// this code is kind of a hack, maybe it can be improved
		// also, it is not guaranteed that reads are selected randomly
		
		#pragma omp flush (sampled_enough)
		if (!sampled_enough)
		{
			//#pragma omp critical
			{
				if (sampled_so_far >= sample_size)
				{
					sampled_enough = true;	
				}
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

	        #pragma omp atomic
	        kmers_tried++;

	        size_t foundN = sample.find('N');
	        if (foundN == std::string::npos) // if sample contains no N
	        {
	        	int sample_abundance = calc_abundance(rlcsa, sample);
	        	assert(sample_abundance > 0);

	        	if (sample_abundance >= abundance)
	        	{
	        		double sample_weight = 1 / (double)sample_abundance;
	        		int in_degree, out_degree;
	        		get_in_out_degrees(sample,rlcsa,abundance,in_degree,out_degree);

	        		//#pragma omp critical
	        		// above line not needed anymore thanks to #pragma omp atomic
	        		{
	        			#pragma omp atomic
	        			kmers_above_abundance += 1 * sample_weight;

	        			if ((out_degree == 1) and (in_degree == 1)) // is internal
		            	{
		            		#pragma omp atomic
		            		n_internal_local += 1 * sample_weight;
		            	}

		            	if ((out_degree > 1) or ((out_degree == 1) and (in_degree != 1))) // is start of some unitigs
		            	{
		            		#pragma omp atomic
		            		n_starts_local += 1 * sample_weight * out_degree;
		            	}

		            	if ((out_degree == 0) and (in_degree == 0)) // is isolated node
		            	{
		            		#pragma omp atomic
		            		n_starts_local += 1 * sample_weight;
		            	}

		            	// update this only when you sample an internal or start node
		            	#pragma omp atomic
		            	sampled_so_far++;
	        		}
	        	}
	        }
    	}
	}
	
	if (not sampled_enough)
	{
		cout << "I sampled only " << sampled_so_far << " out of the " << sample_size << " needed."<< endl;
	}

	//cout << "sampled_so_far" << sampled_so_far << endl;

	n_internal = n_internal_local;
	n_starts = n_starts_local;
	n_nodes = (total_kmers / (double)kmers_tried) * (double)kmers_above_abundance;

}

uint64_t get_sample_size(const double &prop_external_k, 
	const double &delta_avg_unitig_length)
{
    double delta_max = delta_avg_unitig_length / (double)(2 + delta_avg_unitig_length);
    double delta_p_external_k_plus_one = (double)(prop_external_k * delta_max);

    return pow(twosided95p_quantile / delta_p_external_k_plus_one, 2) * (1 - prop_external_k) * prop_external_k;
}

int main(int argc, char** argv)
{
	int opt_index = 0;
	N_THREADS;
    int opt;
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName, outputFileName;
    bool buildindex = false;
	vector<string> reads;
	int abundance = 3;

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
	// parser.add_option("-k", "--kmersize") .type("int") .dest("k") .action("store") .set_default(31) .help("kmer size (default: %default)");
	parser.add_option("-a", "--abundance") .type("int") .dest("a") .action("store") .set_default(3) .help("minimum abundance (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(8) .help("number of threads, in [1..16] (default: %default)");

	
	optparse::Values& options = parser.parse_args(argc, argv);
	buildindex = (options.get("buildindex") ? true : false);
	readFileName = (string) options.get("r");
	outputFileName = (string) options.get("o");
	// kmersize = (size_t) options.get("k");
	abundance = (int) options.get("a");
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
 		cout << "Constructed the index successfully. Now run again the program without -b|--buildindex option.";
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
 	
 	uint64_t sample_size;
    double prop_external_k = 0.5; // initial for proportion of external nodes in sample
    double delta_avg_unitig_length = 0.1; // maximum error 10% of our estimator of average nr of nodes in unitig


 	double average_unitig_length[101];
 	for (int i = 0; i <= 100; i++)
 	{
 		average_unitig_length[i] = 0;
 	}

    ofstream outputFile;
    outputFile.open(outputFileName);
    outputFile << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;

 	for (int kk = 15; kk <= 85; kk++)
 	{
 		double n_internal, n_starts;
 		uint64_t n_nodes;
 		sample_size = get_sample_size(prop_external_k, delta_avg_unitig_length);
 		sample_nodes(rlcsa, kk, abundance, reads, sample_size, n_internal, n_starts, n_nodes);	
 		average_unitig_length[kk] = n_internal / (double)n_starts;

 		//cout << "********* k = " << kk << "****************";

 		outputFile << kk << ",";
 		outputFile << abundance << ",";
 		outputFile << n_nodes << ","; // number of nodes
 		outputFile << ".,"; // number of edges
 		outputFile << (int)average_unitig_length[kk] << ","; // average number of internal nodes in unitigs
		outputFile << (int)average_unitig_length[kk] + kk + 1 << ","; // average length of unitigs
		outputFile << sample_size << ","; // estimated sample size
		outputFile << (int)n_starts << ","; // number of unitigs
		outputFile << "."; // e-size
		outputFile << endl; 
		outputFile.flush();	

 		cout << kk << " avg internal nodes=" << (int)average_unitig_length[kk] << " avg length=" << (int)average_unitig_length[kk] + kk + 1 << " n_nodes=" << n_nodes << " ess=" << sample_size << endl;

 		// update estimator
 		prop_external_k = (2 * n_starts) / (n_internal + 2 * n_starts);
 	}
 	
 	outputFile.close();
 	
 	cout << "Time for sampling: " << readTimer() - startTime << "sec" << endl;

	return EXIT_SUCCESS;
}