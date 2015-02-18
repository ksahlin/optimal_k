// g++ -Wall -O3 -std=c++0x -DMASSIVE_DATA_RLCSA -o main main.cpp rlcsa/rlcsa.a

#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>

#include "rlcsa/rlcsa.h"

using namespace CSA;
using namespace std;

#define twosided95p_quantile 1.96

const char *short_options = "r:b:l:";

static struct option long_options[] = {
// general options
{"readfile",				required_argument,		 0,			 'r'},
{"buildindex",				optional_argument,		 0,			 'b'},
{"loadindex",				optional_argument,		 0,			 'l'},
{0, 0, 0, 0}
};

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

bool is_internal(const string& node, 
	const RLCSA* rlcsa, 
	const int abundance
	)
{
	// at this point, the abundance of 'node' is more than 'abundance'
	int in_nbrs = 0;
	string neighbor;
	pair_type result;
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = nucl + node.substr(0,node.length()-1);
		result = rlcsa->count(neighbor);
	 	if (length(result) >= abundance)
	 	{
	 		in_nbrs++;
	 	}
	 	if (in_nbrs >= 2)
	 	{
	 		return false;
	 	}
	}

    if (in_nbrs != 1)
    {
    	return false;
    }
    else
    { // unitig starts
    	int out_nbrs = 0;
    	for (auto nucl : {'A','C','G','T'})
    	{
    		neighbor = node.substr(1,node.length()-1) + nucl;
    		result = rlcsa->count(neighbor);
    		if (length(result) >= abundance)
    		{
    			out_nbrs++;
    		}
    		if (out_nbrs >= 2)
	 		{
	 			return false;
	 		}
    	}
    	if (out_nbrs != 1)
    	{
    		return false;
    	}
    }
  
  	// normal "inner kmers"
    return true;
}

int calc_abundance(const RLCSA* rlcsa, 
	const string& sample, 
	const int abundance
	)
{
	pair_type result = rlcsa->count(sample);
	return length(result);
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

void sample_internal_nodes(const RLCSA* rlcsa, 
	const int k,
	const int abundance,
	vector<string>& reads,
	const int sample_size,
	double &n_internal,
	double &n_non_internal
	)
{
	double sample_proportion = 0.9;

	uint64_t sampled_so_far = 0;
	uint64_t i = 0;
	string read;
	string sample;

	while (sampled_so_far < sample_size)
	{
		//cout << i << endl;
		if (rand() / (double)RAND_MAX < sample_proportion)
		{
			if (rand() / (double)RAND_MAX < 0.5)
			{
				read = reads[i];
			}
			else
			{
				read = reverse_complement(reads[i]);
			}
			int pos = rand() % (read.length() - k + 1);
            sample = read.substr(pos,k);

            size_t foundN = sample.find('N');
            if (foundN == std::string::npos)
            {
            	int sample_abundance = calc_abundance(rlcsa, sample, abundance);
            	assert(sample_abundance > 0);

            	if (sample_abundance >= abundance)
            	{
            		double sample_weight = 1 / (double)sample_abundance;
				   	if (is_internal(sample,rlcsa,abundance))
	            	{
	            		n_internal += 1 * sample_weight;
	            	}
	            	else
	            	{
	            		n_non_internal += 1 * sample_weight;
	            	}
	            	sampled_so_far++;
            	}
            }
		}
		i++;
		if (i >= reads.size())
		{
			cout << "Reached end of reads" << endl;
			return;
		}
	}
	
	// cout << "Collected " << sampled_so_far << " samples" << endl;
	// cout << "Found " << n_internal << " internal nodes" << endl;
	// cout << "Found " << n_non_internal << " non internal nodes" << endl;

	//return n_internal / (double)(n_non_internal / 2);

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
    int opt;
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName;
    bool buildindex = false;
    bool loadindex = false;
	vector<string> reads;

 	double startTime = readTimer();

	// initializing random number generator
	srand(time(NULL));

	if (argc<3) 
	{
		cerr << "Usage: ./main [--loadindex|--buildindex] --readfile <.fastq file>" << endl;
	 	return 0;
	}

    do 
    {
        opt = getopt_long(argc, argv, short_options, long_options, &opt_index);
        switch (opt) {
			case -1:     /* Done with options. */
				break;
			case 'r':
				readFileName = optarg;
				break;
			case 'b':
				buildindex = true;
				break;
			case 'l':
				loadindex = true;				
				break;
		}
	} while(opt != -1);

	if (((buildindex == false) and (loadindex == false)) or ((buildindex == true) and (loadindex == true)))
	{
		cerr << "You must specify either --buildindex or --loadindex, and not both" << endl;
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
 	get_reads(readFileName, reads);
	
	cout << "Time for loading the index and the reads: " << readTimer() - startTime << "sec" << endl;
 	
 	// we sample internal nodes and check their abundances
 	startTime = readTimer();
 	
 	int abundance = 2;
 	uint64_t sample_size;
    double prop_external_k = 0.5; // initial for proportion of external nodes in sample
    double delta_avg_unitig_length = 0.1; // maximum error 10% of our estimator of average nr of nodes in unitig


 	double average_unitig_length[101];
 	for (int i = 0; i <= 100; i++)
 	{
 		average_unitig_length[i] = 0;
 	}

 	for (int k = 15; k <= 80; k = k + 1)
 	{
 		double n_internal = 0;
 		double n_non_internal = 0;
 		sample_size = get_sample_size(prop_external_k, delta_avg_unitig_length);
 		sample_internal_nodes(rlcsa, k, abundance, reads, sample_size, n_internal, n_non_internal);	
 		average_unitig_length[k] = n_internal / (double)(n_non_internal / 2);

 		cout << k << "," << average_unitig_length[k] << ". ess= " << sample_size << endl;

 		// update estimator
 		prop_external_k = n_non_internal / (n_internal + n_non_internal);
 	}
 	
 	
 	cout << "Time for sampling: " << readTimer() - startTime << "sec" << endl;

	return EXIT_SUCCESS;
}