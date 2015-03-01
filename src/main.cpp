
#include "utils.h"

#define TWOSIDED95P_QUANTILE 1.96

int N_THREADS;

inline int calc_abundance(const RLCSA* rlcsa, 
	const string& sample
	)
{
	pair_type result = rlcsa->count(sample);
	pair_type result_rc = rlcsa->count(reverse_complement(sample));
	int abundance = length(result) + length(result_rc);
	return abundance;
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


void sample_nodes(const RLCSA* rlcsa, 
	const int k,
	const int min_abundance,
	const int max_abundance,
	const vector<string>& reads,
	const vector<uint64_t>& sample_size,
	vector<double> &n_internal,
	vector<double> &n_starts,
	vector<double> &n_nodes,
	vector<double> &n_unitigs
	)
{
	uint64_t total_kmers = reads.size() * (reads[0].length() - k + 1);
	vector<double> kmers_above_abundance(max_abundance + 1, 0);
	vector<double> n_internal_local(max_abundance + 1, 0);
	vector<double> n_starts_local(max_abundance + 1, 0);

	uint64_t kmers_tried = 0;
	uint64_t n_reads = reads.size();
	bool sampled_enough = false;
	vector<uint64_t> sampled_so_far(max_abundance + 1, 0);


	vector<int> shuffle_vector;
	// create a random permutation of [0..n_reads-1]
	for (uint64_t i = 0; i < n_reads; i++)
	{
		shuffle_vector.push_back(i);
	}
	std::random_shuffle(shuffle_vector.begin(), shuffle_vector.end());

	int a_w_max_ss = -1;
	uint64_t max_sample_size = 0;
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
	for (uint64_t i = 0; i < n_reads; i++)
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

			// // MAKE SURE THIS IS OK!
			// if (rand() / (double)RAND_MAX < 0.5)
			// {
			// 	read = reverse_complement(read);
			// }

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
		// fix kmers_tried
		assert(kmers_tried > 0);
		n_nodes[a] = total_kmers / kmers_tried * kmers_above_abundance[a];
		n_unitigs[a] = total_kmers / kmers_tried * n_starts[a];
	}
}

uint64_t get_sample_size(
	const double &prop_external_k, 
	const double &delta_avg_unitig_length)
{
	double p = MIN(prop_external_k,0.5);
	//double delta_max = delta_avg_unitig_length*( p ) / (double)( 1 + (1 - p) * delta_avg_unitig_length );
    //double delta_max = delta_avg_unitig_length / (double)(2 + delta_avg_unitig_length);
    double delta_max = delta_avg_unitig_length;
    double delta_p_external_k_plus_one = (double)(p * delta_max);

    assert(delta_p_external_k_plus_one > 0);
    return pow(TWOSIDED95P_QUANTILE / delta_p_external_k_plus_one, 2) * (1 - p) * p;
}

inline uint64_t get_sample_size_for_proportion(
	double p, 	
	double err)
{
	// MAKE SURE THIS IS OK
	if (p < 0.01)
	{
		p = 0.01;
	}
	if (p > 0.99)
	{
		p = 0.99;
	}
    double absolute_err = p * err;
    return abs((double)(TWOSIDED95P_QUANTILE / absolute_err) * (double)(TWOSIDED95P_QUANTILE / absolute_err) * (double)(1 - p) * (double)p);
}

int main(int argc, char** argv)
{
    int mink, maxk;
    uint64_t char_count;
    uchar *data = NULL;
    string readFileName, outputFileName, indexFileName;
    bool buildindex = false;
	vector<string> reads;
	int min_abundance,max_abundance;
	double relative_error; // maximum relative error 10% of our estimators

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
	parser.add_option("-b", "--buildindex") .action("store_true") .dest("buildindex") .help("force the index to be rebuilt, even though it exists");
	parser.add_option("-a", "--minabundance") .type("int") .dest("a") .action("store") .set_default(3) .help("try all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("int") .dest("A") .action("store") .set_default(3) .help("try all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(8) .help("number of threads; use 0 for all cores (default: %default)");
	parser.add_option("-k", "--mink") .type("int") .dest("k") .action("store") .set_default(5) .help("minimum kmer size to try (default: %default)");
	parser.add_option("-K", "--maxk") .type("int") .dest("K") .action("store") .set_default(85) .help("maximum kmer sizeto try (default: %default)");
	parser.add_option("-e", "--relerror") .type("float") .dest("e") .action("store") .set_default(0.1) .help("relative error of the estimations (default: %default)");
	optparse::Values& options = parser.parse_args(argc, argv);

	buildindex = (options.get("buildindex") ? true : false);
	readFileName = (string) options.get("r");
	indexFileName = get_first_token(readFileName) + "+";
	outputFileName = (string) options.get("o");
	min_abundance = (int) options.get("a");
	max_abundance = (int) options.get("A");
	mink = (int) options.get("k");
	maxk = (int) options.get("K");
	N_THREADS = (int) options.get("t");
	if (N_THREADS == 0)
	{
		N_THREADS = omp_get_num_procs();
	}
	relative_error = (double) options.get("e");


	// if the index does not exist we need to build it
	if ((not is_readable(indexFileName)) or buildindex) 
	{
		if (EXIT_FAILURE == get_data_for_rlcsa(readFileName, data, char_count))
		{
			return EXIT_FAILURE;
		}
		// Build RLCSA and report some information.
		RLCSA rlcsa_built(data, char_count, 32, 0, 1, true); // parameter with value '1' is number of threads; available is compiles with muti-thread support
 		data = 0; // The constructor deleted the data.

		if ( !(rlcsa_built.isOk()) ) 
 		{
 			return EXIT_FAILURE;
 		}
 		rlcsa_built.printInfo();
 		rlcsa_built.reportSize(true);
 		rlcsa_built.writeTo(indexFileName);
	}

	if (outputFileName == "")
	{
		cerr << "The argument -o|--outputfile is needed" << endl;
		return EXIT_FAILURE;
	}

	// we need to load the index
 	const RLCSA* rlcsa = new RLCSA(indexFileName, false);
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

 	// these vectors get re-written for each value of k
 	vector<uint64_t> sample_size(max_abundance + 1, 0);
 	vector<double> n_internal(max_abundance + 1,1), n_starts(max_abundance + 1,1);
 	vector<double> n_nodes(max_abundance + 1,1), n_unitigs(max_abundance + 1,1);

    vector<ofstream> outputFile(max_abundance + 1);
    for (int a = min_abundance; a <= max_abundance; a++)
    {
    	outputFile[a].open((outputFileName + ".a" + int_to_string(a) + ".csv").c_str());
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    } 

 	for (int k = mink; k <= maxk; k++)
 	{
 		// getting the sample size
 		for (int a = min_abundance; a <= max_abundance; a++)
 		{
 			uint64_t SS_n_nodes = get_sample_size_for_proportion(n_nodes[a] / (double)(reads.size() * (reads[0].length() - k + 1)), relative_error);
 			uint64_t SS_n_unitigs = get_sample_size_for_proportion(n_unitigs[a] / (double)(reads.size() * (reads[0].length() - k + 1)), relative_error);
 			uint64_t SS_avg_nodes_unitig = get_sample_size_for_proportion(n_starts[a] / (double)(n_internal[a] + n_starts[a]), 1 / (double)(1 + relative_error) - 1);
 			
 			// cout << "SS_n_nodes = " << SS_n_nodes << endl;
 			// cout << "SS_n_unitigs = " << SS_n_unitigs << endl;
 			// cout << "ratio = " << n_starts[a] / (double)(n_internal[a] + n_starts[a]) << endl;
 			// cout << "SS_avg_nodes_unitig = " << SS_avg_nodes_unitig << endl;

 			uint64_t max_sample_size = MAX(SS_n_nodes,SS_n_unitigs);
 			max_sample_size = MAX(max_sample_size,SS_avg_nodes_unitig);
 			sample_size[a] = max_sample_size;
 		}

 		// sampling
 		sample_nodes(rlcsa, k, min_abundance, max_abundance, reads, sample_size, n_internal, n_starts, n_nodes, n_unitigs);	

 		// printing the results
 		for (int a = min_abundance; a <= max_abundance; a++)
 		{
 			assert(n_starts[a] > 0);
			uint64_t avg_nodes_unitig = n_internal[a] / (double)n_starts[a];

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << (uint64_t)n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << ".,"; // number of edges
	 		outputFile[a] << (uint64_t)avg_nodes_unitig << ","; // average number of internal nodes in unitigs
			outputFile[a] << (uint64_t)avg_nodes_unitig + k + 1 << ","; // average length of unitigs
			outputFile[a] << (uint64_t)sample_size[a] << ","; // estimated sample size
			outputFile[a] << (uint64_t)n_unitigs[a] << ","; // number of unitigs
			outputFile[a] << "."; // e-size
			outputFile[a] << endl; 
			outputFile[a].flush();

	 		cout << k << " " << a << " avg internal nodes=" << (uint64_t)avg_nodes_unitig << " avg length=" << (uint64_t)avg_nodes_unitig + k + 1 << " n_nodes=" << (uint64_t)n_nodes[a] << " n_unitigs=" << (uint64_t)n_unitigs[a] << " ess=" << (uint64_t)sample_size[a] << endl;
	 		cout.flush();
 		}
 	}
 	
 	for (int a = 1; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}
 	
 	cout << "Wall clock time: " << readTimer() - startTime << "sec" << endl;

	return EXIT_SUCCESS;
}