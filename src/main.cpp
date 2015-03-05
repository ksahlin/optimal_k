
#include "utils.h"
#include "esize_estimation.h"

#define TWOSIDED95P_QUANTILE 1.96
#define MIN_SAMPLE_SIZE_UNITIGS 500
#define MAX_SAMPLE_SIZE_UNITIGS 5000
#define RATIO_SAMPLE_SIZE_UNITIGS 0.05

int N_THREADS;

inline void get_in_out_degrees(const string& node, 
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
	const vector<uint64_t>& sample_size_kmers,
	const vector<uint64_t>& sample_size_unitigs,
	vector<double> &n_internal,
	vector<double> &n_starts,
	vector<double> &n_nodes,
	vector<double> &n_unitigs,
	vector<double> &e_size
	)
{
	uint64_t total_kmers = reads.size() * (reads[0].length() - k + 1);
	uint64_t n_reads = reads.size();

	// RANDOM GENERATOR
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::uniform_int_distribution<uint64_t> uniform_read_distribution(0,n_reads - 1);

	vector<double> n_sampled_nodes_weighted(max_abundance + 1, 0);

	vector<double> n_internal_local(max_abundance + 1, 0);
	vector<double> n_starts_local(max_abundance + 1, 0);
	vector<double> e_size_sum_length(max_abundance + 1, 0);
	vector<double> e_size_sum_length_squared(max_abundance + 1, 0);

	uint64_t n_sampled_kmers = 0;
	vector<uint64_t> n_sampled_nodes(max_abundance + 1, 0);
	vector<uint64_t> n_sampled_unitigs(max_abundance + 1, 0);
	bool sampled_enough_kmers = false;
	bool sampled_enough_unitigs = false;
	
	int a_w_max_ss_kmers = -1;
	uint64_t max_sample_size_kmers = 0;
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		if (sample_size_kmers[a] >= max_sample_size_kmers)
		{
			max_sample_size_kmers = sample_size_kmers[a];
			a_w_max_ss_kmers = a;
		}
	}

	int min_abundance_unitigs = min_abundance;
	uint64_t n_total_sampled_unitigs = 0;
	// omp_set_dynamic(0);
	#pragma omp parallel for shared(sampled_enough_kmers,sampled_enough_unitigs)num_threads(N_THREADS)
	for (uint64_t i = 0; i < 10 * n_reads; i++)
	{
		//#pragma omp flush (sampled_enough_kmers)
		if ((not sampled_enough_kmers) or (not sampled_enough_unitigs))
		{
			if (n_sampled_kmers >= sample_size_kmers[a_w_max_ss_kmers])
			{
				sampled_enough_kmers = true;
				//continue;
			}
			
			uint64_t read_index = uniform_read_distribution(generator); // (rand() / (double)RAND_MAX) * reads.size();
			string read = reads[read_index];
			// // MAKE SURE THIS IS OK!
			// if (rand() / (double)RAND_MAX < 0.5)
			// {
			// 	read = reverse_complement(read);
			// }
			std::uniform_int_distribution<uint64_t> uniform_pos_distribution(0,read.length() - k);
			int pos = uniform_pos_distribution(generator); // rand() % (read.length() - k + 1);
	        string sample = read.substr(pos,k);

	        #pragma omp atomic
        	n_sampled_kmers++;	
	        
	        // OPTIMIZE THIS IF POSSIBLE:
        	if (sample.find('N') != string::npos)
        	{
        		continue;
        	}

        	int sample_abundance = calc_abundance(rlcsa, sample);
        	assert(sample_abundance > 0);
        	double sample_weight = 1 / (double)sample_abundance;
        	int for_limit = MIN(max_abundance,sample_abundance);

        	// computing E-size estimates
			vector< vector<uint64_t> > u_length(max_abundance + 1);
    		if (not sampled_enough_unitigs)
    		{	
    			get_unitig_stats_SMART(sample, sample_abundance, rlcsa, min_abundance_unitigs, for_limit, u_length);
				
    			for (int a = min_abundance_unitigs; a <= for_limit; a++)
    			{
    				for (auto length : u_length[a])
    				{
    					double length_to_add = (double)(length + k - 1) / sample_abundance;
    					#pragma omp critical
    					{
    						n_sampled_unitigs[a]++;
    						n_total_sampled_unitigs++;
    						// cout << "n_total_sampled_unitigs = " << n_total_sampled_unitigs << endl;
    						// cout << "n_sampled_unitigs[" << a << "]=" << n_sampled_unitigs[a] << endl;
    						e_size_sum_length[a] += length_to_add;
    						e_size_sum_length_squared[a] += (length + k - 1) * length_to_add;        				
    					}	
    				}
    			}
    			bool thread_sampled_enough_unitigs = true;	
    			for (int a = min_abundance_unitigs; a <= max_abundance; a++)
    			{
    				if (n_sampled_unitigs[a] < sample_size_unitigs[a])
					{
						thread_sampled_enough_unitigs = false;
					} 
					else
					{
						#pragma omp critical
						{
							// PROBABLY THIS IS STILL NOT THREAD SAFE!
							//cout << "abundance " << a << " is sampled for unitigs" << endl;
							min_abundance_unitigs++;	
						}
					}
    			}	
    			if (thread_sampled_enough_unitigs)
    			{
    				#pragma omp critical
    				{
    					sampled_enough_unitigs = true;		
    				}	
    			}
    		}

    		// computing the other estimates
        	vector<int> in_degree(max_abundance + 1), out_degree(max_abundance + 1);
        	get_in_out_degrees(sample,rlcsa,min_abundance,max_abundance,in_degree,out_degree);

        	for (int a = min_abundance; a <= for_limit; a++)
        	{
        		#pragma omp critical
        		{
        			n_sampled_nodes_weighted[a] += 1 * sample_weight;
        		}
        		// is internal
    			if ((out_degree[a] == 1) and (in_degree[a] == 1)) 
            	{
            		#pragma omp critical
            		{
            			n_internal_local[a] += 1 * sample_weight;
            		}
            	} else
            	// is start of some unitigs and not isolated
				if ((out_degree[a] > 1) or ((out_degree[a] == 1) and (in_degree[a] != 1))) 
            	{
            		#pragma omp critical
            		{
            			n_starts_local[a] += 1 * sample_weight * out_degree[a];	
            		}

            	} else
            	// is isolated node
            	if ((out_degree[a] == 0) and (in_degree[a] == 0)) 
            	{
            		#pragma omp critical
            		{
            			n_starts_local[a] += 1 * sample_weight;	
            		}
            	}
        	}
    	} 
	}

	cout << "Sampled " << n_total_sampled_unitigs << " unitigs " << endl;

	if (not sampled_enough_kmers)
	{
		cout << "*** I could sample only " << n_sampled_kmers << " kmers out of " <<  sample_size_kmers[a_w_max_ss_kmers] << endl;
		assert(false);
	}

	for (int a = min_abundance; a <= max_abundance; a++)
	{
		n_internal[a] = n_internal_local[a];
		n_starts[a] = n_starts_local[a];
		
		assert(n_sampled_kmers > 0);
		n_nodes[a] = total_kmers / n_sampled_kmers * n_sampled_nodes_weighted[a];
		n_unitigs[a] = total_kmers / n_sampled_kmers * n_starts[a];
		e_size[a] = e_size_sum_length_squared[a] / e_size_sum_length[a];
	}
}

// uint64_t get_sample_size(
// 	const double &prop_external_k, 
// 	const double &delta_avg_unitig_length)
// {
// 	double p = MIN(prop_external_k,0.5);
// 	//double delta_max = delta_avg_unitig_length*( p ) / (double)( 1 + (1 - p) * delta_avg_unitig_length );
//     //double delta_max = delta_avg_unitig_length / (double)(2 + delta_avg_unitig_length);
//     double delta_max = delta_avg_unitig_length;
//     double delta_p_external_k_plus_one = (double)(p * delta_max);

//     assert(delta_p_external_k_plus_one > 0);
//     return pow(TWOSIDED95P_QUANTILE / delta_p_external_k_plus_one, 2) * (1 - p) * p;
// }

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

	parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("input fastq file (all reads need to have the same length)");
	parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output file");
	parser.add_option("-a", "--minabundance") .type("int") .dest("a") .action("store") .set_default(3) .help("try all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("int") .dest("A") .action("store") .set_default(3) .help("try all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(0) .help("number of threads; use 0 for all cores (default: %default)");
	parser.add_option("-k", "--mink") .type("int") .dest("k") .action("store") .set_default(15) .help("minimum kmer size to try (default: %default)");
	parser.add_option("-K", "--maxk") .type("int") .dest("K") .action("store") .set_default(0) .help("maximum kmer sizeto try (default: read_length - 10)");
	parser.add_option("-e", "--relerror") .type("float") .dest("e") .action("store") .set_default(0.1) .help("relative error of the estimations (default: %default)");
	parser.add_option("-b", "--buildindex") .action("store_true") .dest("buildindex") .help("force the index to be rebuilt, even though it exists");
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
		cout << "*** Running on " << N_THREADS << " cores (change this number with option -t)" << endl;
	}
	relative_error = (double) options.get("e");


	// if the index does not exist we need to build it
	if ((not is_readable(indexFileName + ".rlcsa.array")) or (not is_readable(indexFileName + ".rlcsa.parameters")) or buildindex) 
	{
		cout << "*** Building the RLCSA index on the reads" << endl;
		if (EXIT_FAILURE == get_data_for_rlcsa(readFileName, data, char_count))
		{
			return EXIT_FAILURE;
		}
		// Build RLCSA and report some information.
		RLCSA rlcsa_built(data, char_count, 32, 0, N_THREADS, true);
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
	cout << "*** Loading the RLCSA index for the reads from file '" << indexFileName << "' (force the index to be re-built with option -b)" << endl;
 	const RLCSA* rlcsa = new RLCSA(indexFileName, false);
 	if (!(rlcsa->isOk())) 
 	{
 		return EXIT_FAILURE;
 	}
 	// rlcsa->printInfo();
 	// rlcsa->reportSize(true);

 	// we load the reads
 	if (EXIT_FAILURE == get_reads(readFileName, reads))
	{
		return EXIT_FAILURE;
	}

 	// these vectors get re-written for each value of k
 	vector<uint64_t> sample_size_kmers(max_abundance + 1, 0);
 	vector<uint64_t> sample_size_unitigs(max_abundance + 1, MAX_SAMPLE_SIZE_UNITIGS);
 	vector<double> n_internal(max_abundance + 1,1), n_starts(max_abundance + 1,1);
 	vector<double> n_nodes(max_abundance + 1,1), n_unitigs(max_abundance + 1,1);
 	vector<double> e_size(max_abundance + 1,1);

 	cout << "*** Writing results to files:" << endl;
    vector<ofstream> outputFile(max_abundance + 1);
    for (int a = min_abundance; a <= max_abundance; a++)
    {
    	cout << "***    " << outputFileName + ".mink" + int_to_string(mink) + ".maxk" + int_to_string(maxk) + ".a" + int_to_string(a) + ".metrics.csv" << endl;
    	outputFile[a].open((outputFileName + ".mink" + int_to_string(mink) + ".maxk" + int_to_string(maxk) + ".a" + int_to_string(a) + ".metrics.csv").c_str());
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    } 

    if (maxk == 0)
    {
    	maxk = reads[0].length() - 10;
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
 			sample_size_kmers[a] = max_sample_size;
 			if (k == mink)
 			{
				sample_size_unitigs[a] = MAX_SAMPLE_SIZE_UNITIGS;
 			} 
 			else
 			{
 				sample_size_unitigs[a] = MIN(RATIO_SAMPLE_SIZE_UNITIGS * n_nodes[a], MAX_SAMPLE_SIZE_UNITIGS);
 				sample_size_unitigs[a] = MAX(sample_size_unitigs[a],MIN_SAMPLE_SIZE_UNITIGS);	
 			}
 			
 		}

 		// sampling
 		sample_nodes(rlcsa, k, min_abundance, max_abundance, reads, sample_size_kmers, sample_size_unitigs, n_internal, n_starts, n_nodes, n_unitigs, e_size);	

 		// printing the results
 		for (int a = min_abundance; a <= max_abundance; a++)
 		{
 			assert(n_starts[a] > 0);
			double avg_nodes_unitig = n_internal[a] / (double)n_starts[a];

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << (uint64_t)n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << ".,"; // number of edges
	 		outputFile[a] << avg_nodes_unitig << ","; // average number of internal nodes in unitigs
			outputFile[a] << avg_nodes_unitig + k + 1 << ","; // average length of unitigs
			outputFile[a] << sample_size_kmers[a] << ","; // estimated sample size for kmers
			outputFile[a] << n_unitigs[a] << ","; // number of unitigs
			outputFile[a] << e_size[a]; // e-size
			outputFile[a] << endl; 
			//outputFile[a].flush();

	 		cout << k << " " << a << " avg internal nodes=" << (uint64_t)avg_nodes_unitig << " avg length=" << (uint64_t)avg_nodes_unitig + k + 1 << " n_nodes=" << (uint64_t)n_nodes[a] << " n_unitigs=" << (uint64_t)n_unitigs[a] << " e_size=" << e_size[a] << " ess=" << (uint64_t)sample_size_kmers[a] << endl;
	 		cout.flush();
 		}
 	}
 	
 	for (int a = min_abundance; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}

	return EXIT_SUCCESS;
}