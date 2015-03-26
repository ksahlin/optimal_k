
#include "utils.h"
#include "esize_estimation.h"

#define TWOSIDED95P_QUANTILE 1.96
#define MIN_SAMPLE_SIZE_UNITIGS 500
#define MAX_SAMPLE_SIZE_UNITIGS 5000
#define RATIO_SAMPLE_SIZE_UNITIGS 0.05
#define ESIZE_UPDATE_STATS_STEP 5000
#define ESIZE_MAX_SAMPLED_UNITIGS 100000

uint32_t N_THREADS;

inline void get_in_out_degrees(const string& node, 
	const RLCSA* rlcsa, 
	const uint32_t &min_abundance,
	const uint32_t &max_abundance,
	vector<uint32_t> &in_degree,
	vector<uint32_t> &out_degree
	)
{
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		in_degree[a] = 0;
		out_degree[a] = 0;	
	}
	
	string neighbor;
	pair_type result;
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = nucl + node.substr(0,node.length()-1);
		uint32_t neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (uint32_t a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		in_degree[a]++;
	 	}
	}

	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = node.substr(1,node.length()-1) + nucl;
		uint32_t neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (uint32_t a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		out_degree[a]++;
	 	}
	}
}


void sample_nodes(const RLCSA* rlcsa, 
	const uint32_t k,
	const uint32_t min_abundance,
	const uint32_t max_abundance,
	const vector<compact_read_t>& reads,
	const uint64_t &reads_total_content,
	const uint64_t &reads_number,
	const uint64_t &reads_max_length,
	const vector<uint64_t>& sample_size_start_or_internal_nodes,
	const vector<uint64_t>& sample_size_unitigs,
	vector<double> &n_internal,
	vector<double> &n_starts,
	vector<double> &n_nodes,
	vector<double> &n_unitigs,
	vector<double> &avg_unitig_length,
	vector<double> &avg_unitig_length_error,
	vector<double> &e_size,
	vector<double> &e_size_error,
	double relative_error
	)
{
	uint64_t total_kmers = reads_total_content - (k - 1) * reads_number;
	//uint64_t n_reads = reads_number;
	uint64_t n_sampled_reads = reads.size();

	// initializing the RANDOM GENERATOR
	random_device rd;
	default_random_engine generator(rd());
	//////
	uniform_int_distribution<uint64_t> uniform_read_distribution(0,n_sampled_reads - 1);
	//////
	vector< uniform_int_distribution<int> > uniform_pos_distribution(reads_max_length + 1);
	for (uint32_t read_length = k; read_length <= reads_max_length; read_length++)
	{
		uniform_pos_distribution[read_length] = uniform_int_distribution<int>(0,read_length - k);
	}

	// initializing the vectors needed for counting
	vector<double> n_sampled_nodes_weighted(max_abundance + 1, 0);
	vector<double> n_internal_local(max_abundance + 1, 0);
	vector<double> n_starts_local(max_abundance + 1, 0);
	vector<double> e_size_sum_length(max_abundance + 1, 0);
	vector<double> e_size_sum_length_squared(max_abundance + 1, 0);

	vector<double> e_size_e_x(max_abundance + 1, 0);
	vector<double> e_size_e_x2(max_abundance + 1, 0);
	vector<double> e_size_var_x(max_abundance + 1, 0);
	vector<double> e_size_var_x2(max_abundance + 1, 0);
	vector<double> e_size_cov_x2_x(max_abundance + 1, 0);
	vector< vector<unitig_t> > sampled_unitigs(max_abundance + 1);
	unordered_map<string, vector< vector<uint64_t> > > stored_sampled_unitigs;
	unordered_set<uint32_t> e_size_alive_abundances;
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		e_size_alive_abundances.insert(a);	
	}


	// initializing the variables controlling the sample size
	uint64_t n_sampled_kmers = 0;
	vector<uint64_t> n_sampled_start_or_internal_nodes(max_abundance + 1, 0);
	vector<uint64_t> n_sampled_unitigs(max_abundance + 1, 0);
	bool sampled_enough_start_or_internal_nodes = false;
	bool sampled_enough_unitigs = false;
	
	uint32_t a_w_max_ss_start_or_internal_nodes = -1;
	uint64_t max_sample_size_start_or_internal_nodes = 0;
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		if (sample_size_start_or_internal_nodes[a] >= max_sample_size_start_or_internal_nodes)
		{
			max_sample_size_start_or_internal_nodes = sample_size_start_or_internal_nodes[a];
			a_w_max_ss_start_or_internal_nodes = a;
		}
	}

	uint64_t n_total_sampled_unitigs = 0;
	uint64_t n_total_sampled_unitigs_reset = 0;
	// omp_set_dynamic(0);
	// shared(sampled_enough_start_or_internal_nodes,sampled_enough_unitigs)
	#pragma omp parallel for num_threads(N_THREADS)
	for (uint64_t i = 0; i < 100 * n_sampled_reads; i++)
	{
		if ((not sampled_enough_start_or_internal_nodes) or (not sampled_enough_unitigs))
		{
			if (n_sampled_start_or_internal_nodes[a_w_max_ss_start_or_internal_nodes] >= sample_size_start_or_internal_nodes[a_w_max_ss_start_or_internal_nodes])
			{
				#pragma omp critical
				sampled_enough_start_or_internal_nodes = true;
			}
			uint64_t read_index;
			uint32_t pos;
			string sample;

			#pragma omp critical
			{
				read_index = uniform_read_distribution(generator);	
			}
			compact_read_t cread = reads[read_index];
			if (cread.length < k)
			{
				continue;
			}
			#pragma omp critical
			{
				pos = uniform_pos_distribution[cread.length](generator);	
			}
	       	sample = decode_substring(cread,pos,k);

   			// // MAKE SURE THIS IS OK!
			// if (rand() / (double)RAND_MAX < 0.5)
			// {
			// 		sample = reverse_complement(sample);
			// }

	        #pragma omp atomic
        	n_sampled_kmers++;	
	        
	        // OPTIMIZE THIS IF POSSIBLE:
        	if (sample.find('N') != string::npos)
        	{
        		continue;
        	}

        	uint32_t sample_abundance = calc_abundance(rlcsa, sample);
        	//assert(sample_abundance > 0);
        	double sample_weight = 1 / (double)sample_abundance;
        	uint32_t for_limit = MIN(max_abundance,sample_abundance);

    		if (not sampled_enough_unitigs)
    		{	
	        	// computing E-size estimates
				vector< vector<uint64_t> > u_length(max_abundance + 1);
				if (stored_sampled_unitigs.count(sample) == 0)
				{
					get_unitig_stats_SMART(sample, sample_abundance, rlcsa, e_size_alive_abundances, for_limit, u_length);	
					std::pair<string, vector< vector<uint64_t> > > new_sample (sample,u_length);
					#pragma omp critical
					stored_sampled_unitigs.insert(new_sample);
				}
				else
				{
					u_length = stored_sampled_unitigs[sample];
				}
				
    			for (auto a : e_size_alive_abundances)
    			{
    				if (a > for_limit)
    				{
    					break;
    				}
    				#pragma omp atomic
    				n_sampled_unitigs[a] += u_length[a].size();
    				#pragma omp atomic
    				n_total_sampled_unitigs += u_length[a].size();
    				#pragma omp atomic
    				n_total_sampled_unitigs_reset += u_length[a].size();
    				for (auto length : u_length[a])
    				{
    					unitig_t new_unitig;
    					new_unitig.length = length + k - 1;
    					new_unitig.abundance = sample_abundance;
    					#pragma omp critical 
    					sampled_unitigs[a].push_back(new_unitig);
    				}
    			}
    			//cout << "n_total_sampled_unitigs = " << n_total_sampled_unitigs << endl;
    			// only Master gets to update the unitig stats
    			if (omp_get_thread_num() == 0)
    			{
    				if (n_total_sampled_unitigs_reset > ESIZE_UPDATE_STATS_STEP)
	    			{
	    				bool sampled_enough_unitigs_temp = true;
	    				n_total_sampled_unitigs_reset = 0;

	    				#pragma omp critical
	    				{
	    					for (auto itr = e_size_alive_abundances.begin(); itr != e_size_alive_abundances.end(); )
	    					{
	    						uint32_t a = *itr;
	    						double e_x = 0;
	    						double e_x2 = 0;
	    						double var_x = 0;
	    						double var_x2 = 0;
	    						double cov_x2_x = 0;
	    						double sum_1a = 0;
	    						double sum_1a2 = 0;
	    						for (auto unitig : sampled_unitigs[a])
	    						{
	    							e_x += unitig.length * (double)1 / unitig.abundance;
	    							e_x2 += pow(unitig.length, 2) * (double)1 / unitig.abundance;
	    							sum_1a += (double)1 / unitig.abundance;
	    							sum_1a2 += pow((double)1 / unitig.abundance,2);
	    						}
	    						e_x = e_x * (double)1 / sum_1a;
	    						e_x2 = e_x2 * (double)1 / sum_1a;
	    						for (auto unitig : sampled_unitigs[a])
	    						{
	    							var_x += pow(unitig.length - e_x, 2) * (double)1 / unitig.abundance;
	    							var_x2 += pow(pow(unitig.length,2) - e_x2, 2) * (double)1 / unitig.abundance;
	    							cov_x2_x += (pow(unitig.length,2) - e_x2) * (unitig.length - e_x) * (double)1 / unitig.abundance;
	    						}
	    						var_x = var_x * (double)1 / sum_1a; // biased version
	    						var_x = var_x / (1 - sum_1a2 / pow(sum_1a,2)); // unbiased version
	    						var_x = var_x * sum_1a2 / pow(sum_1a,2); // what we want
	    						
	    						var_x2 = var_x2 * (double)1 / sum_1a; // biased version
	    						var_x2 = var_x2 / (1 - sum_1a2 / pow(sum_1a,2)); // unbiased version
	    						var_x2 = var_x2 * sum_1a2 / pow(sum_1a,2); // what we want

	    						cov_x2_x = cov_x2_x * (double)1 / sum_1a; // biased version
	    						cov_x2_x = cov_x2_x / (1 - sum_1a2 / pow(sum_1a,2)); // unbiased version
	    						cov_x2_x = cov_x2_x * sum_1a2 / pow(sum_1a, 2); // what we want

	    						double esize = e_x2 / e_x; 
	    						double sigma = sqrt(var_x2 / pow(e_x,2) - 2 * e_x2 / pow(e_x,3) * cov_x2_x + pow(e_x2,2) / pow(e_x,4) * var_x);
	    						e_size[a] = esize;
	    						e_size_error[a] = TWOSIDED95P_QUANTILE * sigma / esize;

	    						avg_unitig_length[a] = e_x;
	    						avg_unitig_length_error[a] = TWOSIDED95P_QUANTILE * sqrt(var_x) / e_x;

	    						++itr;
	    						if (e_size_error[a] > relative_error)
	    						{
	    							sampled_enough_unitigs_temp = false;
	    						}
	    						else
	    						{
	    							e_size_alive_abundances.erase(a);
	    							// cout << "-" << a << " ";
	    							// for (auto a2 : e_size_alive_abundances)
	    							// {
	    							// 	cout << a2 << " ";
	    							// }
	    							// cout << endl;
	    						}
	    					}	
	    				}
	    				sampled_enough_unitigs = sampled_enough_unitigs_temp;
	    			}
	    			bool sampled_enough_unitigs_temp = true;
	    			for (auto itr = e_size_alive_abundances.begin(); itr != e_size_alive_abundances.end(); )
	    			{
	    				uint32_t a = *itr;
	    				++itr;
	    				if (n_sampled_unitigs[a] < ESIZE_MAX_SAMPLED_UNITIGS)
	    				{
	    					sampled_enough_unitigs_temp = false;
	    					break;
	    				}
	    				// else
	    				// {
	    				// 	e_size_alive_abundances.erase(a);
	    				// }
	    			}
	    			sampled_enough_unitigs = sampled_enough_unitigs_temp;
    			}
    		}

    		//cout << "new sample " << sample << endl;

    		// computing the other estimates
        	vector<uint32_t> in_degree(max_abundance + 1), out_degree(max_abundance + 1);
        	get_in_out_degrees(sample,rlcsa,min_abundance,max_abundance,in_degree,out_degree);

        	for (uint32_t a = min_abundance; a <= for_limit; a++)
        	{
        		#pragma omp atomic
       			n_sampled_nodes_weighted[a] += sample_weight;
       			#pragma omp atomic
           		n_sampled_start_or_internal_nodes[a]++;

        		// is internal
    			if ((out_degree[a] == 1) and (in_degree[a] == 1)) 
            	{
            		#pragma omp atomic
           			n_internal_local[a] += sample_weight;
            	} 
            	else
            	// is start of some unitigs and not isolated
				if ((out_degree[a] > 1) or ((out_degree[a] == 1) and (in_degree[a] != 1))) 
            	{
            		#pragma omp atomic
            		n_starts_local[a] += sample_weight * out_degree[a];	

            	} else
            	// is isolated node
            	if ((out_degree[a] == 0) and (in_degree[a] == 0)) 
            	{
            		#pragma omp atomic
           			n_starts_local[a] += sample_weight;	
            	}
        	}
    	} 
	}

	cout << "Sampled " << n_total_sampled_unitigs << " unitigs out of which" << endl;
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		cout << n_sampled_unitigs[a] << " for abundance " << a << endl;
	}

	if (not sampled_enough_start_or_internal_nodes)
	{
		cout << "*** I could sample only " << n_sampled_start_or_internal_nodes[a_w_max_ss_start_or_internal_nodes] << " nodes out of " <<  sample_size_start_or_internal_nodes[a_w_max_ss_start_or_internal_nodes] << endl;
		assert(false);
	}

	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		n_internal[a] = n_internal_local[a];
		n_starts[a] = n_starts_local[a];
		
		assert(n_sampled_kmers > 0);
		n_nodes[a] = total_kmers / n_sampled_kmers * n_sampled_nodes_weighted[a];
		n_unitigs[a] = total_kmers / n_sampled_kmers * n_starts[a];
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
    uint32_t mink, maxk;
    string readFileName, outputFileName;
    string buildindex, loadindex;
	uint32_t min_abundance,max_abundance;
	double relative_error; // maximum relative error 10% of our estimators
	bool lowermemory;

	std::cout << std::fixed;
    std::cout << std::setprecision(2);

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

	parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("a file containing a list of FASTA/Q(.gz) file names, one per line");
	parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output file");
	parser.add_option("-a", "--minabundance") .type("uint32_t") .dest("a") .action("store") .set_default(1) .help("try all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("uint32_t") .dest("A") .action("store") .set_default(5) .help("try all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("uint32_t") .dest("t") .action("store") .set_default(0) .help("number of threads; use 0 for all cores (default: %default)");
	parser.add_option("-k", "--mink") .type("uint32_t") .dest("k") .action("store") .set_default(15) .help("try all kmer sizes starting with this value (default: %default)");
	parser.add_option("-K", "--maxk") .type("uint32_t") .dest("K") .action("store") .set_default(0) .help("try all kmer sizes up to this value (default: read_length - 10)");
	parser.add_option("-e", "--relerror") .type("float") .dest("e") .action("store") .set_default(0.1) .help("relative error of the estimations (default: %default)");
	parser.add_option("-b", "--buildindex"). type("string") .dest("buildindex") .action("store") .set_default("") .help("the filename where the index should be saved");
	parser.add_option("-l", "--loadindex"). type("string") .dest("loadindex") .action("store") .set_default("") .help("the filename from where the index should be loaded");
	parser.add_option("-m", "--lowermemory") .dest("lowermemory") .action("store_true") .set_default(false) .help("force the index construction to use less RAM; this slows the construction");
	optparse::Values& options = parser.parse_args(argc, argv);

	lowermemory = (options.get("lowermemory") ? true : false);
	buildindex = (string) options.get("buildindex");
	loadindex = (string) options.get("loadindex");
	readFileName = (string) options.get("r");
	outputFileName = (string) options.get("o");
	min_abundance = (uint32_t) options.get("a");
	max_abundance = (uint32_t) options.get("A");
	mink = (uint32_t) options.get("k");
	maxk = (uint32_t) options.get("K");
	N_THREADS = (uint32_t) options.get("t");
	if (N_THREADS == 0)
	{
		N_THREADS = omp_get_num_procs();
		cout << "*** Running on " << N_THREADS << " cores (change this number with option -t)" << endl;
	}
	relative_error = (double) options.get("e");

	if (readFileName == "")
	{
		cerr << "*** ERROR: -r <read_file> must be given" << endl;	
		return EXIT_FAILURE;
	}
	// checking if we have the index
	if ((buildindex == "") and (loadindex == ""))
	{
		cerr << "*** ERROR: either --buildindex <file_name> or --loadindex <file_name> must be given" << endl;
		return EXIT_FAILURE;
	}
	// checking if we have the index
	if ((buildindex != "") and (loadindex != ""))
	{
		cerr << "*** ERROR: not both of --buildindex <file_name> or --loadindex <file_name> must be given" << endl;
		return EXIT_FAILURE;
	}

	// compact_read_t cread = encode_string("ACTTGGTACAT");
	// for (uint32_t pos = 0; pos < 11; pos++)
	// {
	// 	cout << "pos = " << pos << endl;
	// 	for (uint32_t length = 0; length <= 11; length++)
	// 	{
	// 		cout << length << ": " << decode_substring(cread,pos,length) << endl;
	// 	}
	// }
		
	// cout << decode_substring(cread,4,5) << endl;
	// cout << decode_substring(cread,4,6) << endl;
	// cout << decode_substring(cread,0,5) << endl;
	// cout << decode_substring(cread,0,11) << endl;
	// cout << decode_substring(cread,1,6) << endl;
	// cout << decode_substring(cread,0,8) << endl;
	// cout << decode_substring(cread,0,9) << endl;
	// cout << decode_substring(cread,0,11) << endl;
	// // cout << decode_substring(cread,3,11) << endl;
	// return EXIT_SUCCESS;

	// if we need to build the index
	if (buildindex != "")
	{
		get_data_and_build_rlcsa_iterative(readFileName, buildindex, N_THREADS, lowermemory);
		// get_data_and_build_rlcsa_noniterative(readFileName, buildindex, N_THREADS);
		cout << "*** SUCCESS: now run the program with --loadindex " << buildindex << endl;
		return EXIT_SUCCESS;		
	}

	if (outputFileName == "")
	{
		cerr << "*** ERROR: The argument -o|--outputfile is needed" << endl;
		return EXIT_FAILURE;
	}

	if (lowermemory)
	{
		cerr << "*** Ignoring option --lowermemory because it works only when constructing the index" << endl;
	}
	// we need to load the index
	cout << "*** Loading the RLCSA index for the reads from files: (force the index to be re-built with option -b|--buildindex)" << endl;
	cout << "***    " << loadindex + ".rlcsa.array" << endl;
	cout << "***    " << loadindex + ".rlcsa.parameters" << endl;
	

 	const RLCSA* rlcsa = new RLCSA(loadindex, false);
 	if (!(rlcsa->isOk())) 
 	{
 		return EXIT_FAILURE;
 	}
 	// rlcsa->printInfo();
 	// rlcsa->reportSize(true);
 	cout << "*** Loaded the RLCSA index " << endl;

 	vector<compact_read_t> reads;
 	uint64_t reads_total_content, reads_number;
 	uint32_t reads_max_length, reads_min_length;
 	// we load the reads
 	if (EXIT_FAILURE == get_reads(readFileName, reads, reads_total_content, reads_number, reads_max_length, reads_min_length))
	{
		return EXIT_FAILURE;
	}

 	// these vectors get re-written for each value of k
 	vector<uint64_t> sample_size_start_or_internal_nodes(max_abundance + 1, 0);
 	vector<uint64_t> sample_size_unitigs(max_abundance + 1, MAX_SAMPLE_SIZE_UNITIGS);
 	vector<double> n_internal(max_abundance + 1,1), n_starts(max_abundance + 1,1);

 	uint64_t total_kmers_for_mink = reads_total_content - (mink - 1) * reads_number;
 	vector<double> n_nodes(max_abundance + 1, total_kmers_for_mink / 2), n_unitigs(max_abundance + 1, total_kmers_for_mink / 2);

    if (maxk == 0)
    {
    	// maxk = reads_max_length - 10;
    	// maxk = reads_min_length - 10;
    	maxk = (reads_total_content / reads_number) - 10; // avg read length - 10
    	cout << "*** Setting the maximum kmer size to " << maxk << " (average read length - 10)" << endl;
    }

 	cout << "*** Writing results to files:" << endl;
    vector<ofstream> outputFile(max_abundance + 1);
    for (uint32_t a = min_abundance; a <= max_abundance; a++)
    {
    	cout << "***    " << outputFileName + ".a" + int_to_string(a) + ".csv" << endl;
    	outputFile[a].open((outputFileName + ".a" + int_to_string(a) + ".csv").c_str());
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    	outputFile[a] << std::fixed;
    	outputFile[a] << std::setprecision(2);
    } 

 	for (uint32_t k = mink; k <= maxk; k++)
 	{
 		vector<double> e_size(max_abundance + 1,1);
 		vector<double> e_size_error(max_abundance + 1,1);
 		vector<double> avg_unitig_length(max_abundance + 1,1);
 		vector<double> avg_unitig_length_error(max_abundance + 1,1);

 		// getting the sample size
 		for (uint32_t a = min_abundance; a <= max_abundance; a++)
 		{
 			uint64_t total_kmers = reads_total_content - (k - 1) * reads_number;
 			uint64_t SS_n_nodes = get_sample_size_for_proportion(n_nodes[a] / (double)(total_kmers), relative_error);
 			uint64_t SS_n_unitigs = get_sample_size_for_proportion(n_unitigs[a] / (double)(total_kmers), relative_error);
 			uint64_t SS_n_start_or_internal_nodes = get_sample_size_for_proportion(n_starts[a] / (double)(n_internal[a] + n_starts[a]), 1 / (double)(1 + relative_error) - 1);
 			
 			// cout << "SS_n_nodes = " << SS_n_nodes << endl;
 			// cout << "SS_n_unitigs = " << SS_n_unitigs << endl;
 			// cout << "ratio = " << n_starts[a] / (double)(n_internal[a] + n_starts[a]) << endl;
 			// cout << "SS_n_start_or_internal_nodes = " << SS_n_start_or_internal_nodes << endl;

 			uint64_t max_sample_size = MAX(SS_n_nodes,SS_n_unitigs);
 			max_sample_size = MAX(max_sample_size,SS_n_start_or_internal_nodes);
 			sample_size_start_or_internal_nodes[a] = max_sample_size;
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
 		sample_nodes(rlcsa, 
 			k, 
 			min_abundance, 
 			max_abundance, 
 			reads, 
 			reads_total_content, 
 			reads_number, 
 			reads_max_length, 
 			sample_size_start_or_internal_nodes, 
 			sample_size_unitigs, 
 			n_internal, 
 			n_starts, 
 			n_nodes, 
 			n_unitigs, 
 			avg_unitig_length,
 			avg_unitig_length_error,
 			e_size,
 			e_size_error,
 			relative_error);	

 		// pruint32_ting the results
 		for (uint32_t a = min_abundance; a <= max_abundance; a++)
 		{
 			assert(n_starts[a] > 0);
			double avg_nodes_unitig = n_internal[a] / (double)n_starts[a];

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << (uint64_t)n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << ".,"; // number of edges
	 		outputFile[a] << avg_nodes_unitig << ","; // average number of internal nodes in unitigs
			outputFile[a] << avg_nodes_unitig + k + 1 << ","; // average length of unitigs
			outputFile[a] << sample_size_start_or_internal_nodes[a] << ","; // estimated sample size for kmers
			outputFile[a] << (uint64_t)n_unitigs[a] << ","; // number of unitigs
			outputFile[a] << e_size[a]; // e-size
			outputFile[a] << endl; 
			//outputFile[a].flush();

	 		cout << k << " " << a << " avg internal nodes=" << (uint64_t)avg_nodes_unitig << " avg length=" << (uint64_t)avg_nodes_unitig + k + 1 << " = " << avg_unitig_length[a] << " rel_err=" << avg_unitig_length_error[a] << " n_nodes=" << (uint64_t)n_nodes[a] << " n_unitigs=" << (uint64_t)n_unitigs[a] << " e_size=" << e_size[a] << " rel_err=" << e_size_error[a] << " ess=" << (uint64_t)sample_size_start_or_internal_nodes[a] << endl;
	 		cout.flush();
 		}
 	}
 	
 	for (uint32_t a = min_abundance; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}

 	delete rlcsa;

	return EXIT_SUCCESS;
}