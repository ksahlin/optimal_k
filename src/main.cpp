
#include "utils.h"
#include "esize_estimation.h"

#define TWOSIDED95P_QUANTILE 1.96
#define KMER_UPDATE_STATS_STEP 5000
#define KMER_MAX_SAMPLE 500000
#define ESIZE_UPDATE_STATS_STEP 9000
// #define ESIZE_MAX_SAMPLED_UNITIGS_NONUNIQUEKMERS 100

#define LARGE_NUMBER 1048576

#define _first_k 15
#define _last_k 60
#define _first_a 3
#define _last_a 5



uint64_t ESIZE_MAX_SAMPLED_UNITIGS;
uint32_t N_THREADS;
double _n_nodes_proportion;

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

inline void update_unitig_stats_from_sample(const RLCSA* rlcsa, 
	const uint32_t &k,
	const string& sample,
	const uint32_t &sample_abundance,
	const uint32_t &min_abundance,
	const uint32_t &max_abundance,
	vector<double> &avg_unitig_length,
	vector<double> &avg_unitig_length_error,
	vector<double> &e_size,
	vector<double> &e_size_error,
	const double &relative_error,
	vector< vector<unitig_t> > &sampled_unitigs,
	vector<uint64_t> &n_unique_sampled_unitigs,
	vector<uint64_t> &n_total_sampled_unitigs_reset,
	unordered_map<string, vector< vector<uint64_t> > > &stored_sampled_unitigs,
	unordered_set<uint32_t> &e_size_alive_abundances,
	unordered_set<uint32_t> &e_size_alive_abundances_thread,
	bool &sampled_enough_unitigs
	)
{
	uint32_t for_limit = MIN(max_abundance,sample_abundance);
	// computing E-size estimates
	vector< vector<uint64_t> > u_length(max_abundance + 1);
	if (stored_sampled_unitigs.count(sample) == 0)
	{
		get_unitig_stats_SMART(sample, sample_abundance, rlcsa, e_size_alive_abundances_thread, for_limit, u_length);	
		std::pair<string, vector< vector<uint64_t> > > new_sample (sample,u_length);
		#pragma omp critical
		stored_sampled_unitigs.insert(new_sample);
		
		// updating the unique kmers stats
		for (auto a : e_size_alive_abundances_thread)
		{
			if (a > for_limit) break;
			#pragma omp atomic 
			n_unique_sampled_unitigs[a] += u_length[a].size();
		}
	}
	else
	{
		u_length = stored_sampled_unitigs[sample];
	}
	
	for (auto a : e_size_alive_abundances_thread)
	{
		if (a > for_limit) break;
		#pragma omp atomic
		n_total_sampled_unitigs_reset[a] += u_length[a].size();
		for (auto length : u_length[a])
		{
			unitig_t new_unitig;
			new_unitig.length = length + k - 1;
			new_unitig.abundance = sample_abundance;
			#pragma omp critical 
			sampled_unitigs[a].push_back(new_unitig);
		}
	}

	// only Master gets to update the unitig stats
	if (omp_get_thread_num() == 0)
	{
		#pragma omp critical
		{
			bool sampled_enough_unitigs_temp = true;
			for (auto itr = e_size_alive_abundances_thread.begin(); itr != e_size_alive_abundances_thread.end(); )
			{
				uint32_t a = *itr;
				++itr;

				if (n_total_sampled_unitigs_reset[a] > ESIZE_UPDATE_STATS_STEP)
				{
					n_total_sampled_unitigs_reset[a] = 0;
					
					
					double e_x = 0;
					double e_x2 = 0;
					double var_x = 0;
					double var_x2 = 0;
					double cov_x2_x = 0;
					double A1 = 0;
					double A2 = 0;
					for (auto &unitig : sampled_unitigs[a])
					{
						e_x += unitig.length / (double)unitig.abundance;
						e_x2 += pow(unitig.length, 2) / (double)unitig.abundance;
						A1 += (double)1 / unitig.abundance;
						A2 += (double)1 / pow(unitig.abundance,2);
					}
					e_x = e_x * (double)1 / A1;
					e_x2 = e_x2 * (double)1 / A1;
					for (auto unitig : sampled_unitigs[a])
					{
						var_x += pow(unitig.length - e_x, 2) / (double)unitig.abundance;
						var_x2 += pow(pow(unitig.length,2) - e_x2, 2) / (double)unitig.abundance;
						cov_x2_x += (pow(unitig.length,2) - e_x2) * (unitig.length - e_x) / (double)unitig.abundance;
					}
					var_x = var_x / A1; // biased version
					var_x = var_x / (1 - A2 / pow(A1,2)); // unbiased version
					var_x = var_x * A2 / pow(A1,2); // what we want
					
					var_x2 = var_x2 * (double)1 / A1; // biased version
					var_x2 = var_x2 / (1 - A2 / pow(A1,2)); // unbiased version
					var_x2 = var_x2 * A2 / pow(A1,2); // what we want

					cov_x2_x = cov_x2_x * (double)1 / A1; // biased version
					cov_x2_x = cov_x2_x / (1 - A2 / pow(A1,2)); // unbiased version
					cov_x2_x = cov_x2_x * A2 / pow(A1, 2); // what we want

					double esize = e_x2 / e_x; 
					double sigma = sqrt(var_x2 / pow(e_x,2) - 2 * e_x2 / pow(e_x,3) * cov_x2_x + pow(e_x2,2) / pow(e_x,4) * var_x);
					e_size[a] = esize;
					e_size_error[a] = TWOSIDED95P_QUANTILE * sigma / esize;

					avg_unitig_length[a] = e_x;
					avg_unitig_length_error[a] = TWOSIDED95P_QUANTILE * sqrt(var_x) / e_x;

					if (e_size_error[a] > relative_error)
					{
						sampled_enough_unitigs_temp = false;
					}
					else
					{
						e_size_alive_abundances.erase(a);
						// cout << "-" << a << " REMAINING: ";
						// for (auto a : e_size_alive_abundances)
						// {
						// 	cout << a << " ";
						// }
						// cout << endl;
					}
				}	
			}
			sampled_enough_unitigs = sampled_enough_unitigs_temp;
		}

		bool sampled_enough_unitigs_temp = true;
		#pragma omp critical 
		{
			for (auto itr = e_size_alive_abundances_thread.begin(); itr != e_size_alive_abundances_thread.end(); )
			{
				uint32_t a = *itr;
				++itr;
				if (n_unique_sampled_unitigs[a] < ESIZE_MAX_SAMPLED_UNITIGS)
				{
					sampled_enough_unitigs_temp = false;
					break;
				}
				else
				{
					//e_size_alive_abundances.erase(a);
				}
			}
		}
		sampled_enough_unitigs = sampled_enough_unitigs_temp;
	}
}

inline void update_node_stats_from_sample(const uint32_t &sample_abundance,
	const uint32_t &min_abundance,
	const uint32_t &for_limit,
	const uint32_t &max_abundance,
	const uint64_t &total_kmers,
	const uint64_t &n_sampled_kmers,
	uint64_t &n_sampled_kmers_reset_for_nodes,
	vector<double> &n_nodes_weighted_sum,
	vector<double> &n_nodes,
	vector<double> &n_nodes_error,
	bool &sampled_enough_kmers,
	const double &relative_error
	)
{
	double sample_weight = 1 / (double)sample_abundance;
	for (uint32_t a = min_abundance; a <= for_limit; a++)
	{
		#pragma omp atomic
		n_nodes_weighted_sum[a] += sample_weight;
	}

	// only Master gets to update the node stats
	if (omp_get_thread_num() == 0)
	{
		if (n_sampled_kmers_reset_for_nodes > KMER_UPDATE_STATS_STEP)
		{
			n_sampled_kmers_reset_for_nodes = 0;

			for (uint32_t a = min_abundance; a <= max_abundance; a++)
			{
				double prop_n_nodes_a = n_nodes_weighted_sum[a] / n_sampled_kmers;
				n_nodes[a] = prop_n_nodes_a * total_kmers;
				n_nodes_error[a] = TWOSIDED95P_QUANTILE * sqrt(prop_n_nodes_a * (1 - prop_n_nodes_a) / n_sampled_kmers) / prop_n_nodes_a;
			}

			if (not sampled_enough_kmers)
			{
				bool sampled_enough_kmers_temp = true;
				for (uint32_t a = min_abundance; a <= max_abundance; a++)
				{
					if ((n_sampled_kmers < KMER_MAX_SAMPLE) and (n_nodes_error[a] > relative_error))
					{
						sampled_enough_kmers_temp = false;
						break;
					}
				}
				sampled_enough_kmers = sampled_enough_kmers_temp;	
			}
		}
	}
}

inline void update_n_unitig_stats_from_sample(const uint32_t &sample_abundance,
	const uint32_t &min_abundance,
	const uint32_t &for_limit,
	const uint32_t &max_abundance,
	const vector<uint32_t> &in_degree,
	const vector<uint32_t> &out_degree,
	const uint64_t &total_kmers,
	const uint64_t &n_sampled_kmers,
	uint64_t &n_sampled_kmers_reset_for_unitigs,
	vector<double> &n_unitigs_weighted_sum,
	vector<double> &n_unitigs,
	vector<double> &n_unitigs_error,
	bool &sampled_enough_for_n_unitigs,
	const double &relative_error
	)
{
	double sample_weight = 1 / (double)sample_abundance;
	for (uint32_t a = min_abundance; a <= for_limit; a++)
	{
		// is source of some proper unitig
		if ((out_degree[a] > 1) or ((out_degree[a] == 1) and (in_degree[a] != 1))) 
    	{
    		#pragma omp atomic
    		n_unitigs_weighted_sum[a] += sample_weight * out_degree[a];	

    	} else
    	// is isolated node
    	if ((out_degree[a] == 0) and (in_degree[a] == 0)) 
    	{
    		#pragma omp atomic
   			n_unitigs_weighted_sum[a] += sample_weight;	
    	}
	}

	// only Master gets to update the n_unitig stats
	if (omp_get_thread_num() == 0)
	{
		if (n_sampled_kmers_reset_for_unitigs > KMER_UPDATE_STATS_STEP)
		{
			n_sampled_kmers_reset_for_unitigs = 0;

			for (uint32_t a = min_abundance; a <= max_abundance; a++)
			{
				double prop_n_unitigs_a = n_unitigs_weighted_sum[a] / n_sampled_kmers;
				n_unitigs[a] = prop_n_unitigs_a * total_kmers;
				n_unitigs_error[a] =  TWOSIDED95P_QUANTILE * sqrt(prop_n_unitigs_a * (1 - prop_n_unitigs_a) / n_sampled_kmers) / prop_n_unitigs_a;
			}

			if (not sampled_enough_for_n_unitigs)
			{
				bool sampled_enough_for_n_unitigs_temp = true;
				for (uint32_t a = min_abundance; a <= max_abundance; a++)
				{
					if ((n_sampled_kmers < KMER_MAX_SAMPLE) and (n_unitigs_error[a] > relative_error))
					{
						sampled_enough_for_n_unitigs_temp = false;
						break;
					}
				}
				sampled_enough_for_n_unitigs = sampled_enough_for_n_unitigs_temp;	
			}
		}
	}
}

void sample_nodes(const RLCSA* rlcsa, 
	const uint32_t &k,
	uint32_t min_abundance,
	uint32_t max_abundance,
	const vector<compact_read_t>& reads,
	const uint64_t &reads_total_content,
	const uint64_t &reads_number,
	const uint64_t &reads_max_length,
	vector<double> &n_nodes,
	vector<double> &n_nodes_error,
	vector<double> &n_unitigs,
	vector<double> &n_unitigs_error,
	vector<double> &avg_unitig_length,
	vector<double> &avg_unitig_length_error,
	vector<double> &e_size,
	vector<double> &e_size_error,
	double relative_error,
	const uint64_t &n_nodes_h,
	const bool &verbose
	)
{

	uint64_t total_kmers = reads_total_content - (k - 1) * reads_number;
	//uint64_t n_reads = reads_number;
	uint64_t n_sampled_reads = reads.size();

	// initializing the RANDOM GENERATOR
	random_device rd;
	default_random_engine generator(rd());
	/*here*/ uniform_int_distribution<uint64_t> uniform_read_distribution(0,n_sampled_reads - 1);
	vector< uniform_int_distribution<int> > uniform_pos_distribution(reads_max_length + 1);
	for (uint32_t read_length = k; read_length <= reads_max_length; read_length++)
	{
		uniform_pos_distribution[read_length] = uniform_int_distribution<int>(0,read_length - k);
	}

	/////////// ESIZE_UNITIG AND AVG_UNITIG_LENGTH STATS ////////////
	vector< vector<unitig_t> > sampled_unitigs(max_abundance + 1);
	vector<uint64_t> n_unique_sampled_unitigs(max_abundance + 1, 0);
	unordered_map<string, vector< vector<uint64_t> > > stored_sampled_unitigs;
	unordered_set<uint32_t> e_size_alive_abundances;
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		e_size_alive_abundances.insert(a);
	}
	vector<uint64_t> n_total_sampled_unitigs_reset(max_abundance + 1, 0);
	bool sampled_enough_unitigs = false;

	/////////// NODES STATS ////////////
	uint64_t n_sampled_kmers = 0;
	uint64_t n_sampled_kmers_reset_for_nodes = 0;
	bool sampled_enough_kmers = false;
	vector<double> n_nodes_weighted_sum(max_abundance + 1, 0);

	/////////// N_UNITIGS STATS ////////////
	bool sampled_enough_for_n_unitigs = false;
	uint64_t n_sampled_kmers_reset_for_unitigs = 0;
	vector<double> n_unitigs_weighted_sum(max_abundance + 1, 0);

	/////////// SOME STORED VALUES /////////
	unordered_map<string, uint32_t> stored_abundance;
	unordered_map<string, vector<uint32_t> > stored_in_degree;
	unordered_map<string, vector<uint32_t> > stored_out_degree;

	// heuristic for pruning the search space on (k,a)
	double relative_error_n_unitigs = relative_error;
	if (n_nodes_h == 0)
	{
		sampled_enough_unitigs = true;
		relative_error_n_unitigs = LARGE_NUMBER;
	}

		
	#pragma omp parallel for num_threads(N_THREADS)
	for (uint64_t i = 0; i < 100 * n_sampled_reads; i++)
	{
		if ((not sampled_enough_unitigs) or (not sampled_enough_for_n_unitigs) or (not sampled_enough_kmers))
		{
			unordered_set<uint32_t> e_size_alive_abundances_thread = e_size_alive_abundances;
			uint64_t read_index;
			uint32_t pos;
			string sample;

			#pragma omp critical
			{
				read_index = uniform_read_distribution(generator);	
			}
			compact_read_t cread = reads[read_index];
			if (cread.length < k) continue;

			#pragma omp critical
			{
				pos = uniform_pos_distribution[cread.length](generator);	
			}
	       	sample = decode_substring(cread,pos,k);

	        #pragma omp atomic
        	n_sampled_kmers++;	
        	#pragma omp atomic
        	n_sampled_kmers_reset_for_nodes++;
        	#pragma omp atomic
        	n_sampled_kmers_reset_for_unitigs++;
	        
        	if (sample.find('N') != string::npos) continue;

			uint32_t sample_abundance, for_limit;
			vector<uint32_t> in_degree(max_abundance + 1), out_degree(max_abundance + 1);

        	if (stored_abundance.count(sample) == 0)
        	{
        		sample_abundance = calc_abundance(rlcsa, sample);
        		for_limit = MIN(max_abundance,sample_abundance);
        		get_in_out_degrees(sample,rlcsa,min_abundance,for_limit,in_degree,out_degree);	
        		
        		std::pair<string, uint32_t > new_abundance(sample,sample_abundance);
        		#pragma omp critical
        		stored_abundance.insert(new_abundance);

        		std::pair<string, vector<uint32_t> > new_in_degree(sample,in_degree);
        		#pragma omp critical
        		stored_in_degree.insert(new_in_degree);

        		std::pair<string, vector<uint32_t> > new_out_degree(sample,out_degree);
        		#pragma omp critical
        		stored_out_degree.insert(new_out_degree);
        	}
        	else
        	{
        		sample_abundance = stored_abundance[sample];
        		for_limit = MIN(max_abundance,sample_abundance);
        		in_degree = stored_in_degree[sample];
        		out_degree = stored_out_degree[sample];
        	}

    		if (not sampled_enough_unitigs)
    		{	
    			update_unitig_stats_from_sample(rlcsa, 
    				k,
    				sample,
    				sample_abundance,
    				min_abundance,
    				max_abundance,
    				avg_unitig_length,
    				avg_unitig_length_error,
    				e_size,
    				e_size_error,
    				relative_error,
    				sampled_unitigs,
    				n_unique_sampled_unitigs,
    				n_total_sampled_unitigs_reset,
    				stored_sampled_unitigs,
    				e_size_alive_abundances,
    				e_size_alive_abundances_thread,
    				sampled_enough_unitigs);

    			// heuristic criterion for abandoning this k and a
    			// only master thread gets to update sampled_enough_unitigs
    			if ((omp_get_thread_num() == 0) and (n_nodes_h != 0))
    			{
    				for (uint32_t a = min_abundance; a <= max_abundance; a++)	
    				{
    					// if we have a good estimate of the number of nodes
    					// and the estimated number of nodes is less than _n_nodes_proportion (e.g. 70%) of n_nodes_h
    					// then we abort this abundance and all greater abundances
    					if ((n_nodes_error[a] <= relative_error) and
    						(n_nodes[a] <= _n_nodes_proportion * n_nodes_h))
    					{
    						#pragma omp critical
    						{
    							for (uint32_t a2 = a; a2 <= max_abundance; a2++)
    							{
    								n_unitigs_error[a] = LARGE_NUMBER;
    								e_size_error[a] = LARGE_NUMBER;
    								avg_unitig_length_error[a] = LARGE_NUMBER;
    								e_size_alive_abundances.erase(a2);
    							}	
    												
	    						max_abundance = a - 1;
	    						if (verbose)
	    						{
	    							cout << "Setting max_abundance = " << max_abundance << " for k = " << k << endl;
	    						}
	    						if (max_abundance < min_abundance)
	    						{
	    							sampled_enough_unitigs = true;
	    							sampled_enough_for_n_unitigs = true;
	    						}
    						}
    						break;
    					}
    				}
    			}
    		}

    		update_node_stats_from_sample(sample_abundance,
    			min_abundance,
    			for_limit,
    			max_abundance,
    			total_kmers,
    			n_sampled_kmers,
    			n_sampled_kmers_reset_for_nodes,
    			n_nodes_weighted_sum,
    			n_nodes,
    			n_nodes_error,
    			sampled_enough_kmers,
    			relative_error);

    		update_n_unitig_stats_from_sample(sample_abundance,
    			min_abundance,
    			for_limit,
    			max_abundance,
    			in_degree,
    			out_degree,
    			total_kmers,
				n_sampled_kmers,
				n_sampled_kmers_reset_for_unitigs,
        		n_unitigs_weighted_sum,
    			n_unitigs,
    			n_unitigs_error,
    			sampled_enough_for_n_unitigs,
    			relative_error_n_unitigs);
    	} 
	}

	if (verbose and (n_nodes_h != 0))
	{
		cout << "Sampled:" << endl;
		for (uint32_t a = min_abundance; a <= max_abundance; a++)
		{
			cout << sampled_unitigs[a].size() << " (" << n_unique_sampled_unitigs[a] << " unique) unitigs for abundance " << a << endl;
		}
		cout << "Sampled " << n_sampled_kmers << " k-mers" << endl;	
	}
}

int main(int argc, char** argv)
{
    uint32_t mink, maxk;
    string readFileName, outputFileName;
    string buildindex, loadindex;
	uint32_t min_abundance,max_abundance;
	double relative_error; // maximum relative error
	bool lowermemory, verbose;

	std::cout << std::fixed;
    std::cout << std::setprecision(2);

	// command line argument parser
	string usage = "\n  %prog OPTIONS";
	const string version = "%prog 0.2\nCopyright (C) 2014-2015\n"
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
	parser.add_option("-a", "--minabundance") .type("int") .dest("a") .action("store") .set_default(1) .help("try all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("int") .dest("A") .action("store") .set_default(5) .help("try all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("int") .dest("t") .action("store") .set_default(0) .help("number of threads; use 0 for all cores (default: %default)");
	parser.add_option("-k", "--mink") .type("int") .dest("k") .action("store") .set_default(15) .help("try all kmer sizes starting with this value (default: %default)");
	parser.add_option("-K", "--maxk") .type("int") .dest("K") .action("store") .set_default(0) .help("try all kmer sizes up to this value (default: read_length - 10)");
	parser.add_option("-b", "--buildindex"). type("string") .dest("buildindex") .action("store") .set_default("") .help("the filename where the index should be saved");
	parser.add_option("-l", "--loadindex"). type("string") .dest("loadindex") .action("store") .set_default("") .help("the filename from where the index should be loaded");
	parser.add_option("-m", "--lowermemory") .dest("lowermemory") .action("store_true") .set_default(false) .help("force the index construction to use less RAM; this slows the construction (default %default)");	
	parser.add_option("-e", "--relerror") .type("float") .dest("e") .action("store") .set_default(0.1) .help("relative error of the estimations (default: %default)");
	parser.add_option("-h", "--heuristic") .type("float") .dest("h") .action("store") .set_default(0.90) .help("abandon sampling unitigs for a pair (k,a) if the estimated number of nodes for (k,a) is less than <h> * the estimated number of nodes of the genome (default: %default)");
	parser.add_option("-s", "--samples") .type("int") .dest("s") .action("store") .set_default(15000) .help("maximum number of unique k-mers that are sampled (default: %default)");
	parser.add_option("-v", "--nonverbose") .dest("nonverbose") .action("store_true") .set_default(false) .help("print some stats to stdout (default %default)");

	optparse::Values& options = parser.parse_args(argc, argv);


	lowermemory = (options.get("lowermemory") ? true : false);
	verbose = (options.get("nonverbose") ? false : true);
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
	ESIZE_MAX_SAMPLED_UNITIGS = (uint64_t) options.get("s");
	if (ESIZE_MAX_SAMPLED_UNITIGS < ESIZE_UPDATE_STATS_STEP)
	{
		ESIZE_MAX_SAMPLED_UNITIGS = ESIZE_UPDATE_STATS_STEP;
	}
	relative_error = (double) options.get("e");
	_n_nodes_proportion = (double) options.get("h");

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

	// if we need to build the index
	if (buildindex != "")
	{
		get_data_and_build_rlcsa_iterative(readFileName, buildindex, N_THREADS, lowermemory, verbose);
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
 	if (EXIT_FAILURE == get_reads(readFileName, reads, reads_total_content, reads_number, reads_max_length, reads_min_length, verbose))
	{
		return EXIT_FAILURE;
	}

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
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_internal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size,rel_err n_nodes,rel_err avg_length_unitigs,rel_err n_unitigs,rel_err e_size" << endl;
    	outputFile[a] << std::fixed;
    	outputFile[a] << std::setprecision(2);
    } 

    uint64_t sum_n_nodes_h = 0;
    uint32_t values_n_nodes_h = 0;
    uint64_t n_nodes_h = 0;

    // some initial sampling in order to determina the "BEST" number of nodes of the graph
    if (max_abundance >= _last_a)
    {
    	if (verbose)
    	{
    		cout << "*** Sampling for a few values of k to determine the 'best' number of nodes of the graph" << endl;	
    	}
    	
	    for (uint32_t k = 21; k <= ((reads_total_content / reads_number) - 40); k = k + 10)
	    {
	    	cout << "***    k = " << k << endl;
	    	vector<double> n_nodes(max_abundance + 1,0); 
			vector<double> n_nodes_error(max_abundance + 1,LARGE_NUMBER);
			vector<double> n_unitigs(max_abundance + 1,0);
			vector<double> n_unitigs_error(max_abundance + 1,LARGE_NUMBER);
	 		vector<double> avg_unitig_length(max_abundance + 1,0);
	 		vector<double> avg_unitig_length_error(max_abundance + 1,LARGE_NUMBER);
	 		vector<double> e_size(max_abundance + 1,0);
	 		vector<double> e_size_error(max_abundance + 1,LARGE_NUMBER);

	 		// sampling
	 		sample_nodes(rlcsa, 
	 			k, 
	 			_first_a, 
	 			_last_a,
	 			reads, 
	 			reads_total_content, 
	 			reads_number, 
	 			reads_max_length, 
	 			n_nodes, 
	 			n_nodes_error,
	 			n_unitigs, 
	 			n_unitigs_error,
	 			avg_unitig_length,
	 			avg_unitig_length_error,
	 			e_size,
	 			e_size_error,
	 			relative_error,
	 			n_nodes_h,
	 			verbose);

	 		for (uint32_t a = _first_a; a <= _last_a; a++)
	 		{
	 			sum_n_nodes_h += n_nodes[a];
				values_n_nodes_h++;
				cout << k << " " << a << " n_nodes=" << n_nodes[a] << "±" << n_nodes_error[a] << "%" << endl;
	 		}
	
	    }
	    n_nodes_h = (double)sum_n_nodes_h / values_n_nodes_h;
	    if (verbose)
	    {
	    	cout << "Number of nodes used in pruning the search space over pairs (k,a): " << n_nodes_h << endl;
	    	cout << "If for a pair (k,a) the estimated number of nodes is less than " << _n_nodes_proportion << " * " << n_nodes_h << " = " << _n_nodes_proportion * n_nodes_h << ", then that pair is abandoned" << endl;
	    }
    }
    ///////////////////////////////////////////////////////////////    


 	for (uint32_t k = mink; k <= maxk; k++)
 	{
 		vector<double> n_nodes(max_abundance + 1,0); 
		vector<double> n_nodes_error(max_abundance + 1,LARGE_NUMBER);
		vector<double> n_unitigs(max_abundance + 1,0);
		vector<double> n_unitigs_error(max_abundance + 1,LARGE_NUMBER);
 		vector<double> avg_unitig_length(max_abundance + 1,0);
 		vector<double> avg_unitig_length_error(max_abundance + 1,LARGE_NUMBER);
 		vector<double> e_size(max_abundance + 1,0);
 		vector<double> e_size_error(max_abundance + 1,LARGE_NUMBER);

 		// sampling
 		sample_nodes(rlcsa, 
 			k, 
 			min_abundance, 
 			max_abundance, 
 			reads, 
 			reads_total_content, 
 			reads_number, 
 			reads_max_length, 
 			n_nodes, 
 			n_nodes_error,
 			n_unitigs, 
 			n_unitigs_error,
 			avg_unitig_length,
 			avg_unitig_length_error,
 			e_size,
 			e_size_error,
 			relative_error,
 			n_nodes_h,
 			verbose);

 		if (verbose)
 		{
 			cout << "*** " << currentDateTime() << endl;
 		}
 		// printing the results
 		for (uint32_t a = min_abundance; a <= max_abundance; a++)
 		{
 			if (n_nodes[a] < _n_nodes_proportion * n_nodes_h)
 			{
 				e_size_error[a] = LARGE_NUMBER;
 			}

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << (uint64_t)n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << "." << ","; // number of edges
	 		outputFile[a] << "." << ","; // average number of internal nodes in unitigs
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? double_to_string(avg_unitig_length[a]) : ".") << ","; // average length of unitigs // OLD: avg_unitig_length + k + 1 << ","; 
			outputFile[a] << "." << ","; // estimated sample size for kmers
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? int_to_string(n_unitigs[a]) : "." ) << ","; // number of unitigs: we don't print it if the number of nodes of the graph is too small
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? double_to_string(e_size[a]) : ".") << ","; // e-size	
			outputFile[a] << n_nodes_error[a] << ","; // rel_err for n_nodes
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? double_to_string(avg_unitig_length_error[a]) : ".") << ","; // rel_err avg_length_unitigs
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? double_to_string(n_unitigs_error[a]) : ".") << ","; // rel_err n_unitigs
			outputFile[a] << (e_size_error[a] != LARGE_NUMBER ? double_to_string(e_size_error[a]) : "."); // rel_err e_size				
			outputFile[a] << endl; 
			//outputFile[a].flush();


			if (verbose)
			{
				cout << k << " " << a << " ";
				cout << "n_nodes=" << (uint64_t)n_nodes[a] << "±" << (uint32_t)ceil(100 * n_nodes_error[a]) << "% ";
				if (e_size_error[a] != LARGE_NUMBER)
				{
					cout << "n_unitigs=" << (uint64_t)n_unitigs[a] << "±" << (uint32_t)ceil(100 * n_unitigs_error[a]) << "% ";
				}
				if (e_size_error[a] != LARGE_NUMBER)
				{
					cout << "avg_length=" << avg_unitig_length[a] << "±" << (uint32_t)ceil(100 * avg_unitig_length_error[a]) << "% ";	
				}
				if (e_size_error[a] != LARGE_NUMBER)
				{
					cout << "e_size=" << e_size[a] << "±" << (uint32_t)ceil(100 * e_size_error[a]) << "%";	
				}
				cout << endl;
				cout.flush();	
			}
 		}
 		if (not verbose)
 		{
 			cout << "*** Processed k = " << k << endl;
 		}
 	}
 	
 	for (uint32_t a = min_abundance; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}

 	delete rlcsa;

	return EXIT_SUCCESS;
}