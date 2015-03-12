
#include "utils.h"
#include "esize_estimation.h"

#define TWOSIDED95P_QUANTILE 1.96
#define MIN_SAMPLE_SIZE_UNITIGS 500
#define MAX_SAMPLE_SIZE_UNITIGS 5000
#define RATIO_SAMPLE_SIZE_UNITIGS 0.05

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
	const vector<string>& reads,
	const vector<uint64_t>& sample_size_kmers,
	const vector<uint64_t>& sample_size_unitigs,
	vector<double> &n_uint32_ternal,
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
	std::uniform_int_distribution<int> uniform_pos_distribution(0,reads[0].length() - k);

	vector<double> n_sampled_nodes_weighted(max_abundance + 1, 0);

	vector<double> n_uint32_ternal_local(max_abundance + 1, 0);
	vector<double> n_starts_local(max_abundance + 1, 0);
	vector<double> e_size_sum_length(max_abundance + 1, 0);
	vector<double> e_size_sum_length_squared(max_abundance + 1, 0);

	uint64_t n_sampled_kmers = 0;
	vector<uint64_t> n_sampled_nodes(max_abundance + 1, 0);
	vector<uint64_t> n_sampled_unitigs(max_abundance + 1, 0);
	bool sampled_enough_kmers = false;
	bool sampled_enough_unitigs = false;
	
	uint32_t a_w_max_ss_kmers = -1;
	uint64_t max_sample_size_kmers = 0;
	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		if (sample_size_kmers[a] >= max_sample_size_kmers)
		{
			max_sample_size_kmers = sample_size_kmers[a];
			a_w_max_ss_kmers = a;
		}
	}

	uint32_t min_abundance_unitigs = min_abundance;
	uint64_t n_total_sampled_unitigs = 0;
	// omp_set_dynamic(0);
	// shared(sampled_enough_kmers,sampled_enough_unitigs)
	#pragma omp parallel for num_threads(N_THREADS)
	for (uint64_t i = 0; i < 100 * n_reads; i++)
	{
		if ((not sampled_enough_kmers) or (not sampled_enough_unitigs))
		{
			if (n_sampled_kmers >= sample_size_kmers[a_w_max_ss_kmers])
			{
				#pragma omp critical
				sampled_enough_kmers = true;
			}	
			
			uint64_t read_index = uniform_read_distribution(generator);
			string read = reads[read_index];
			// // MAKE SURE THIS IS OK!
			// if (rand() / (double)RAND_MAX < 0.5)
			// {
			// 	read = reverse_complement(read);
			// }
			uint32_t pos = uniform_pos_distribution(generator);
	        string sample = read.substr(pos,k);

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
    			get_unitig_stats_SMART(sample, sample_abundance, rlcsa, min_abundance_unitigs, for_limit, u_length);
				
    			for (uint32_t a = min_abundance_unitigs; a <= for_limit; a++)
    			{
    				#pragma omp atomic
    				n_sampled_unitigs[a] += u_length[a].size();
    				#pragma omp atomic
    				n_total_sampled_unitigs += u_length[a].size();
    				for (auto length : u_length[a])
    				{
    					#pragma omp atomic
   						e_size_sum_length[a] += (double)(length + k - 1) / sample_abundance;;
   						#pragma omp atomic
    					e_size_sum_length_squared[a] += (length + k - 1) * (double)(length + k - 1) / sample_abundance;;
    				}
    			}
    			// only Master gets to update which is the min_abundance_unitigs
    			if (omp_get_thread_num() == 0)
    			{
    				if (n_sampled_unitigs[min_abundance_unitigs] >= sample_size_unitigs[min_abundance_unitigs])
	    			{
	    				min_abundance_unitigs++;
	    				// cout << "min_abundance_unitigs = " << min_abundance_unitigs << endl;
	    				// cout << "n_total_sampled_unitigs = " << n_total_sampled_unitigs << endl;
	    				// for (uint32_t a = min_abundance; a <= max_abundance; a++)
	    				// {
	    				// 	cout << "n_sampled_unitigs[" << a << "]=" << n_sampled_unitigs[a] << endl;
	    				// }
	    			}	
    			}
    			if (min_abundance_unitigs > max_abundance)
	    		{
	    			#pragma omp critical
	    			sampled_enough_unitigs = true;		
	    		}
    		}

    		// computing the other estimates
        	vector<uint32_t> in_degree(max_abundance + 1), out_degree(max_abundance + 1);
        	get_in_out_degrees(sample,rlcsa,min_abundance,max_abundance,in_degree,out_degree);

        	for (uint32_t a = min_abundance; a <= for_limit; a++)
        	{
        		#pragma omp atomic
       			n_sampled_nodes_weighted[a] += sample_weight;

        		// is uint32_ternal
    			if ((out_degree[a] == 1) and (in_degree[a] == 1)) 
            	{
            		#pragma omp atomic
           			n_uint32_ternal_local[a] += sample_weight;
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

	cout << "Sampled " << n_total_sampled_unitigs << " unitigs " << endl;

	if (not sampled_enough_kmers)
	{
		cout << "*** I could sample only " << n_sampled_kmers << " kmers out of " <<  sample_size_kmers[a_w_max_ss_kmers] << endl;
		assert(false);
	}

	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		n_uint32_ternal[a] = n_uint32_ternal_local[a];
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
    uint32_t mink, maxk;
    string readFileName, outputFileName;
    string buildindex, loadindex;
	uint32_t min_abundance,max_abundance;
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

	parser.add_option("-r", "--readfile") .type("string") .dest("r") .set_default("") .help("a file containing a list of FASTA/Q(.gz) file names, one per line (all reads need to have the same length)");
	parser.add_option("-o", "--outputfile") .type("string") .dest("o") .set_default("") .help("output file");
	parser.add_option("-a", "--minabundance") .type("uint32_t") .dest("a") .action("store") .set_default(1) .help("try all abundances starting with this value (default: %default)");
	parser.add_option("-A", "--maxabundance") .type("uint32_t") .dest("A") .action("store") .set_default(5) .help("try all abundances up to this value (default: %default)");
	parser.add_option("-t", "--threads") .type("uint32_t") .dest("t") .action("store") .set_default(0) .help("number of threads; use 0 for all cores (default: %default)");
	parser.add_option("-k", "--mink") .type("uint32_t") .dest("k") .action("store") .set_default(15) .help("try all kmer sizes starting with this value (default: %default)");
	parser.add_option("-K", "--maxk") .type("uint32_t") .dest("K") .action("store") .set_default(0) .help("try all kmer sizes up to this value (default: read_length - 10)");
	parser.add_option("-e", "--relerror") .type("float") .dest("e") .action("store") .set_default(0.1) .help("relative error of the estimations (default: %default)");
	parser.add_option("-b", "--buildindex"). type("string") .dest("buildindex") .action("store") .set_default("") .help("the filename where the index should be saved");
	parser.add_option("-l", "--loadindex"). type("string") .dest("loadindex") .action("store") .set_default("") .help("the filename from where the index should be loaded");
	optparse::Values& options = parser.parse_args(argc, argv);

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
		cout << "*** ERROR: -r <read_file> must be given" << endl;	
		return EXIT_FAILURE;
	}
	// checking if we have the index
	if ((buildindex == "") and (loadindex == ""))
	{
		cout << "*** ERROR: either --buildindex <file_name> or --loadindex <file_name> must be given" << endl;
		return EXIT_FAILURE;
	}
	// checking if we have the index
	if ((buildindex != "") and (loadindex != ""))
	{
		cout << "*** ERROR: not both of --buildindex <file_name> or --loadindex <file_name> must be given" << endl;
		return EXIT_FAILURE;
	}

	// compact_read cread = encode_string("AAA");
	// return EXIT_SUCCESS;

	// if we need to build the index
	if (buildindex != "")
	{
		get_data_and_build_rlcsa_iterative(readFileName, buildindex, N_THREADS);
		// get_data_and_build_rlcsa_noniterative(readFileName, buildindex, N_THREADS);
		cout << "*** SUCCESS: now run the program with --loadindex " << buildindex << endl;
		return EXIT_SUCCESS;		
	}

	if (outputFileName == "")
	{
		cerr << "The argument -o|--outputfile is needed" << endl;
		return EXIT_FAILURE;
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

 	vector<string> reads;
 	// we load the reads
 	if (EXIT_FAILURE == get_reads_using_Bank(readFileName, reads))
	{
		return EXIT_FAILURE;
	}

 	// these vectors get re-written for each value of k
 	vector<uint64_t> sample_size_kmers(max_abundance + 1, 0);
 	vector<uint64_t> sample_size_unitigs(max_abundance + 1, MAX_SAMPLE_SIZE_UNITIGS);
 	vector<double> n_uint32_ternal(max_abundance + 1,1), n_starts(max_abundance + 1,1);
 	vector<double> n_nodes(max_abundance + 1,1), n_unitigs(max_abundance + 1,1);
 	vector<double> e_size(max_abundance + 1,1);

    if (maxk == 0)
    {
    	maxk = reads[0].length() - 10;
    }

 	cout << "*** Writing results to files:" << endl;
    vector<ofstream> outputFile(max_abundance + 1);
    for (uint32_t a = min_abundance; a <= max_abundance; a++)
    {
    	cout << "***    " << outputFileName + ".mink" + int_to_string(mink) + ".maxk" + int_to_string(maxk) + ".a" + int_to_string(a) + ".metrics.csv" << endl;
    	outputFile[a].open((outputFileName + ".mink" + int_to_string(mink) + ".maxk" + int_to_string(maxk) + ".a" + int_to_string(a) + ".metrics.csv").c_str());
    	outputFile[a] << "k,a,nr_nodes,nr_edges,avg_uint32_ternal_nodes,avg_length_unitigs,est_sample_size,nr_unitigs,e_size" << endl;
    } 

 	for (uint32_t k = mink; k <= maxk; k++)
 	{
 		// getting the sample size
 		for (uint32_t a = min_abundance; a <= max_abundance; a++)
 		{
 			uint64_t SS_n_nodes = get_sample_size_for_proportion(n_nodes[a] / (double)(reads.size() * (reads[0].length() - k + 1)), relative_error);
 			uint64_t SS_n_unitigs = get_sample_size_for_proportion(n_unitigs[a] / (double)(reads.size() * (reads[0].length() - k + 1)), relative_error);
 			uint64_t SS_avg_nodes_unitig = get_sample_size_for_proportion(n_starts[a] / (double)(n_uint32_ternal[a] + n_starts[a]), 1 / (double)(1 + relative_error) - 1);
 			
 			// cout << "SS_n_nodes = " << SS_n_nodes << endl;
 			// cout << "SS_n_unitigs = " << SS_n_unitigs << endl;
 			// cout << "ratio = " << n_starts[a] / (double)(n_uint32_ternal[a] + n_starts[a]) << endl;
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
 		sample_nodes(rlcsa, k, min_abundance, max_abundance, reads, sample_size_kmers, sample_size_unitigs, n_uint32_ternal, n_starts, n_nodes, n_unitigs, e_size);	

 		// pruint32_ting the results
 		for (uint32_t a = min_abundance; a <= max_abundance; a++)
 		{
 			assert(n_starts[a] > 0);
			double avg_nodes_unitig = n_uint32_ternal[a] / (double)n_starts[a];

	 		outputFile[a] << k << ",";
	 		outputFile[a] << a << ",";
	 		outputFile[a] << (uint64_t)n_nodes[a] << ","; // number of nodes
	 		outputFile[a] << ".,"; // number of edges
	 		outputFile[a] << avg_nodes_unitig << ","; // average number of uint32_ternal nodes in unitigs
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
 	
 	for (uint32_t a = min_abundance; a <= max_abundance; a++)
 	{
 		outputFile[a].close();	
 	}

	return EXIT_SUCCESS;
}