
#include "utils.h"

struct assembled_string_t
{
    uint64_t length;
    string last_node;
};

inline void get_in_out_degrees_and_unique_out_neighbors(const string& node, 
	const RLCSA* rlcsa, 
	const uint32_t &min_abundance,
	const uint32_t &max_abundance,
	char in_degree[],
	char out_degree[],
	char out_neighbor_char[]
)
{
	uint32_t a;
	string neighbor = 'X' + node.substr(0,node.length()-1);
	uint32_t neighbor_abundance;

	for (a = min_abundance; a <= max_abundance; a++)
	{
		in_degree[a] = 0;
		out_degree[a] = 0;
	}

	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[0] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		in_degree[a]++;
	 	}
	}

	neighbor = node.substr(1,node.length()-1) + 'X';
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[neighbor.length()-1] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		out_degree[a]++;
	 		if (out_degree[a] == 1)
	 		{
	 			out_neighbor_char[a] = nucl;		
	 		}
	 	}
	}
}

inline void get_in_out_neighbors_with_abundances(const string& node, 
	const RLCSA* rlcsa, 
	const uint32_t &min_abundance,
	const uint32_t &max_abundance,
	vector< vector<string> > &in_neighbors,
	vector< vector<string> > &out_neighbors,
	vector<uint32_t> &out_abundances
	)
{
	string neighbor;
	uint32_t neighbor_abundance;
	uint32_t a;

	for (a = min_abundance; a <= max_abundance; a++)
	{
		in_neighbors[a].clear();
		in_neighbors[a].reserve(4);
		out_neighbors[a].clear();
		out_neighbors[a].reserve(4);
	}
	out_abundances.clear();
	out_abundances.reserve(4);

	neighbor = 'X' + node.substr(0,node.length()-1);
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[0] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		in_neighbors[a].push_back(neighbor);
	 	}
	}

	neighbor = node.substr(1,node.length()-1) + 'X';
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[neighbor.length()-1] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		out_neighbors[a].push_back(neighbor);
	 	}
	 	if (neighbor_abundance >= min_abundance)
	 	{
	 		out_abundances.push_back(neighbor_abundance);	
	 	}
	}
}

inline void get_in_out_neighbors(const string& node, 
	const RLCSA* rlcsa, 
	const uint32_t &min_abundance,
	const uint32_t &max_abundance,
	vector< vector<string> > &in_neighbors,
	vector< vector<string> > &out_neighbors
	)
{
	string neighbor = 'X' + node.substr(0,node.length()-1);
	uint32_t neighbor_abundance;
	uint32_t a;

	for (a = min_abundance; a <= max_abundance; a++)
	{
		in_neighbors[a].clear();
		in_neighbors[a].reserve(4);
		out_neighbors[a].clear();
		in_neighbors[a].reserve(4);
	}
	
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[0] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		in_neighbors[a].push_back(neighbor);
	 	}
	}

	neighbor = node.substr(1,node.length()-1) + 'X';
	for (auto nucl : {'A','C','G','T'})
	{
		neighbor[neighbor.length()-1] = nucl;
		neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (a = min_abundance; a <= neighbor_abundance; a++)
	 	{
	 		out_neighbors[a].push_back(neighbor);
	 	}
	}
}

inline void extend_unitig_from_node_OLD(const string& node,
	const RLCSA* rlcsa, 
	const uint32_t abundance,
	uint64_t &u_length,
	double &u_abundance,
	const char direction
)
{
	u_length = 0;
	u_abundance = 0;

	string current_node = node;

	vector< vector<string> > in_neighbors(abundance + 1), out_neighbors(abundance + 1);
	while (true)
	{
		get_in_out_neighbors(current_node, rlcsa, abundance, abundance, in_neighbors, out_neighbors);

		// if unary
		if ((in_neighbors[abundance].size() == 1) and (out_neighbors[abundance].size() == 1))
		{
			u_length += 1;
			u_abundance += calc_abundance(rlcsa, current_node);
			if (direction == 'o')
			{
				current_node = out_neighbors[abundance][0];
			} 
			else if (direction == 'i')
			{
				current_node = in_neighbors[abundance][0];
			}
		} 
		else
		{
			u_length += 1;
			u_abundance += calc_abundance(rlcsa, current_node) / (double)(out_neighbors[abundance].size() + in_neighbors[abundance].size());
			break;
		}
	}
}

inline void get_unitig_stats_OLD(const string& node,
	const uint32_t& node_abundance,
	const RLCSA* rlcsa, 
	const uint32_t abundance,
	uint64_t &u_length,
	double &u_abundance,
	default_random_engine& generator
)
{
	u_length = 0;
	u_abundance = 0;
	vector< vector<string> > in_neighbors(abundance + 1), out_neighbors(abundance + 1);
	get_in_out_neighbors(node, rlcsa, abundance, abundance, in_neighbors, out_neighbors);

	// check if node is isolated
	if ((in_neighbors[abundance].size() == 0) and (out_neighbors[abundance].size() == 0))
	{
		u_length = 1;
		u_abundance = node_abundance; //calc_abundance(rlcsa, node);
		return;
	}
	
	// check if node is unary
	if ((in_neighbors[abundance].size() == 1) and (out_neighbors[abundance].size() == 1))
	{
		// check if self-loop
		if (node == out_neighbors[abundance][0])
		{
			u_length = 2;
			u_abundance = node_abundance; //calc_abundance(rlcsa, node);
			return;	
		}

		uint64_t out_length, in_length;
		double out_abundance, in_abundance;
		extend_unitig_from_node_OLD(out_neighbors[abundance][0], rlcsa, abundance, out_length, out_abundance, 'o');
		extend_unitig_from_node_OLD(in_neighbors[abundance][0], rlcsa, abundance, in_length, in_abundance, 'i');

		u_length = 1 + out_length + in_length;
		u_abundance = node_abundance + out_abundance + in_abundance; // calc_abundance(rlcsa, node)
		return;
	}
	
	// node is not unary and not isolated
	
	char direction;
	string neighbor;
	// choose a random out-neighbor
	// uint32_t random_index = (rand() / (double)RAND_MAX) * (in_neighbors[abundance].size() + out_neighbors[abundance].size());
	std::uniform_int_distribution<char> uniform_neighbor_distribution(0,in_neighbors[abundance].size() + out_neighbors[abundance].size() - 1);
	uint32_t random_index = uniform_neighbor_distribution(generator);

	if (random_index < in_neighbors[abundance].size()) // chosen an in-neighbor
	{
		direction = 'i';
		neighbor = in_neighbors[abundance][random_index];
	}
	else // chosen an out-neighbor
	{
		direction = 'o';
		neighbor = out_neighbors[abundance][random_index - in_neighbors[abundance].size()];
	}

	// compute stats for the path starting with node, neighbor
	uint64_t temp_length;
	double temp_abundance;
	extend_unitig_from_node_OLD(neighbor, rlcsa, abundance, temp_length, temp_abundance, direction);

	u_length = 1 + temp_length;
	u_abundance = temp_abundance + node_abundance / (double)(out_neighbors[abundance].size() + in_neighbors[abundance].size());

}

/***************************************************/
// This Function is simpler, but is is very slightly slower 
/***************************************************/
inline void extend_unitig_from_node_SMART_simpler_but_slightly_slower(const string& node,
	const RLCSA* rlcsa, 
	const uint32_t min_abundance,
	const uint32_t max_abundance,
	uint64_t u_length[]
)
{
	vector<bool> alive_abunance(max_abundance + 1, true);
	bool exists_alive_abundance = true;

	// ALREADY INITIALIZED IN THE CALLING FUNCTION
	// u_length[a] = 0 for all a
	string current_node = node;
	char out_neighbor_char[max_abundance + 1];
	char in_degree[max_abundance + 1], out_degree[max_abundance + 1];
	// vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);

	// assert(calc_abundance(rlcsa,node) >= max_abundance);

	while (exists_alive_abundance)
	{
		// get_in_out_neighbors(current_node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors);
		get_in_out_degrees_and_unique_out_neighbors(current_node, rlcsa, min_abundance, max_abundance, in_degree, out_degree, out_neighbor_char);

		exists_alive_abundance = false;
		for (uint32_t a = min_abundance; a <= max_abundance; a++)
		{
			// assert(in_degree[a] == in_neighbors[a].size());
			// assert(out_degree[a] == out_neighbors[a].size());
			if (alive_abunance[a])
			{
				u_length[a] += 1;
				// if not unary
				// if ((in_neighbors[a].size() != 1) and (out_neighbors[a].size() != 1))
				if ((in_degree[a] != 1) or (out_degree[a] != 1))
				{
					alive_abunance[a] = false;					
				}
				else if (not exists_alive_abundance)
				{
					current_node = current_node.substr(1,current_node.length()-1) + out_neighbor_char[a];
					exists_alive_abundance = true;	
					// assert(current_node == out_neighbors[a][0]);
				}				
			}
		}
	}
}

inline void extend_unitig_from_node_SMART(const string& node,
	const RLCSA* rlcsa, 
	const uint32_t min_abundance,
	const uint32_t max_abundance,
	uint64_t u_length[],
	string u_last_node[]
)
{
	uint32_t min_alive_abundance = min_abundance;
	uint32_t max_alive_abundance = max_abundance;
	// ALREADY INITIALIZED IN THE CALLING FUNCTION
	// u_length[a] = 0; for all a
	string current_node = node;
	char out_neighbor_char[max_abundance + 1];
	char in_degree[max_abundance + 1], out_degree[max_abundance + 1];
	bool updated_current_node;
	uint32_t a;

	while (min_alive_abundance <= max_alive_abundance)
	{
		get_in_out_degrees_and_unique_out_neighbors(current_node, rlcsa, min_alive_abundance, max_alive_abundance, in_degree, out_degree, out_neighbor_char);
		updated_current_node = false;
		
		for (a = min_alive_abundance; a <= max_alive_abundance; a++)
		{
			u_length[a]++;
			if (out_degree[a] == 0)
			{
				for (uint16_t a2 = a; a2 <= max_alive_abundance; a2++)
				{
					u_last_node[a2] = current_node;	
				} 
				max_alive_abundance = a - 1;
				break;
			}
			// if unary
			if ((out_degree[a] == 1) and (in_degree[a] == 1))
			{
				if ((u_length[a] >= 3000) and (current_node == node))
				{
					return;
				}
				if (not updated_current_node)
				{
					current_node = current_node.substr(1,current_node.length()-1) + out_neighbor_char[a];
					updated_current_node = true;
				}
			} 
			else 
			{
				if ((out_degree[a] > 1) or (in_degree[a] > 1))
				{
					u_last_node[a] = current_node;
					min_alive_abundance = a + 1;
				}
			}
		}
	}
}

inline void get_unitig_stats_SMART(const string& node,
	const uint32_t& node_abundance,
	const RLCSA* rlcsa, 
	const unordered_set<uint32_t> &alive_abundances,
	const uint32_t &max_abundance_where_node_is_present,
	vector< vector<uint64_t> > &u_length,
	bool flag_unitigs
)
{
	if (alive_abundances.size() == 0)
	{
		return;
	}
	uint32_t min_alive_abundance = *std::min_element(alive_abundances.begin(), alive_abundances.end());
	if (min_alive_abundance > max_abundance_where_node_is_present)
	{
		return;
	}
	uint32_t max_alive_abundance = *std::max_element(alive_abundances.begin(), alive_abundances.end());	

	// ALREADY INITIALIZED IN THE CALLING FUNCTION: u_length[a].clear(); for all al
	vector< vector<string> > in_neighbors(max_alive_abundance + 1), out_neighbors(max_alive_abundance + 1);
	vector<uint32_t> out_abundances, abundances_for_which_node_is_start;
	get_in_out_neighbors_with_abundances(node, rlcsa, min_alive_abundance, max_alive_abundance, in_neighbors, out_neighbors, out_abundances);

	uint32_t min_abundance_for_which_node_is_start = max_alive_abundance + 1;
	uint32_t max_abundance_for_which_node_is_start = min_alive_abundance - 1;

	for (auto a : alive_abundances)
	{
		if (a > max_abundance_where_node_is_present)
		{
			break;
		}
		// if is truly start of some unitig
		if ((out_neighbors[a].size() > 1) or ((out_neighbors[a].size() == 1) and (in_neighbors[a].size() != 1)))
		{
			abundances_for_which_node_is_start.push_back(a);
			if (a < min_abundance_for_which_node_is_start)
			{
				min_abundance_for_which_node_is_start = a;
			} 
			if (a > max_abundance_for_which_node_is_start)
			{
				max_abundance_for_which_node_is_start = a;
			}	
		} 
		else // if isolated node
		{
			// adding an isolated node only if we are in 'unitigs' mode
			if (flag_unitigs)
			{
				if ((out_neighbors[a].size() == 0) and (in_neighbors[a].size() == 0))
				{
					u_length[a].push_back(1);
				}	
			}	
		}
	}

	// if there is some abundance for which 
	// node is the start node of some unitig

	vector< vector<assembled_string_t> > assembled_strings(max_abundance_for_which_node_is_start + 1);

	if (min_abundance_for_which_node_is_start <= max_abundance_for_which_node_is_start)
	{
		uint64_t extended_length[max_abundance_for_which_node_is_start + 1];
		string u_last_node[max_abundance_for_which_node_is_start + 1];
		uint32_t nbr_idx = 0;

		// extend unitig towards each possible out-neighbor
		for (auto neighbor : out_neighbors[min_abundance_for_which_node_is_start])
		{
			uint32_t neighbor_abundance = out_abundances[nbr_idx]; // calc_abundance(rlcsa,neighbor); 
			uint32_t for_limit = MIN(max_abundance_for_which_node_is_start,neighbor_abundance);
			
			for (uint32_t a = min_abundance_for_which_node_is_start; a <= for_limit; a++)
			{
				extended_length[a] = 0;
			}
			extend_unitig_from_node_SMART(neighbor, rlcsa, min_abundance_for_which_node_is_start, for_limit, extended_length, u_last_node);

			for (auto a : abundances_for_which_node_is_start)
			{
				if (a > for_limit)
				{
					break;
				}

				if (flag_unitigs or (1 + extended_length[a] >= node.length()) )
				{
					assembled_string_t assembled_string;
					assembled_string.length = 1 + extended_length[a];
					assembled_string.last_node = u_last_node[a];
					assembled_strings[a].push_back(assembled_string);
				}
			}
			nbr_idx++;
		}
	}

	for (uint32_t a = min_abundance_for_which_node_is_start; a <= max_abundance_for_which_node_is_start; a++)
	{
		// checking whether we have a bubble in the CONTIGS mode
		if ((not flag_unitigs) and
			(out_neighbors[a].size() == 2) and (assembled_strings[a].size() == 2)    and
			(assembled_strings[a][0].length    == assembled_strings[a][1].length)    and
			(assembled_strings[a][0].last_node == assembled_strings[a][1].last_node)
			)
		{
			// getting the in/out-degrees of the last node
			vector< vector<string> > in_neighbors_last_node(a + 1), out_neighbors_last_node(a + 1);
			get_in_out_neighbors(assembled_strings[a][0].last_node, rlcsa, a, a, in_neighbors_last_node, out_neighbors_last_node);
			if (in_neighbors_last_node[a].size() == 2)
			{
				// we have found a bubble
				// cout << "-------- found a bubble of length " << assembled_strings[a][0].length << endl;

				// checking whether the bubble is extendable
				if ((in_neighbors[a].size() == 1) and (out_neighbors_last_node[a].size() == 1))
				{
					// the bubble is extendable
					uint64_t extended_length_right[a + 1];
					uint64_t extended_length_left[a + 1];
					extended_length_right[a] = 0;
					extended_length_left[a] = 0;
					string u_last_node[a + 1];
					
					extend_unitig_from_node_SMART(assembled_strings[a][0].last_node, rlcsa, a, a, extended_length_right, u_last_node);
					extend_unitig_from_node_SMART(reverse_complement(node), rlcsa, a, a, extended_length_left, u_last_node);
					
					uint64_t final_length = assembled_strings[a][0].length + extended_length_right[a] + extended_length_left[a];
					// cout << "-------- the bubble is extendable to length " << final_length << endl;
					u_length[a].push_back(final_length);
				}
				else
				{
					// the bubble is not extendable
					// we can report only the first path (once)
					u_length[a].push_back(assembled_strings[a][0].length);
				}
			}


		}
		else
		{
			// no bubble, we can report all assembled strings 
			for (auto assembled_string : assembled_strings[a])
			{
				u_length[a].push_back(assembled_string.length);	
			}	
		}
		
	}
}
