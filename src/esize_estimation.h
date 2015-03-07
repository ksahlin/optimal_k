
#include "utils.h"

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
		extend_unitig_from_node(out_neighbors[abundance][0], rlcsa, abundance, out_length, out_abundance, 'o');
		extend_unitig_from_node(in_neighbors[abundance][0], rlcsa, abundance, in_length, in_abundance, 'i');

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
	extend_unitig_from_node(neighbor, rlcsa, abundance, temp_length, temp_abundance, direction);

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
	uint64_t u_length[]
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
				max_alive_abundance = a - 1;
				break;
			}
			// if unary
			if ((out_degree[a] == 1) and (in_degree[a] == 1))
			{
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
					min_alive_abundance = a + 1;
				}
			}
		}
	}
}

inline void get_unitig_stats_SMART(const string& node,
	const uint32_t& node_abundance,
	const RLCSA* rlcsa, 
	const uint32_t& min_abundance,
	const uint32_t& max_abundance,
	vector< vector<uint64_t> > &u_length
)
{
	if (min_abundance > max_abundance)
	{
		return;
	}
	// ALREADY INITIALIZED IN THE CALLING FUNCTION: u_length[a].clear(); for all al
	vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);
	vector<uint32_t> out_abundances;
	get_in_out_neighbors_with_abundances(node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors, out_abundances);

	uint32_t min_abundance_for_which_node_is_start = max_abundance + 1;
	uint32_t max_abundance_for_which_node_is_start = min_abundance - 1;

	for (uint32_t a = min_abundance; a <= max_abundance; a++)
	{
		// if is truly start of some unitig
		if ((out_neighbors[a].size() > 1) or ((out_neighbors[a].size() == 1) and (in_neighbors[a].size() != 1)))
		{
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
		if ((out_neighbors[a].size() == 0) and (in_neighbors[a].size() == 0))
		{
			u_length[a].push_back(1);
		}
	}

	// if there is some abundance for which 
	// node is the start node of some unitig
	if (min_abundance_for_which_node_is_start <= max_abundance_for_which_node_is_start)
	{
		uint64_t extended_length[max_abundance_for_which_node_is_start + 1];
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
			extend_unitig_from_node_SMART(neighbor, rlcsa, min_abundance_for_which_node_is_start, for_limit, extended_length);

			for (uint32_t a = min_abundance_for_which_node_is_start; a <= for_limit; a++)
			{
				u_length[a].push_back(1 + extended_length[a]);	
			}
			nbr_idx++;
		}
	}
}
