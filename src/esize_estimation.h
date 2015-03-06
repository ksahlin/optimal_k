
#include "utils.h"

inline void get_in_out_degrees_and_unique_out_neighbors(const string& node, 
	const RLCSA* rlcsa, 
	const int &min_abundance,
	const int &max_abundance,
	int in_degree[],
	int out_degree[],
	char out_neighbor_char[]
)
{
	int a;
	string neighbor = 'X' + node.substr(0,node.length()-1);
	int neighbor_abundance;

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


inline void get_in_out_neighbors(const string& node, 
	const RLCSA* rlcsa, 
	const int &min_abundance,
	const int &max_abundance,
	vector< vector<string> > &in_neighbors,
	vector< vector<string> > &out_neighbors
	)
{
	string neighbor = 'X' + node.substr(0,node.length()-1);
	int neighbor_abundance;
	int a;

	for (a = min_abundance; a <= max_abundance; a++)
	{
		in_neighbors[a].clear();
		out_neighbors[a].clear();	
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

inline void extend_unitig_from_node(const string& node,
	const RLCSA* rlcsa, 
	const int abundance,
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

inline void get_unitig_stats(const string& node,
	const int& node_abundance,
	const RLCSA* rlcsa, 
	const int abundance,
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
	// uint random_index = (rand() / (double)RAND_MAX) * (in_neighbors[abundance].size() + out_neighbors[abundance].size());
	std::uniform_int_distribution<char> uniform_neighbor_distribution(0,in_neighbors[abundance].size() + out_neighbors[abundance].size() - 1);
	uint random_index = uniform_neighbor_distribution(generator);

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

inline void extend_unitig_from_node_SMART_very_slightly_slower(const string& node,
	const RLCSA* rlcsa, 
	const int min_abundance,
	const int max_abundance,
	uint64_t u_length[]
)
{
	vector<bool> alive_abunance(max_abundance + 1, true);
	bool exists_alive_abundance = true;

	// ALREADY INITIALIZED IN THE CALLING FUNCTION
	// u_length[a] = 0 for all a
	string current_node = node;
	char out_neighbor_char[max_abundance + 1];
	int in_degree[max_abundance + 1], out_degree[max_abundance + 1];
	// vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);

	assert(calc_abundance(rlcsa,node) >= max_abundance);

	while (exists_alive_abundance)
	{
		// get_in_out_neighbors(current_node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors);
		get_in_out_degrees_and_unique_out_neighbors(current_node, rlcsa, min_abundance, max_abundance, in_degree, out_degree, out_neighbor_char);

		exists_alive_abundance = false;
		for (int a = min_abundance; a <= max_abundance; a++)
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

// THIS FUNCTION IS NOT USED
// this function is extending the unitig starting at node
// by keeping an interval of alive abundances
// but in practice it is only slightly faster
// *******************************************************************
inline void extend_unitig_from_node_SMART(const string& node,
	const RLCSA* rlcsa, 
	const int min_abundance,
	const int max_abundance,
	uint64_t u_length[]
)
{
	int min_alive_abundance = min_abundance;
	int max_alive_abundance = max_abundance;
	// ALREADY INITIALIZED IN THE CALLING FUNCTION
	// u_length[a] = 0; for all a
	string current_node = node;
	char out_neighbor_char[max_abundance + 1];
	int in_degree[max_abundance + 1], out_degree[max_abundance + 1];
	bool updated_current_node;

	while (min_alive_abundance <= max_alive_abundance)
	{
		get_in_out_degrees_and_unique_out_neighbors(current_node, rlcsa, min_alive_abundance, max_alive_abundance, in_degree, out_degree, out_neighbor_char);
		//assert(in_neighbors[min_alive_abundance].size() > 0);
		//assert(in_degree[min_alive_abundance] > 0);
		updated_current_node = false;

		for (int a = min_alive_abundance; a <= max_alive_abundance; a++)
		{
			// cout << "aaa" << endl;
			u_length[a] += 1;
			// if unary
			if ((out_degree[a] == 1) and (in_degree[a] == 1))
			{
				if (not updated_current_node)
				{
					current_node = current_node.substr(1,current_node.length()-1) + out_neighbor_char[a];
					updated_current_node = true;
					// cout << current_node << endl;	
				}
			} 
			else 
			{
				if ((out_degree[a] > 1) or (in_degree[a] > 1))
				{
					min_alive_abundance = a + 1;
				}
				if (out_degree[a] == 0)
				{
					max_alive_abundance = a - 1;
					break;
				}
			}
		}
	}
}

inline void get_unitig_stats_SMART(const string& node,
	const int& node_abundance,
	const RLCSA* rlcsa, 
	const int& min_abundance,
	const int& max_abundance,
	vector< vector<uint64_t> > &u_length
)
{
	// ALREADY INITIALIZED IN THE CALLING FUNCTION
	// u_length[a].clear(); for all al
	
	vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);
	get_in_out_neighbors(node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors);

	int min_abundance_for_which_node_is_start = max_abundance + 1;
	vector<int> abundances_for_which_node_is_start;
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		// if is truly start of some unitig
		if ((out_neighbors[a].size() > 1) or ((out_neighbors[a].size() == 1) and (in_neighbors[a].size() != 1)))
		{
			abundances_for_which_node_is_start.push_back(a);
			if (a < min_abundance_for_which_node_is_start)
			{
				min_abundance_for_which_node_is_start = a;
			}
		} 
		else // if isolated node
		if ((out_neighbors[a].size() == 0) and (in_neighbors[a].size() == 0))
		{
			u_length[a].push_back(1);
		}
	}

	if (abundances_for_which_node_is_start.size() > 0)
	{
		uint64_t extended_length[max_abundance + 1];
		int neighbor_abundance;
		// extend unitig towards each possible out-neighbor
		for (auto neighbor : out_neighbors[min_abundance])
		{
			for (int a = min_abundance; a <= max_abundance; a++)
			{
				extended_length[a] = 0;
			}
			neighbor_abundance = calc_abundance(rlcsa,neighbor);
			if (neighbor_abundance >= min_abundance_for_which_node_is_start)
			{
				extend_unitig_from_node_SMART(neighbor, rlcsa, min_abundance, MIN(max_abundance,neighbor_abundance), extended_length);

				for (auto a : abundances_for_which_node_is_start)
				{
					if (neighbor_abundance >= a)
					{
						u_length[a].push_back(1 + extended_length[a]);	
					}
				}
			}
		}
	}
}
