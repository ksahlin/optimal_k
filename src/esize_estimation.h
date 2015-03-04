
#include "utils.h"

inline void get_in_out_neighbors(const string& node, 
	const RLCSA* rlcsa, 
	const int &min_abundance,
	const int &max_abundance,
	vector< vector<string> > &in_neighbors,
	vector< vector<string> > &out_neighbors
	)
{
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		in_neighbors[a].clear();
		out_neighbors[a].clear();	
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
	 		in_neighbors[a].push_back(neighbor);
	 	}
	}

	for (auto nucl : {'A','C','G','T'})
	{
		neighbor = node.substr(1,node.length()-1) + nucl;
		int neighbor_abundance = calc_abundance(rlcsa,neighbor);
		neighbor_abundance = MIN(max_abundance, neighbor_abundance);
	 	for (int a = min_abundance; a <= neighbor_abundance; a++)
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


inline void extend_unitig_from_node_SMART(const string& node,
	const RLCSA* rlcsa, 
	const int min_abundance,
	const int max_abundance,
	vector<uint64_t> &u_length
)
{
	vector<bool> alive_abunance(max_abundance + 1, true);
	bool exists_alive_abundance = true;

	for (int a = min_abundance; a <= max_abundance; a++)
	{
		u_length[a] = 0;	
	}
	string current_node = node;
	vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);

	while (exists_alive_abundance)
	{
		get_in_out_neighbors(current_node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors);
		exists_alive_abundance = false;

		for (int a = min_abundance; a <= max_abundance; a++)
		{
			if (alive_abunance[a])
			{
				// if unary
				if ((in_neighbors[a].size() == 1) and (out_neighbors[a].size() == 1))
				{
					u_length[a] += 1;
					current_node = out_neighbors[a][0];
					exists_alive_abundance = true;
				} 
				else 
				{
					u_length[a] += 1;
					alive_abunance[a] = false;					
				}
			}
		}
	}
}

inline void get_unitig_stats_SMART(const string& node,
	const RLCSA* rlcsa, 
	const int min_abundance,
	const int max_abundance,
	vector< vector<uint64_t> > &u_length
)
{
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		u_length[a].clear();	
	}
	
	vector< vector<string> > in_neighbors(max_abundance + 1), out_neighbors(max_abundance + 1);
	get_in_out_neighbors(node, rlcsa, min_abundance, max_abundance, in_neighbors, out_neighbors);

	vector<int> abundances_for_which_start;
	for (int a = min_abundance; a <= max_abundance; a++)
	{
		// if is truly start of some unitig
		if ((out_neighbors[a].size() > 1) or ((out_neighbors[a].size() == 1) and (in_neighbors[a].size() != 1)))
		{
			abundances_for_which_start.push_back(a);
		}
		if ((out_neighbors[a].size() == 0) and (in_neighbors[a].size() == 0))
		{
			u_length[a].push_back(1);
		}
	}

	if (abundances_for_which_start.size() > 0)
	{
		// extend unitig to the each possible out-neighbor

		for (auto neighbor : out_neighbors[min_abundance])
		{
			vector<uint64_t> extended_length(max_abundance + 1, 0);
			extend_unitig_from_node_SMART(neighbor, rlcsa, min_abundance, max_abundance, extended_length);

			for (auto a : abundances_for_which_start)
			{
				u_length[a].push_back(1 + extended_length[a]);
			}
		}
	}
}
