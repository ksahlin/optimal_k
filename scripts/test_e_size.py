
import networkx as nx
import random
import sys
SAMPLE_SIZE=1000

def Generate_samples(G):
	sample_vector = []

	#add weighted sample vector
	for n in G.nodes():
		# print n
		# print 'lol',G[n]
		# print G.node[n]['a']
		for i in range(G.node[n]['a']):
			sample_vector.append(n)

	samples = []
	for sample in range(SAMPLE_SIZE):
		samples.append(random.choice(sample_vector))
	print 'Samples:', samples
	return samples

def get_path_stats(G,G_rev_comp,n,start,seen_start= 0):
	# use this to check for simple loop
	length = 0
	tot_abundance = 0
	while True:
		#print n
		if n == start:
			length += 1
			seen_start[0] += 1
			nbr_count = len(G.neighbors(n)) + len(G_rev_comp.neighbors(n)) 
			tot_abundance += G.node[n]['a'] / float(nbr_count)
			return length ,tot_abundance
		# print G.neighbors(n), G_rev_comp.neighbors(n)
		if len(G.neighbors(n)) == 1 and len(G_rev_comp.neighbors(n)) == 1:
			length += 1
			tot_abundance += G.node[n]['a']
			n = G.neighbors(n)[0]
		else:
			length += 1
			nbr_count = len(G.neighbors(n)) + len(G_rev_comp.neighbors(n)) 
			tot_abundance += G.node[n]['a'] / float(nbr_count)
			return length, tot_abundance

def compute_path_length(G,G_compl, start_node):
	#print start_node

	# start_node is isolated	
	if not G.neighbors(start_node) and not G_compl.neighbors(start_node):
		#print 'no nbrs:',start_node
		return 1, G.node[start_node]['a'] 

	# start_node is unary (internal)
	elif  len(G.neighbors(start_node)) == 1 and len(G_compl.neighbors(start_node)) == 1:
		seen_start = [0]
		#print 'unary:', start_node
		right_nbr = G.neighbors(start_node)[0]
		left_nbr = G_compl.neighbors(start_node)[0]
		# right path length
		right_length, right_abundance =  get_path_stats(G,G_compl,right_nbr,start_node,seen_start)
		# left path length
		left_length, left_abundance =  get_path_stats(G_compl,G,left_nbr, start_node, seen_start)
		if seen_start[0] < 2:
			a =  left_abundance + right_abundance + G.node[start_node]['a'] # the start node plus the two paths
			l = 1 + right_length + left_length # the node itself plus the length of its paths in each direction
			return l, a
		elif seen_start[0] == 2:
			a = G.node[start_node]['a'] + right_abundance  - G.node[start_node]['a']/2
			l = right_length 
			return l, a


	# start_node has at least two in or out neighbors (is etxtemity) 
	else:
		# chose extension randomly with equal weight
		nbr = random.choice( G.neighbors(start_node) + G_compl.neighbors(start_node))
		nbr_count = len(G.neighbors(start_node) + G_compl.neighbors(start_node))
		# the nbr we chose in in G
		if G.has_edge(start_node, nbr):
			l, a = get_path_stats(G, G_compl, nbr, start_node)

		# the nbr we chose in in G_compl
		elif G_compl.has_edge(start_node, nbr):
			l, a = get_path_stats(G_compl, G, nbr, start_node)
		else:
			print 'There is a bug'
			sys.exit()
		a =  a + G.node[start_node]['a']/float(nbr_count) # the start node plus the two paths
		l = 1 + l # the node itself plus the length of its paths in each direction

		return l,a
		
		

def main():
	G_initial, G_compl_initial = generate_initial_example()
	G_loop, G_compl_loop = generate_simple_loop()
	G_ab, G_compl_ab = generate_initial_diff_ab()

	for G,G_compl in [(G_initial, G_compl_initial), (G_loop, G_compl_loop), (G_ab, G_compl_ab)]:
		samples = Generate_samples(G)
		unitig_stats = []
		for sample in samples:
			path_length, abundance = compute_path_length(G,G_compl, sample)
			unitig_stats.append( (path_length, abundance) )

		lengths = map(lambda x: x[0], unitig_stats)
		abundances = map(lambda x: x[1], unitig_stats)
		#print lengths, abundances
		sum_lengths_squared = 0
		sum_lengths = 0
		for l,a in zip(lengths,abundances):
			sum_lengths += l/float(a)
			sum_lengths_squared += l**2/ float(a)
		e_size = sum_lengths_squared/sum_lengths
		print e_size

		## WHAT WE DERIVED
		lengths = map(lambda x: x[0], unitig_stats)
		abundances = map(lambda x: x[1], unitig_stats)
		#print lengths, abundances
		inner_sum = 0
		for l,a in zip(lengths,abundances):
			inner_sum += l*(l/float(a))

		normalization_const = sum(map(lambda (x,y): x/float(y) , zip(lengths, abundances)))
		e_size = inner_sum/normalization_const
		print 'New:',e_size



###EXAMPLES####

def generate_initial_example():

	G=nx.DiGraph()
	G.add_edges_from([(0,1),(1,2),(2,3),(2,4)])
	G.add_node(5)
	G.node[0]['a'] = 1
	G.node[1]['a'] = 2
	G.node[2]['a'] = 2
	G.node[3]['a'] = 1
	G.node[4]['a'] = 1
	G.node[5]['a'] = 3

	G_compl = nx.DiGraph()
	for n1,n2 in G.edges():
		G_compl.add_edge(n2,n1)
		G_compl.node[n1]['a'] = G.node[n1]['a']
		G_compl.node[n2]['a'] = G.node[n2]['a']
	for n in G.nodes():
		G_compl.add_node(n)
		G_compl.node[n]['a'] = G.node[n]['a']

	# print G.nodes(data=True)
	# print G_compl.nodes(data=True)

	#print G.edges(data=True)
	#print G_compl.edges(data=True)
	return G, G_compl

def generate_initial_diff_ab():

	G=nx.DiGraph()
	G.add_edges_from([(0,1),(1,2),(2,3),(2,4)])
	G.add_node(5)
	G.node[0]['a'] = 50
	G.node[1]['a'] = 2
	G.node[2]['a'] = 2
	G.node[3]['a'] = 1
	G.node[4]['a'] = 2
	G.node[5]['a'] = 1

	G_compl = nx.DiGraph()
	for n1,n2 in G.edges():
		G_compl.add_edge(n2,n1)
		G_compl.node[n1]['a'] = G.node[n1]['a']
		G_compl.node[n2]['a'] = G.node[n2]['a']
	for n in G.nodes():
		G_compl.add_node(n)
		G_compl.node[n]['a'] = G.node[n]['a']

	# print G.nodes(data=True)
	# print G_compl.nodes(data=True)

	#print G.edges(data=True)
	#print G_compl.edges(data=True)
	return G, G_compl

def generate_simple_loop():

	G=nx.DiGraph()
	G.add_edges_from([(0,1),(1,2),(2,0)])
	G.node[0]['a'] = 1
	G.node[1]['a'] = 1
	G.node[2]['a'] = 1

	G_compl = nx.DiGraph()
	for n1,n2 in G.edges():
		G_compl.add_edge(n2,n1)
		G_compl.node[n1]['a'] = G.node[n1]['a']
		G_compl.node[n2]['a'] = G.node[n2]['a']
	for n in G.nodes():
		G_compl.add_node(n)
		G_compl.node[n]['a'] = G.node[n]['a']

	# print G.nodes(data=True)
	# print G_compl.nodes(data=True)

	#print G.edges(data=True)
	#print G_compl.edges(data=True)
	return G, G_compl

main()
