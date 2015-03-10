# optimal_k
Calculates optimal k for a DBG assembler given one or more read libraries

# Installation for optimal-k

In *src directory*, run 

	make

This places the executable *optimal-k* in the directory *bin*.

# Installation for Unitiger

In *Unitiger* directory, run

	make

The places the executable *Unitiger* and the Python wrapper 
*Unitiger_wrapper.py* in the directory *bin/Unitiger*

# Usage

	optimal-k OPTIONS

	python Unitiger_wrapper.py OPTIONS

# Simple example

Suppose you are in data/human/

	../../bin/optimal-k -r frag_1.fastq,frag_2.fastq -o metrics_file

	python ../../bin/Unitiger/Unitiger_wrapper.py -r frag_1.fastq,frag_2.fastq -o metrics_file

# Test run

We force the program to build the index and then exit (notice that we set k > K)

	../../bin/optimal-k -r frag_1.fastq,frag_2.fastq -o metrics_file -k 2 -K 1 -b

We now run it for more abunances (the index is now loaded from disk)

	../../bin/optimal-k -r frag_1.fastq,frag_2.fastq -o metrics_file -a 1 -A 5