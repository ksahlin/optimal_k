# optimal_k
Calculates optimal k for a DBG assembler given one or more read libraries

# Installation for optimal_k

In src directory, run 

	make

This places the executable *optimal-k* in the directory *bin*.

In Unitiger directory, run

	make

The places the executable *Unitiger* and the 
Python wrapper *Python_wrapper.py* in the directory *bin/Unitiger*

# Usage

Run 

	optimal-k OPTIONS

	python Unitiger_wrapper.py OPTIONS

# Simple example

	optimal-k -r frag_1.fastq,frag_2.fastq -o metrics_file

	python Unitiger_wrapper.py -r frag_1.fastq,frag_2.fastq -o metrics_file