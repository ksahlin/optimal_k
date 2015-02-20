#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <fstream>
#include <iostream>
#include <vector>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>
#include <omp.h>
#include <algorithm>

#include "../bin/rlcsa/rlcsa.h"
#include "OptionParser.h"

using namespace CSA;
using namespace std;

#define MIN(a,b) (a <= b) ? a : b

inline char reverse_complement_char(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

inline string reverse_complement(string& s)
{
    string reverse;

    for (int i = s.length()-1; i >= 0; i--)
    {
        reverse += reverse_complement_char(s[i]);
    }

    return reverse;
}

inline string int_to_string(size_t x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}

int get_reads(const string readFileName, 
	vector<string>& reads
	);


int get_data_for_rlcsa(const string& readFileName, 
	uchar*& data, 
	uint64_t& char_count
	);
#endif // UTILS_H_INCLUDED