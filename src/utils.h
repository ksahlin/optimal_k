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
#include <string>
#include <random>
#include <zlib.h>

#include "rlcsa/rlcsa.h"
#include "rlcsa/rlcsa_builder.h"
// #include "rlcsa/misc/utils.h"

#include "minia/Bank.h"
#include "OptionParser.h"

using namespace CSA;
using namespace std;

#define MIN(a,b) (a <= b) ? a : b
#define MAX(a,b) (a <= b) ? b : a

#define MEGABYTE 1000000

inline char reverse_complement_char(char c)
{
    if (c == 'A') return 'T';
    if (c == 'T') return 'A';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    return c;
}

inline string reverse_complement(const string& s)
{
    string reverse = s;

    for (int i = s.length()-1; i >= 0; i--)
    {
        reverse[s.length() - 1 - i] = reverse_complement_char(s[i]);
    }

    return reverse;
}

inline int calc_abundance(const RLCSA* rlcsa, 
    const string& sample
    )
{
    pair_type result = rlcsa->count(sample);
    pair_type result_rc = rlcsa->count(reverse_complement(sample));
    int abundance = length(result) + length(result_rc);
    return abundance;
    // return length(result);
}

inline string int_to_string(size_t x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}

inline void make_upper_case(string& s)
{
    for (uint i = 0; i < s.length(); i++)
    {
        s[i] = toupper(s[i]);
    }
}

inline string get_first_token(string s)
{
    istringstream ss(s);
    string token;

    getline(ss, token, ',');

    return token;
}

// int get_reads(const string readFileName, 
// 	vector<string>& reads
// 	);

int get_reads_using_Bank(const string readFileName, 
    vector<string>& reads
    );

// int get_data_for_rlcsa(const string& readFileName, 
//     uchar*& data, 
//     uint64_t& char_count
//     );

int get_data_for_rlcsa_using_Bank(const string& readFileName, 
	uchar*& data, 
	uint64_t& char_count
	);

int get_data_and_build_rlcsa_using_Bank(const string& readFileName, 
    const string& indexFileName,
    const uint N_THREADS
    );

// Check if a file is readable
inline bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 

#endif // UTILS_H_INCLUDED