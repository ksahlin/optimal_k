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
#include <math.h>

#include "rlcsa/rlcsa.h"
#include "rlcsa/rlcsa_builder.h"
#include "rlcsa/misc/definitions.h"

#include "minia/Bank.h"
#include "OptionParser.h"

#include "lut.h"

using namespace CSA;
using namespace std;

#define MIN(a,b) (a <= b) ? a : b
#define MAX(a,b) (a <= b) ? b : a

struct compact_read 
{
    unsigned char *read; // the encoded read
    uint16_t length; // the length of the original read
};

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

inline compact_read encode_string(const string& s)
{
    compact_read cread;
    cread.read = new unsigned char[(uint32_t)ceil(s.length()/3)];
    cread.length = s.length();

    char first, second, third;
    uint32_t cread_index = 0;
    for (uint32_t i = 0; i < s.length(); i = i + 3)
    {
        first = base_to_number[(int)s[i]];
        if (i + 1 < s.length())
        {
            second = base_to_number[(int)s[i + 1]];
            if (i + 2 < s.length())
            {
                third = base_to_number[(int)s[i + 2]];
            }
            else
            {
                third = 0;
            }
        }
        else
        {
            second = 0;
            third = 0;
        }
        cread.read[cread_index] = 25 * first + 5 * second + third;
        cread_index++; 
    }

    return cread;
}

inline string decode_substring(const compact_read& cread, const uint32_t &start, const uint32_t &length)
{

    // This conditions must be satisfied!!!
    // assert(start + length <= cread.length);
    // assert(length >= 3);

    string read = "";
    
    uint32_t start_triplet = start / 3;  // i.e. floor(start / 3)
    uint32_t end_triplet = (start + length - 1) / 3;  // i.e. floor((start + length - 1) / 3)

    // adding the characters from the first triplet
    for (uint32_t i = start % 3; i < 3; i++)
    {
        read += number_to_basetriplet[cread.read[start_triplet]][i];
    }

    // adding all the characters from the middle triplets
    for (uint32_t triplet = start_triplet + 1; triplet < end_triplet; triplet++)
    {
        for (uint32_t i = 0; i < 3; i++)
        {
            read += number_to_basetriplet[cread.read[triplet]][i];
        }
    }

    // adding the characters from the last triplet
    for (uint32_t i = 0; i <= (start + length - 1) % 3; i++)
    {
        read += number_to_basetriplet[cread.read[end_triplet]][i];
    }    

    return read;
}

int get_reads(const string readFileName, 
    vector<compact_read>& reads,
    uint64_t &reads_total_content,
    uint32_t &reads_max_length,
    uint32_t &reads_min_length
    );

int get_data_and_build_rlcsa_noniterative(const string& readFileName, 
    const string& indexFileName,
    const uint N_THREADS
    );

int get_data_and_build_rlcsa_iterative(const string& readFileName, 
    const string& indexFileName,
    const uint N_THREADS,
    const bool lower_memory_construction
    );

// Check if a file is readable
inline bool is_readable( const std::string & file ) 
{ 
    std::ifstream f( file.c_str() ); 
    return !f.fail(); 
} 




#endif // UTILS_H_INCLUDED