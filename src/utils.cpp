
#include "utils.h"

// int get_reads(const string readFileName, 
// 	vector<string>& reads
// 	)
// {
// 	std::istringstream ssFileName(readFileName);
// 	string tokenFileName, line;

// 	while(getline(ssFileName, tokenFileName, ',')) 
// 	{
// 		try 
// 		{
// 			ifstream readFile;
// 			readFile.open(tokenFileName);

// 			cout << "*** Reading from file '" << tokenFileName << "'" << endl;
// 		   	while (getline(readFile , line)) // this is the comment line
// 		   	{
// 		    	getline(readFile , line); // the actual read
// 		    	make_upper_case(line);
// 		    	reads.push_back(line);
// 		    	//reads.push_back(reverse_complement(line));
// 		    	assert(getline(readFile , line)); // the +/- sign
// 		    	assert(getline(readFile , line)); // the quality values
// 		   	}
// 		   	readFile.close();

// 		} catch (exception& error) 
// 		{ // check if there was any error
// 			std::cerr << "Error: " << error.what() << std::endl;
// 			return EXIT_FAILURE;
// 		}
// 	}

// 	cout << "*** Input file(s) contain(s) " << reads.size() << " reads." << endl;

//    	return EXIT_SUCCESS;
// }

int get_reads_using_Bank(const string readFileName, 
	vector<string>& reads
	)
{
	Bank *Reads = new Bank(const_cast<char*>(readFileName.c_str()));
	int readlen;
	char *rseq;
	string line;

	while( Reads->get_next_seq(&rseq,&readlen) )
    {
    	line = rseq;
        make_upper_case(line);
        reads.push_back(line);
    }

	cout << "*** The file(s) listed in " << readFileName << " contain(s) " << reads.size() << " reads." << endl;

	delete Reads;

   	return EXIT_SUCCESS;
}


// int get_data_for_rlcsa(const string& readFileName, 
// 	uchar*& data, 
// 	uint64_t& char_count
// 	)
// {
// 	std::istringstream ssFileName(readFileName);
// 	string tokenFileName, line;
// 	char_count = 0;
// 	vector<string> reads;

// 	while(getline(ssFileName, tokenFileName, ',')) 
// 	{
// 		try
// 		{
// 	    	ifstream readFile;
// 			readFile.open(tokenFileName);

// 			cout << "*** Reading from file '" << tokenFileName << "'" << endl;
// 		   	// counting total number of reads and their length
// 		   	while (getline(readFile , line)) // this is the comment line
// 		   	{
// 		    	getline(readFile , line); // the actual read
// 		    	make_upper_case(line);
// 		    	reads.push_back(line);
// 		    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
// 		 		// reads.push_back(reverse_complement(line));
// 		   		// char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'   	
// 		    	assert(getline(readFile , line)); // the +/- sign
// 		    	assert(getline(readFile , line)); // the quality values
// 		   	}
// 		   	readFile.close();	
// 		} catch (exception& error) 
// 		{ // check if there was any error
// 			std::cerr << "Error: " << error.what() << std::endl;
// 			return EXIT_FAILURE;
// 		}
// 	}

//    	cout << "*** Input file(s) " << readFileName << " contain(s) " << reads.size() << " reads." << endl;
//    	cout << "*** The temporary data array will have size " << (double)char_count/1000000000 << "GB." << endl;

//    	uint64_t i = 0;
// 	data = new uchar[char_count];
// 	for (auto read : reads)
// 	{
//     	for (uint64_t j = 0; j < read.length(); j++)
//     	{
//     		data[i] = (uchar)read[j];
//     		i++;
//     	}
//     	data[i] = '\0';
//     	i++;		
// 	}

//    	return EXIT_SUCCESS;
// }

int get_data_for_rlcsa_using_Bank(const string& readFileName, 
	uchar*& data, 
	uint64_t& char_count
	)
{
	Bank *Reads = new Bank(const_cast<char*>(readFileName.c_str()));
	int readlen;
	char *rseq;
	string line;
	vector<string> reads;
	char_count = 0;	

	while( Reads->get_next_seq(&rseq,&readlen) )
    {
    	line = rseq;
    	make_upper_case(line);
    	reads.push_back(line);
    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
    }

   	cout << "*** The file(s) listed in " << readFileName << " contain(s) " << reads.size() << " reads." << endl;
   	cout << "*** The temporary data array will have size " << (double)char_count/1000000000 << "GB." << endl;

   	uint64_t i = 0;
	data = new uchar[char_count];
	for (auto read : reads)
	{
    	for (uint64_t j = 0; j < read.length(); j++)
    	{
    		data[i] = (uchar)read[j];
    		i++;
    	}
    	data[i] = '\0';
    	i++;		
	}

	delete Reads;

   	return EXIT_SUCCESS;
}