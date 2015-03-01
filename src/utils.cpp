
#include "utils.h"

int get_reads(const string readFileName, 
	vector<string>& reads
	)
{
	std::istringstream ssFileName(readFileName);
	string tokenFileName, line;

	while(getline(ssFileName, tokenFileName, ',')) 
	{
		try 
		{

			ifstream readFile;
			readFile.open(tokenFileName);

		   	while (getline(readFile , line)) // this is the comment line
		   	{
		    	getline(readFile , line); // the actual read
		    	make_upper_case(line);
		    	reads.push_back(line);
		    	//reads.push_back(reverse_complement(line));
		    	assert(getline(readFile , line)); // the +/- sign
		    	assert(getline(readFile , line)); // the quality values
		   	}
		   	readFile.close();

		   	cout << "Finished reading from " << tokenFileName << endl;
		} catch (exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}

	cout << "Input file(s) contain(s) " << reads.size() << " reads." << endl;

   	return EXIT_SUCCESS;
}


int get_data_for_rlcsa(const string& readFileName, 
	uchar*& data, 
	uint64_t& char_count
	)
{
	std::istringstream ssFileName(readFileName);
	string tokenFileName, line;
	char_count = 0;
	vector<string> reads;

	while(getline(ssFileName, tokenFileName, ',')) 
	{
		try
		{
	    	ifstream readFile;
			readFile.open(tokenFileName);

		   	// counting total number of reads and their length
		   	while (getline(readFile , line)) // this is the comment line
		   	{
		    	getline(readFile , line); // the actual read
		    	make_upper_case(line);
		    	reads.push_back(line);
		    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
		 		// reads.push_back(reverse_complement(line));
		   //  	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'   	
		    	assert(getline(readFile , line)); // the +/- sign
		    	assert(getline(readFile , line)); // the quality values
		   	}
		   	readFile.close();	
		   	cout << "Finished reading from " << tokenFileName << endl;		
		} catch (exception& error) 
		{ // check if there was any error
			std::cerr << "Error: " << error.what() << std::endl;
			return EXIT_FAILURE;
		}
	}


   	cout << "Input file(s) " << readFileName << " contain(s) " << reads.size() << " reads." << endl;
   	cout << "The temporary data array will have size " << (double)char_count/1000000000 << "GB." << endl;

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
   	cout << "Created the data array" << endl;

   	return EXIT_SUCCESS;
}