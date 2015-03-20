
#include "utils.h"

int get_reads(const string readFileName, 
	vector<compact_read>& reads,
	uint64_t &reads_total_content,
	uint64_t &reads_number,
	uint32_t &reads_max_length,
	uint32_t &reads_min_length
	)
{
	double start_reading = readTimer();
	Bank *reads_bank = new Bank(const_cast<char*>(readFileName.c_str()));
	int readlen;
	char *rseq;
	string line;
	uint64_t estimated_nb_reads = reads_bank->estimate_nb_reads();

	random_device rd;
	default_random_engine generator(rd());
	uniform_real_distribution<double> uniform_read_distribution(0.0,1.0);

   	try
   	{
   		reads.reserve(estimated_nb_reads * READ_PROPORTION);
   	}
   	catch (exception& error) 
	{
		cout << "*** ERROR: could not allocate memory for storing reads:" << endl;
		cout << "***    " << error.what() << std::endl;
		return EXIT_FAILURE;
	}
	cout << "*** Reserved memory for a sample of " << reads.capacity() << " reads " << endl;
	
	reads_total_content = 0;
	reads_max_length = 0;
	reads_min_length = 1000000;
	reads_number = 0;

	while( reads_bank->get_next_seq(&rseq,&readlen) )
    {
    	line = rseq;
        // not needed because the encoding make the read uppercase
        //make_upper_case(line);
        reads_total_content += line.length();
        reads_number++;
        if (line.length() > reads_max_length)
        {
        	reads_max_length = line.length();
        }
        if (line.length() < reads_min_length)
        {
        	reads_min_length = line.length();
        }
        if (uniform_read_distribution(generator) < READ_PROPORTION)
    	{
    		reads.push_back(encode_string(line));	
    	}
    }

	cout << "*** The file(s) listed in " << readFileName << " contain(s) " << reads_number << " reads" << endl;
	cout << "*** They were read in " << readTimer() - start_reading << " seconds" << endl;

	delete reads_bank;

   	return EXIT_SUCCESS;
}

int get_data_and_build_rlcsa_noniterative(const string& readFileName, 
	const string& indexFileName,
	const uint N_THREADS
	)
{
	Bank *reads_bank = new Bank(const_cast<char*>(readFileName.c_str()));
	int readlen;
	char *rseq;
	string line;
	vector<string> reads;
	uint64_t char_count = 0;
	uchar *data = NULL;

	cout << "*** Building the RLCSA index on the reads and saving it to files:" << endl;
	cout << "***    " << indexFileName + ".rlcsa.array" << endl;
	cout << "***    " << indexFileName + ".rlcsa.parameters" << endl;

	cout << "*** Initialized the Bank object from the file(s) listed in " << readFileName << endl;

	// getting the size of the temporary data array
	while( reads_bank->get_next_seq(&rseq,&readlen) )
    {
    	line = rseq;
    	make_upper_case(line);
    	reads.push_back(line);
    	char_count = char_count + (line.length() + 1); // +1 for the \0 which will terminate each read in uchar array 'data'
    }

    delete reads_bank;

   	cout << "*** The file(s) listed in " << readFileName << " contain(s) " << reads.size() << " reads." << endl;
   	cout << "*** The temporary data array will have size " << (double)char_count/1000000000 << "GB." << endl;

   	uint64_t i = 0;
   	try
   	{
   		data = new uchar[char_count];
   	}
   	catch (exception& error) 
	{
		cout << "*** ERROR: could not allocate memory for the temporary data array:" << endl;
		std::cerr << "***    " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

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
	cout << "*** Successfully created the temporary data array" << endl;

	cout << "*** Now creating the index" << endl;

	try
	{
		RLCSA rlcsa_built(data, char_count, 32, 0, N_THREADS, true);
 		data = 0; // The constructor deleted the data.

		if ( !(rlcsa_built.isOk()) ) 
 		{
 			cout << "*** ERROR: could not create the index" << endl;
 			return EXIT_FAILURE;
 		}
 		rlcsa_built.printInfo();
 		rlcsa_built.reportSize(true);
 		try
		{
 		rlcsa_built.writeTo(indexFileName);	
		}
		catch (exception& error)
		{
			cout << "*** ERROR: could not write the index to file:" << endl;
			std::cerr << "***    " << error.what() << std::endl;
			return EXIT_FAILURE;
		}	 		
	}
	catch (exception& error) 
	{ // check if there was any error
		std::cerr << "*** ERROR: " << error.what() << std::endl;
		return EXIT_FAILURE;
	}


   	return EXIT_SUCCESS;
}

int get_data_and_build_rlcsa_iterative(const string& readFileName, 
	const string& indexFileName,
	const uint N_THREADS,
	const bool lower_memory_construction
	)
{
	uint64_t DATA_SIZE;
	if (lower_memory_construction)
	{
		DATA_SIZE = 1000;
	}
	else
	{
		DATA_SIZE = 3999;
	}
	uint64_t INSERT_SIZE = DATA_SIZE - 1;

	Bank *reads_bank = new Bank(const_cast<char*>(readFileName.c_str()));
	int readlen;
	char *read_seq;
	string line;
	char* data = NULL;
	uint64_t char_count, seq_count;	

	cout << "*** Building the RLCSA index on the reads and saving it to files:" << endl;
	cout << "***    " << indexFileName + ".rlcsa.array" << endl;
	cout << "***    " << indexFileName + ".rlcsa.parameters" << endl;

	cout << "*** Initialized the Bank object from the file(s) listed in " << readFileName << endl;

	try
   	{
   		data = new char[DATA_SIZE * MEGABYTE];
   	}
   	catch (exception& error) 
	{
		cerr << "*** ERROR: could not allocate memory for the temporary data array:" << endl;
		cerr << "***    " << error.what() << std::endl;
		return EXIT_FAILURE;
	}

  	RLCSABuilder builder(32, 0, DATA_SIZE * MEGABYTE, N_THREADS); 

	// getting the size of the temporary data array
	char_count = 0;
	seq_count = 0;
	uint64_t n_insertions = 0;
	while( reads_bank->get_next_seq(&read_seq,&readlen) )
    {
    	for (int i = 0; i < readlen; i++)
    	{
    		data[char_count] = std::toupper(read_seq[i]);
    		char_count++;
    		seq_count++;
    	}

    	data[char_count] = '\0';
    	char_count++; // +1 for the \0 which will terminate each read in uchar array 'data'
		
		if (char_count > INSERT_SIZE * MEGABYTE)
		{
			n_insertions++;
  			cout << "*** "<< n_insertions << ": Inserting " << (double)seq_count / MEGABYTE << "MB of sequence into the RLCSA index" << endl;
  			builder.insertSequence(data, char_count - 1, false); // -1 because the last \0 should not count
  			char_count = 0;
  			seq_count = 0;
		}
    }
    // inserting the remaining sequence
	if (char_count > 2 * MEGABYTE)	
	{
		n_insertions++;
		cout << "*** "<< n_insertions << ": Inserting " << (double)seq_count / MEGABYTE << "MB of sequence into the RLCSA index" << endl;
		// Insert the sequence:
		builder.insertSequence(data, char_count - 1, false); // -1 because the last \0 should not count
		char_count = 0;
	}

	cout << "*** Now trying to write the index to disk" << endl;

    // delete reads_bank;
    // delete[] data;

    // If successful, write the index to disk.
	if (builder.isOk())
	{
		RLCSA* rlcsa = builder.getRLCSA();
		rlcsa->printInfo();
	 	rlcsa->reportSize(true);
		rlcsa->writeTo(indexFileName);
	 	cout << "*** Created the RLCSA index and saved it to disk" << endl;
		// delete rlcsa;
	}
	else
	{
		cerr << "*** ERROR: the RLCSA index could not be created" << endl;
		return EXIT_FAILURE;
	}

   	return EXIT_SUCCESS;
}

