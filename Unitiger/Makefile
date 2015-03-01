build: 
	cd ../bin/Unitiger; cmake ../../Unitiger; make; cp ../../Unitiger/src/Unitiger_wrapper.py Unitiger_wrapper.py

fast:
	cd ../bin/Unitiger; g++ ../../Unitiger/src/Unitiger.cpp ../../Unitiger/src/OptionParser.cpp -o Unitiger -I ./ext/gatb-core/include -I ~/bin/boost_1_57_0 -L ./ext/gatb-core/lib -lgatbcore -llibgatbcore -lhdf5 -ldl -lz -lpthread -std=c++0x -O3 -DBOOST_NO_CXX11_RVALUE_REFERENCES=1