CC = g++
FILES = main.cpp esize_estimation.h OptionParser.cpp utils.cpp minia/Bank.cpp
OUT_EXE = optimal-k

FILES_VALIDATE = validate.cpp OptionParser.cpp utils.cpp
OUT_EXE_VALIDATE = validate

CPPFLAGS = -std=c++0x -O3 -fopenmp -DMASSIVE_DATA_RLCSA
# -Wall

COMPILE_RLCSA = rm -rf ../bin/rlcsa; mkdir -p ../bin/rlcsa; cp -r rlcsa/ ../bin/; cd ../bin/rlcsa; make build_rlcsa

build: $(FILES)
	mkdir -p ../bin; $(COMPILE_RLCSA); cd ../../src/; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) ../bin/rlcsa/rlcsa.a -Wall -lz

quick: $(FILES)
	mkdir -p ../bin; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) ../bin/rlcsa/rlcsa.a -Wall -lz

# validate: $(FILES_VALIDATE)
# 	cd ../bin/rlcsa; make build_rlcsa; cd ../../src/; mkdir ../bin; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE_VALIDATE) $(FILES_VALIDATE) ../bin/rlcsa/rlcsa.a

clean:
	rm -f ../bin/$(OUT_EXE); rm -rf ../bin/rlcsa;

rebuild: clean build