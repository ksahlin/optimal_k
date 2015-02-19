CC = g++
FILES = main.cpp OptionParser.cpp
OUT_EXE = optimal-k

FILES_VALIDATE = validate.cpp OptionParser.cpp
OUT_EXE_VALIDATE = validate

CPPFLAGS = -std=c++0x -O3 -fopenmp -DMASSIVE_DATA_RLCSA
# -Wall

build: $(FILES)
	cd rlcsa; make build_rlcsa; cd ..; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE) $(FILES) ../bin/rlcsa/rlcsa.a

validate: $(FILES_VALIDATE)
	cd rlcsa; make build_rlcsa; cd ..; $(CC) $(CPPFLAGS) -o ../bin/$(OUT_EXE_VALIDATE) $(FILES_VALIDATE) ../bin/rlcsa/rlcsa.a

clean:
	rm -f $(OUT_EXE); cd ../bin/rlcsa; make clean

rebuild: clean build