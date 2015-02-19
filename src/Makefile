CC = g++
FILES = main.cpp OptionParser.cpp
OUT_EXE = optimal-k

FILES_VALIDATE = validate.cpp OptionParser.cpp
OUT_EXE_VALIDATE = validate

CPPFLAGS = -std=c++0x -O3 -fopenmp -DMASSIVE_DATA_RLCSA
# -Wall

build: $(FILES)
	cd rlcsa; make; cd ..; $(CC) $(CPPFLAGS) -o bin/$(OUT_EXE) $(FILES) rlcsa/rlcsa.a

validate: $(FILES_VALIDATE)
	cd rlcsa; make; cd ..; $(CC) $(CPPFLAGS) -o $(OUT_EXE_VALIDATE) $(FILES_VALIDATE) rlcsa/rlcsa.a

clean:
	rm -f $(OUT_EXE)

rebuild: clean build