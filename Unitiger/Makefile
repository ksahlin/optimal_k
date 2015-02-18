CC = g++
FILES = main.cpp OptionParser.cpp
OUT_EXE = unitiger

CPPFLAGS = -std=c++0x -O3 -Igatb-core/include -I ~/bin/boost_1_57_0 -Lgatb-core/lib -lgatbcore -lhdf5 -ldl -lz -lpthread -DBOOST_NO_CXX11_RVALUE_REFERENCES=1 -stdlib=libstdc++
# -Wall

build: $(FILES)
	$(CC) $(CPPFLAGS) -o $(OUT_EXE) $(FILES)

clean:
	rm -f $(OUT_EXE)

rebuild: clean build