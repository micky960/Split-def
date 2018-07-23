CC=g++ 
CXXFLAGS=-std=c++11 -g -O3 -fopenmp
DEPS = lefparser.h
OBJ = split_def.o lefparser.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS)

split_def: $(OBJ)
	$(CC) $(CXXFLAGS) -o split_def split_def.o lefparser.o -I.

clean:
	rm -f $(binaries) *.o

