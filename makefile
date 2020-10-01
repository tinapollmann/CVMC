CC = c++ -Wall -std=c++11
FLAGS = $$(root-config --cflags --libs)
HEADERS = $(wildcard *.h) $(wildcard *.hpp)
CXXES = $(wildcard *.cxx)
EXE = runCVMC

$(EXE) : $(CXXES) $(HEADERS)
	$(CC) $(FLAGS) $(CXXES) -o $@

clean:
	rm -f $(EXE)
	rm -f *.dot
	rm -f *.txt
	rm -f *.gnu
	rm -f *.root
	rm -f *.png	
	
.PHONY: clean
