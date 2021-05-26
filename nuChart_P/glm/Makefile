CXX=g++
CXXFLAGS=-Wall -Wextra -g

all: GLM.o Fit.o
	$(CXX) $(CXXFLAGS) GLM.o Fit.o #-lm -lgsl -lgslcblas -L /usr/local/lib
	
GLM.o: GLM.cpp GLM.hpp
	$(CXX) $(CXXFLAGS) -c GLM.cpp
Fit.o: Fit.cpp Fit.hpp GLM.hpp
	$(CXX) $(CXXFLAGS) -c Fit.cpp
clean:
	rm -f *.o *~
