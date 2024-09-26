all: serial life-nonblocking

CXX = g++
MPICXX = mpicxx
CXXFLAGS = -std=c++11 -O2

serial: serial.cpp
	$(CXX) $(CXXFLAGS) -o life $<

life-nonblocking: life-nonblocking.cpp
	$(CXX) $(CXXFLAGS) -o life-nonblocking $<

clean:
	rm -f life life-nonblocking
