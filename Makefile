# adjust to your needs:
CXX=clang++-4.0
EIGEN = /usr/include/eigen3
MKLROOT = /opt/intel/mkl


CXXFLAGS = -O3 -Wall -pedantic -std=c++11
MKLINCLUDE = -I${MKLROOT}/include
MKL = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl


all:
	$(CXX) $(CXXFLAGS) src/aw2giw.cpp -o aw2giw.o
	$(CXX) $(CXXFLAGS) src/aw2gtau.cpp -o aw2gtau.o
	$(CXX) $(CXXFLAGS) $(MKLINCLUDE) -I$(EIGEN) $(MKL) src/hk2gwfw.cpp -o hk2gwfw.o

giw:
	$(CXX) $(CXXFLAGS) src/aw2giw.cpp -o aw2giw.o

gtau:
	$(CXX) $(CXXFLAGS) src/aw2gtau.cpp -o aw2gtau.o

gwfw:
	$(CXX) $(CXXFLAGS) $(MKLINCLUDE) -I$(EIGEN) $(MKL) src/hk2gwfw.cpp -o hk2gwfw.o

clean:
	rm -f *.o
