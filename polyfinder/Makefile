CXXFLAGS=-O3 -Wall -Wextra -pedantic -pthread -lboost_system -std=c++17 -march=native -mtune=native

#g++ square.cpp -O3 -Wall -pthread -lboost_system -std=c++17 -o square_mt


.PHONY: all clean

all: square_mt

square_mt: square.cpp gf2_monomial.cpp common.cpp stringops.cpp parallel_hashmap/phmap.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

clean:
	@rm -rf square_mt
