CXXFLAGS=-O3 -Wall -Wextra -pedantic -pthread -lboost_system -std=gnu++17 -march=native -mtune=native

#g++ square.cpp -O3 -Wall -pthread -lboost_system -std=c++17 -o square_mt


.PHONY: all clean masks.hpp

all: square_mt

masks.hpp: masks
	./masks > masks.hpp

masks: masks.cpp config.h utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

square_mt: square.cpp masks.hpp gf2_monomial.cpp stringops.cpp parallel_hashmap/phmap.h config.h utils.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

clean:
	@rm -rf square_mt masks masks.hpp
