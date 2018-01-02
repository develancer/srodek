CXX = mpic++
CXXFLAGS = -W -Wall -Wextra -O2 --std=c++11
LDLIBS = -lm -lshp -lGeographic

.PHONY : default

default : srodek
