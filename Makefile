
CXX = c++ -std=c++11
RM = rm -f

CXXFLAGS += -pedantic -Wall -Wextra -Weffc++
CXXFLAGS += -DNDEBUG -O3 -fomit-frame-pointer

CXXFLAGS += `pkg-config --cflags htslib`
LIBS += `pkg-config --libs htslib`

CXXFLAGS += `pkg-config --cflags zlib`
LIBS += `pkg-config --libs zlib`

CXXFLAGS += -g

all: vcf-pileup

clean:
	$(RM) vcf-pileup

%: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

.PHONY: all clean
.SUFFIXES:
