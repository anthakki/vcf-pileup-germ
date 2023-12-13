
CC = cc -ansi
CXX = c++ -std=c++11
RM = rm -f

CFLAGS += -pedantic -Wall -Wextra
CFLAGS += -DNDEBUG -O3 -fomit-frame-pointer

CFLAGS += `pkg-config --cflags htslib`
LIBS += `pkg-config --libs htslib`

CFLAGS += `pkg-config --cflags zlib`
LIBS += `pkg-config --libs zlib`

# CFLAGS += -g

CXXFLAGS += $(CFLAGS)
CXXFLAGS += -Weffc++

LIBS += -lpthread

targets = vcf-call-fs vcf-pileup vcf-pileup-mt

all: $(targets)

clean:
	$(RM) $(targets)

%: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

fishertest.o: fishertest.c fishertest.h
	$(CC) $(CFLAGS) -c -o $@ $<

vcf-call-fs: vcf-call-fs.cc fishertest.o

vcf-pileup: vcf-pileup.cc

vcf-pileup-mt: vcf-pileup-mt.cc

.PHONY: all clean
.SUFFIXES:
