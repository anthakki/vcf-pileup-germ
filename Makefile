
CC = cc -std=c99
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

targets := $(basename $(wildcard *.cc))
objects := fishertest.o gtbeta.o jacobi_rule.o $(addsuffix .o,$(targets))

all: $(targets)

cleanobj:
	$(RM) $(objects)

clean: cleanobj
	$(RM) $(targets)

%: %.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

fishertest.o: fishertest.c fishertest.h
	$(CC) $(CFLAGS) -c -o $@ $<

gtbeta.o: gtbeta.c gtbeta.h
	$(CC) $(CFLAGS) -c -o $@ $<

jacobi_rule.o: jacobi_rule.c jacobi_rule.h
	$(CC) $(CFLAGS) -c -o $@ $<

vcf-call-fs.o: fishertest.h gzfile.hh string_view.hh vcfreader.hh
vcf-call-fs: vcf-call-fs.o fishertest.o

vcf-call-germ.o: gtbeta.h vcfreader.hh
vcf-call-germ: vcf-call-germ.o gtbeta.o jacobi_rule.o

vcf-pileup.o: gzfile.hh string_view.hh ordereddict.hh tempfile.hh vcfreader.hh
vcf-pileup: vcf-pileup.o

vcf-pileup-mt.o: gzfile.hh string_view.hh vcfreader.hh
vcf-pileup-mt: vcf-pileup-mt.o

vcf-sum-format.o: vcfreader.hh
vcf-sum-format: vcf-sum-format.o

.PHONY: all cleanobj clean
.SUFFIXES:
