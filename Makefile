#
# RepeatAfterMe
#   A package to extend alignments anchored on a core alignment.  Currently this
#   only includes RAMExtend a tool for treating an existing MSA as an aligned
#   core.
#
#  -Robert Hubley, Sep/2022
#
#
#
VERSION = 0.0.6

# Installation Directory
INSTDIR = /usr/local/RepeatAfterMe-$(VERSION)

CFLAGS = -std=gnu11 -fpic -mtune=native -g -O3 -Wall -Iminunit -I.
LIBS = -lm 
LDFLAGS =

EXTOBJ = ram_extend.o common.o report.o score_system.o cmd_line_opts.o version.o sequence.o bnw_extend.o
EXTSRC = ram_extend.c common.c report.c score_system.c cmd_line_opts.c version.c sequence.c bnw_extend.c
EXTINC = ram_extend.h common.h report.h score_system.h cmd_line_opts.h sequence.h bnw_extend.h
TEST_OBJ = cmd_line_opts.o common.o report.o version.o score_system.o sequence.o bnw_extend.o test/test_suite.o    

BUILD_NUMBER_FILE=build.dat
BUILD_NUMBER = $$(cat $(BUILD_NUMBER_FILE))
#BUILD_NUMBER = 436
BUILD_DATE = $$(date +'%Y%m%d')
#BUILD_INSTALL_DIR = RepeatAfterMe-$(VERSION)-build-$(BUILD_NUMBER)-$$(date +'%Y%m%d')

all: RAMExtend test_suite

RAMExtend: $(EXTOBJ) $(EXTINC) kentsrc/libTwoBit.a
	$(CC) $(EXTOBJ) -o $@ $(LIBS) -Lkentsrc -lTwoBit -fopenmp

kentsrc/libTwoBit.a: kentsrc/twoBitNew.c
	$(MAKE) -C kentsrc

test_suite: $(TEST_OBJ) bnw_extend.h 
	$(CC) $(TEST_OBJ) -o $@ $(LIBS) -lpthread -Lkentsrc -lTwoBit $(LIBFMIDX_LIBS) -fopenmp $(LDFLAGS)

version.c: Makefile $(BUILD_NUMBER_FILE)
	echo "char const* Version = \"$(VERSION)\";" > version.c
	echo "char const* BuildNumber = \"$(BUILD_NUMBER)\";" >> version.c
	echo "char const* BuildDate = \"$(BUILD_DATE)\";" >> version.c

# Build number file.  Increment if any object file changes.
EXTSRC_SANS_VERSION := $(filter-out version.c, $(EXTSRC))
$(BUILD_NUMBER_FILE): $(EXTSRC_SANS_VERSION) $(EXTINC)
	@if ! test -f $(BUILD_NUMBER_FILE); then echo 0 > $(BUILD_NUMBER_FILE); fi
	@echo $$(($$(cat $(BUILD_NUMBER_FILE)) + 1)) > $(BUILD_NUMBER_FILE)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o $(CCINCLUDES)

beautify:
	indent -bap -cdb -bl -bli0 -npcs -nut -lp build_kmer_table.c
	indent -bap -cdb -bl -bli0 -npcs -nut -lp build_repeat_families.c
	indent -bap -cdb -bl -bli0 -npcs -nut -lp cmd_line_opts.c
	indent -bap -bl -bli0 -npcs -nut -lp ram_extend.c

try:
	./RAMExtend -twobit test/extension-test2.2bit -ranges test/extension-test2.tsv

try2:
	./RAMExtend -twobit ../satellite-1/genome.2bit -ranges ../satellite-1/linup.tsv -matrix 14p43g

install: all
	@mkdir $(INSTDIR)
	@mkdir $(INSTDIR)/bin
	cp RAMExtend $(INSTDIR)/bin
	cp README.md $(INSTDIR)

clean:
	-@rm *.o *~ RAMExtend test_suite
	cd kentsrc; make clean; cd ..
