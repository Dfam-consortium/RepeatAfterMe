#
# RepeatAfterMe
#   A package to extend alignments anchored on a core alignment.  Currently this
#   only includes RAMExtend a tool for treating an existing MSA as an aligned
#   core.
#
#  -Robert Hubley, Sep/2022
#

# Set the version here
VERSION = 0.0.4

# Installation Directory
INSTDIR = /usr/local/RepeatAfterMe-$(VERSION)

LIBFMIDX_DIR = AvxWindowFmIndex
LIBFMIDX_LIBS = ${LIBFMIDX_DIR}/build/libawfmindex.a ${LIBFMIDX_DIR}/build/libdivsufsort64.a
LIBFMIDX_INCDIRS = -I${LIBFMIDX_DIR}/build -I${LIBFMIDX_DIR}/src

CFLAGS = -std=gnu11 -fpic -mtune=native -g -O3 -mavx2 -Wall ${LIBFMIDX_INCDIRS}
LIBS = -lm 
LDFLAGS = 

EXTOBJ = extend_align.o common.o report.o score_system.o cmd_line_opts.o version.o sequence.o bnw_extend.o
EXTINC = extend_align.h common.h report.h score_system.h cmd_line_opts.h sequence.h bnw_extend.h

#BUILD_NUMBER_FILE=build.dat
#BUILD_NUMBER = $$(cat $(BUILD_NUMBER_FILE))
BUILD_NUMBER = 436
BUILD_DATE = $$(date +'%Y%m%d')
BUILD_INSTALL_DIR = rs-$(VERSION)-build-$(BUILD_NUMBER)-$$(date +'%Y%m%d')

all: RAMExtend

RAMExtend: $(EXTOBJ) $(EXTINC) $(BUILD_NUMBER_FILE) kentsrc/libTwoBit.a
	$(CC) $(EXTOBJ) -o $@ $(LIBS) -Lkentsrc -lTwoBit -fopenmp

kentsrc/libTwoBit.a: kentsrc/twoBitNew.c
	$(MAKE) -C kentsrc
	
version.c: Makefile $(BUILD_NUMBER_FILE)
	echo "char const* Version = \"$(VERSION)\";" > version.c
	echo "char const* BuildNumber = \"$(BUILD_NUMBER)\";" >> version.c
	echo "char const* BuildDate = \"$(BUILD_DATE)\";" >> version.c

.c.o:
	$(CC) $(CFLAGS) -c $< -o $*.o $(CCINCLUDES)

clean:
	-@rm *.o *~ RAMExtend
	cd kentsrc; make clean; cd ..

#$(BUILD_NUMBER_FILE): extend_align.o
#	@if ! test -f $(BUILD_NUMBER_FILE); then echo 0 > $(BUILD_NUMBER_FILE); fi
#	@echo $$(($$(cat $(BUILD_NUMBER_FILE)) + 1)) > $(BUILD_NUMBER_FILE)
