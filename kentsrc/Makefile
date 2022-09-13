
CC=gcc
CFLAGS=-I.


SRCS := twoBitNew.c common.c memalloc.c hash.c errAbort.c obscure.c localmem.c dnautil.c dlist.c dystring.c linefile.c options.c verbose.c fa.c dnaseq.c bits.c
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

%.o: %.c
	$(CC) $(CPPFLAGS) -g $(CFLAGS) -c $< -o $@

all: $(SRCS) libTwoBit.a twoBitReadTest

# Generate a static lib
libTwoBit.a: twoBitNew.o common.o memalloc.o hash.o errAbort.o obscure.o localmem.o dnautil.o dlist.o dystring.o linefile.o options.o verbose.o fa.o dnaseq.o bits.o
	ar rcs libTwoBit.a dystring.o dlist.o dnautil.o localmem.o obscure.o errAbort.o hash.o memalloc.o common.o twoBitNew.o linefile.o options.o verbose.o fa.o dnaseq.o bits.o
#	ranlib libTwoBit.a

twoBitReadTest: twoBitReadTest.o libTwoBit.a
	$(CC) -lm -o twoBitReadTest twoBitReadTest.o -L. -lTwoBit

twoBitToFa: twoBitToFa.o libTwoBit.a
	$(CC) -lm -o twoBitToFa twoBitToFa.o -L. -lTwoBit

clean:
	-@rm *.o *~ *.a twoBitReadTest twoBitToFa
