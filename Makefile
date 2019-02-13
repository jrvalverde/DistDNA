# your C compiler
CC = cc

# use this for 64 bit machines and GCC
# CFLAGS = -O -m64
#and this for agnostic builds
CFLAGS = -O

LDFLAGS = -lm

# default is debugging executable
ddna:	distdna.c
	$(CC) -g -o ddna distdna.c $(LDFLAGS)

distdna:	distdna.c
	$(CC) $(CFLAGS) -o distdna distdna.c $(LDFLAGS)

