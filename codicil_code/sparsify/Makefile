include Makefile.in

INCLUDES = -I. 
CFLAGS = $(COPTIONS) $(INCLUDES)
LIBS = -lm 

LD = $(CC) $(LDOPTIONS) 

default: sparsify
#default: sparsify_debug

sparsify_debug: *.c *.h
	$(CC) $(CFLAGS) $(GDBFLAGS) -c *.c
	$(LD) $(LIBS) *.o -o $@
	chmod 755 $@

sparsify: *.c *.h
	$(CC) $(CFLAGS) $(OPTFLAGS) -c *.c
	$(LD) $(LIBS) *.o -o $@
	chmod 755 $@

clean:
	rm -f *.o sparsify sparsify_debug
