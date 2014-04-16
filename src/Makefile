exe:	lib
	gcc -L../lib -Wall -lmw mwiggle.c -o ../bin/mwiggle

lib: 	mw.o
	gcc -shared -Wl,-soname,libmw.so -o ../lib/libmw.so *.o

mw.o: 	mw.c mw.h
	gcc -Wall -fPIC -c mw.c

clean: 
	@rm -f *.o 
