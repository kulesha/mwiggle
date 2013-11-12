mwiggle: lib mwiggle.c
	LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
	gcc -L. -Wall -lmw mwiggle.c -o mwiggle

lib: mw.o
	gcc -shared -Wl,-soname,libmw.so -o libmw.so *.o

mw.o: mw.c mw.h
	gcc -Wall -fPIC -c mw.c

clean: 
	rm *.o *.so mwiggle
