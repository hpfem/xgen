CC=g++
CFLAGS= -O2 -I/usr/X11R6/include
LDFLAGS= -L/usr/X11R6/lib -lXm -lX11 -lm

init: xgen
	@echo xgen ready.

clean:
	/bin/rm -f obj/*.o *~ */*~

xgen: obj/disc.o obj/xgen.o obj/xmain.o
	$(CC) -o xgen obj/xmain.o obj/disc.o obj/xgen.o $(LDFLAGS)

obj/disc.o: src/disc.cpp src/disc.h
	$(CC) -o obj/disc.o -c src/disc.cpp $(CFLAGS)

obj/xgen.o: src/xgen.cpp src/xgen.h
	$(CC) -o obj/xgen.o -c src/xgen.cpp $(CFLAGS)

obj/xmain.o: src/xmain.cpp src/disc.cpp src/xgen.cpp
	$(CC) -o obj/xmain.o -c src/xmain.cpp $(CFLAGS)


