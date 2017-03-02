CC=gcc
CFLAGS=-m32 -O2 -Wall
LDLIBS=-lm

EXEC=mk_MatSmoothing mk_MatPaths


all: $(EXEC)

mk_MatSmoothing: mk_MatSmoothing.o

mk_MatPaths: mk_MatPaths.o


.PHONY: clean cleanall save


cleanall: clean
	rm -f $(EXEC)
	

clean:
	rm -f *.o 


SAVEFILE=ANSWTOMO_TB_$(shell date +'%Y%m%d').tgz
save:
	tar -cvzf $(SAVEFILE) *.m *.c Makefile DATAFiles m_map
