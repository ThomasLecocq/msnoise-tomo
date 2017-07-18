
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

vg_fta.so: configparser.cpp fft_NR.cpp fta_param.cpp libfta.cpp readsac.cpp vg_fta.cpp
    g++ -c configparser.cpp
    g++ -c fft_NR.cpp
    g++ -c fta_param.cpp
    g++ -c libfta.cpp
    g++ -c readsac.cpp
    g++ -c vg_fta.cpp
    g++ vg_fta.cpp -shared -o vg_fta.so libfta.o readsac.o configparser.o fta_param.o fft_NR.o -lm