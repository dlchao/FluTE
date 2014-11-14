### make changes accordingly ###
CC       = gcc
CPP      = g++
CLINKER  = gcc
CCLINKER = g++
MPICC       = mpicc
MPICPP      = mpicxx
MPICCLINKER = mpicxx
MAKE     = make --no-print-directory
SHELL    = /bin/sh
CFLAGS		= -Wall -pedantic
OPTI            = -O3 -march=native
#OPTI = -pg # for profiling
#the following two lines are for compiling mpiflute with OpenMPI
MPICFLAGS	= -Wall -I/app/openmpi/include -I/app/openmpi/include/openmpi/ompi -pthread
MPILDFLAGS	= -L. -L/app/openmpi/lib -lmpi -lopen-rte -lopen-pal --export-dynamic -lm -lutil -lnsl -ldl -Wl
#uncomment the following two lines for compiling mpiflute with MPICH2 
#MPICFLAGS       = -Wall -I/opt/mpich2/include  -pthread
#MPILDFLAGS      = -L. -L/opt/mpich2/lib -lmpich -lrt --export-dynamic -lm -lutil -lnsl -ldl -Wl
LDFLAGS	= -lm
INCLUDES	= 
LIBS	= 
OBJS	= flute.o epimodel.o params.o epimodelparameters.o SFMT19937.o bnldev.o
MPIOBJS	= mpiflute.o mpiepimodel.o mpiparams.o mpiepimodelparameters.o mpiSFMT19937.o mpibnldev.o
DEFINES = -DVERBOSE
MPIDEFINES = -DPARALLEL

default: flute

flute: $(OBJS) Makefile
	$(CCLINKER) -o flute $(OBJS) $(OPTI) $(LDFLAGS) $(LIBS) $(DEFINES)

R0flute: $(OBJS) Makefile R0model.o R0model.h
	$(CCLINKER) -o R0flute R0model.o epimodel.o params.o epimodelparameters.o SFMT19937.o bnldev.o $(LDFLAGS) $(LIBS) $(DEFINES)

mpiflute: $(MPIOBJS) Makefile
	$(MPICCLINKER) -o mpiflute $(MPIOBJS) $(MPILDFLAGS) $(LIBS)

%.o: %.cpp epimodel.h epimodelparameters.h params.h Makefile
	$(CPP) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c $<

mpiflute.o: flute.cpp Makefile
	$(MPICPP) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -c flute.cpp -o mpiflute.o

mpiepimodel.o: epimodel.cpp epimodel.h Makefile
	$(MPICPP) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -c epimodel.cpp -o mpiepimodel.o

mpiparams.o: params.cpp params.h Makefile
	$(MPICPP) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -c params.cpp -o mpiparams.o

mpiepimodelparameters.o: epimodelparameters.cpp epimodelparameters.h epimodel.h Makefile
	$(MPICPP) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -c epimodelparameters.cpp -o mpiepimodelparameters.o

mpiSFMT19937.o: SFMT.c SFMT.h SFMT-params19937.h Makefile
	$(MPICC) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -std=c99 --param max-inline-insns-single=1800 -fno-strict-aliasing -Wmissing-prototypes -msse2 -DHAVE_SSE2 -DMEXP=19937 -DNDEBUG -c SFMT.c -o mpiSFMT19937.o

SFMT19937.o: SFMT.c SFMT.h SFMT-params19937.h Makefile
	$(CC) $(CFLAGS) $(OPTI)  -std=c99 --param max-inline-insns-single=1800 -fno-strict-aliasing -Wmissing-prototypes -msse2 -DHAVE_SSE2 -DMEXP=19937 -DNDEBUG $(INCLUDES) $(DEFINES) -c SFMT.c -o SFMT19937.o

bnldev.o: bnldev.c bnldev.h Makefile
	$(CC) $(CFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) -c bnldev.c -o bnldev.o

mpibnldev.o: bnldev.c bnldev.h Makefile
	$(MPICC) $(MPICFLAGS) $(OPTI) $(INCLUDES) $(DEFINES) $(MPIDEFINES) -c bnldev.c -o mpibnldev.o

zip: *.c *.cpp *.h Makefile README LICENSE HISTORY examples/*
	cd ..; zip flute/flute.zip flute/README flute/LICENSE flute/gpl.txt flute/HISTORY flute/Makefile flute/*.cpp flute/*.c flute/*.h flute/one-*dat flute/seattle-*dat flute/la-*dat flute/usa-*dat flute/examples/*

emacs:
	emacs Makefile *.R *.h *.c *.cpp README HISTORY &

clean:
	rm -f *.o flute mpiflute R0flute *~ output.mpiflute errormsg.mpiflute
	rm -f Summary? Tracts? Log? Individuals?
	echo 0 > run-number
