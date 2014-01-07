	PGM=HydroC
HEADER=$(wildcard *.h)
SRC=cmpflx.c conservar.c equation_of_state.c hydro_godunov.c main.c parametres.c riemann.c trace.c vtkfile.c compute_deltat.c constoprim.c hydro_funcs.c hydro_utils.c make_boundary.c qleftright.c slope.c utils.c SplitSurface.c

BCK_SRC=hydro_godunov.c riemann.c hggodunov_cuda.cu
HOST=$(shell uname -n | sed 's/[0-9]//g')

OBJ = $(patsubst %.c, %.o, ${SRC})
DEP = $(patsubst %.c, %.d, ${SRC}) $(DEP2)

OPT=-O3 -DNDEBUG -DFAST -std=c99 -fopenmp
CC=mpicc

CFLAGS=$(OPT) -Wall -Wno-unknown-pragmas

$(PGM): $(DEP) $(OBJ)
	$(CC) -o $(PGM) $(OPT) $(OBJ) -lm $(HMPPENDFLAGS)


clean   :
	-/bin/rm -f *.d *.o *.so *~ *.vts  *.bic *.bak *.cu ${PGM} 
	-/bin/rm -f *.hdpp.* *.preproc.* *.translated.* *.inline.c *.extracted.c
	-/bin/rm -f *.lst *.ptx *.cub *.fatbin
	-/bin/rm -f *.linkinfo
	-/bin/rm -rf Dep/ *.pvd


-include $(DEP)

.SUFFIXES:  .o .cu .d .c

.c.d:
	$(CC) ${CPPFLAGS} ${CFLAGS} $(HFLAGS) -M $< -o $@

.c.o    :
	${CC} ${CFLAGS} -c $< -o $@ $(HMPPENDFLAGS)
.cpp.o  :
	${CC} ${CFLAGS} -c $<
.cu.o:
	${NVCC} ${CFLAGS} -c $<


FORCE:
#EOF
