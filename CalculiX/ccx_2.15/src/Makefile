
CFLAGS = -Wall -O3 -pthread -I ../../../SPOOLES.2.2 -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT -DUSE_MT
FFLAGS = -Wall -O3 

CC=gcc-8
FC=gfortran-8
AR=ar

.c.o :
	$(CC) $(CFLAGS) -c $<
.f.o :
	$(FC) $(FFLAGS) -c $<

include Makefile.inc

SCCXMAIN = ccx_2.15.c

OCCXF = $(SCCXF:.f=.o)
OCCXC = $(SCCXC:.c=.o)
OCCXMAIN = $(SCCXMAIN:.c=.o)

DIR=../../../SPOOLES.2.2

LIBS = \
       $(DIR)/spooles.a \
       ../../../ARPACK/libarpack_MAC.a

ccx_2.15: $(OCCXMAIN) ccx_2.15.a  $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c ccx_2.15.c; $(FC)  -Wall -O3 -o $@ $(OCCXMAIN) ccx_2.15.a $(LIBS)

ccx_2.15.a: $(OCCXF) $(OCCXC)
	$(AR) vr $@ $?
