LIB = mpi.a
LIB2 = ../share/share.a

OBJS = commpi.o dimset.o mpbcasti.o mpgetspd.o mpnk2nn.o mpnx2x.o mpsynrow.o \
equdiv.o mpbcastr.o mpbcastl.o mpgf2p.o mpnk2nx.o mpny2nl.o mpx2nx.o \
equlw12.o equdis.o mpdimset.o mpgp2f.o mpnl2ny.o mpsf2p.o spcshfli.o \
mpfine.o mpinit.o mpnn2n.o mpsp2f.o spcshflo.o preduceg.o \
mpabort.o mpgetken.o mplatall.o mpnn2nk.o mpsynall.o \
mpbcastc.o mpgetlat.o mpn2nn.o mpnx2nk.o mpsyncol.o \
mpgfpk2fk.o mpgf2xk.o mpxk2f.o mpgf2yk.o mpyk2f.o farg.o \
mpgp2f_box.o mpgf2p_box.o

SRCS = ${OBJS:.o=.F90}

MACHINE = linux_intel
SHELL = /bin/sh
AR = ar
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = icc
CPP = cpp -traditional-cpp
F77 = mpiifort
FORT_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
LOAD_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
INCLUDE_DIR = -I../.. -I../include 
LIB_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/libs_gmp/lib
INSTALLDLIB = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin

EXP_PARENT = ..

ALL = ${LIB}

all : ${ALL}

${SRCS} :
	ln -fs ${EXP_PARENT}/fcst/$@ $@

${LIB} : ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

yftest.x : yftest.o ${LIB} ${LIB2}
	${F77} ${LOAD_FLAGS} -o $@ yftest.o ${LIB} ${LIB2}

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	${CPP} -P ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i

.F90.o :
	${CPP} -P ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o core ${ALL}
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
