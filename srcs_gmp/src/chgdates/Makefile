EXECUTABLE = chgdates.x
OBJS = chgdates.o

MACHINE = linux_intel
SHELL = /bin/sh
CP = cp
RM = rm
CD = cd
CHMOD = chmod

CC = icc
CPP = cpp -traditional-cpp
F77 = ifort
FORT_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
LOAD_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
INCLUDE_DIR = -I../.. -I../include 
LIB_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/libs_gmp/lib
INSTALL_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXTRA_LIBS = 

LIBS = ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a ../share/share.a ../common/common.a

${EXECUTABLE} : ${OBJS} ${LIBS}
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS}
	${CP} $@ ${INSTALL_DIR}

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90

.F90.o :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o ${EXECUTABLE} core

dist_clean :
	- ${RM} -f *.f90 *.o *.x core Makefile
	- ${RM} -rf rii_files
