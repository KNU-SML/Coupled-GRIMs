EXECUTABLE = co2.x
OBJS = co2_temp.o co2o3_step3.o co2o3_step4.o co2o3_step2.o \
co2_ozone_main.o co2o3_step1.o co2_tau_interp.o co2_mixing_ratio.o \
co2_new_sigma.o co2_ozone_data.o path.o co2_quard_interp.o \
co2_quardratic.o co2_trans_funct.o co2_sigma_interp.o co2_interp.o

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

LIBS = ${INSTALL_DIR}/share.a ${INSTALL_DIR}/common.a ${INSTALL_DIR}/sfcl.a

${EXECUTABLE} : ${OBJS} ${LIBS}
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS}
	${CP} $@ ${INSTALL_DIR}

${INSTALL_DIR}/.co2 : ${EXECUTABLE}
	sh ./co2.sh || exit 8

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TMP.F90
	${CPP} -P ${INCLUDE_DIR} TMP.F90 >$*.i ; rm TMP.F90
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i

.F90.o :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TMP.F90
	${CPP} -P ${INCLUDE_DIR} TMP.F90 >$*.i ; rm TMP.F90
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o *.x core ${INSTALL_DIR}/.co2

dist_clean :
	- ${RM} -f *.f90 *.o *.x core Makefile ${INSTALL_DIR}/.co2
	- ${RM} -rf rii_files
