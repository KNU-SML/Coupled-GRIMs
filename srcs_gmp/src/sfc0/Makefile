EXECUTABLE = sfc0.x
OBJ_SFC0 = noa1toosu1.o noa1toosu2.o noa1tovic1.o noa1tovic2.o \
           osu1tonoa1.o osu1toosu2.o osu1tovic1.o osu1tovic2.o \
           osu2tonoa1.o osu2toosu1.o osu2tovic1.o osu2tovic2.o \
           vic1tonoa1.o vic1toosu1.o vic1toosu2.o vic1tovic2.o \
           vic2tonoa1.o vic2toosu1.o vic2toosu2.o vic2tovic1.o \
                 sfc0.o      frh2o.o     getslc.o \
             nodepara.o    sfcfcnv.o    vic1ini.o

OBJS = ${OBJ_SFC0}

SRCS = ${OBJS:.o=.F90}

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
INSTALLDLIB = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
INSTALL_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXTRA_LIBS = 
EXTRA_LIBS_SFC0 = 
EXP_PARENT = ..

LIBS = ${INSTALLDLIB}/sfcl.a   ${INSTALLDLIB}/common.a ${INSTALLDLIB}/share.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a   ${EXTRA_LIBS_SFC0}

${EXECUTABLE} : ${SRCS} ${OBJS} ${LIBS}
	if [ linux_intel = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}

${SRCS} :
	ln -fs ${EXP_PARENT}/fcst/$@ $@
 
depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90

.F90.o :
	sed 's/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i ; rm TEMP.F90
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o ${EXECUTABLE} core
	- find ./ -type l -exec rm -f {} \;

dist_clean : clean
	- ${RM} -rf rii_files Makefile
