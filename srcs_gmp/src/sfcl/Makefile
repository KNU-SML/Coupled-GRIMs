LIB = sfcl.a

OBJ_SFCL = sfc_check_file.o   sfc_grib_file.o    sfc_read_file_driver.o    sfc_read_file.o   \
	sfc_cycle_driver.o    sfc_merge_file.o   sfc_interp_setup.o        sfc_interp_grid.o \
	sfc_quick_print.o     sfc_usgs_rough.o   sfc_usgs_albedo.o         sfc_land_type.o   \
	sfc_file_subst.o      sfc_tsfc_correc.o  dyn_lon_lat.o                               \
	sfc_grib_sum.o        sfc_get_latlon.o   sfc_max_unfwat.o

OBJ_SMP = sfc_hintp_smp.o     sfc_gausst62_smp.o 

OBJS=${OBJ_SFCL} ${OBJ_SMP}

SRC_SFCL = ${OBJ_SFCL:.o=.F90}
SRC_SMP  = ${OBJ_SMP:.o=.F90}
SRCS     = ${SRC_SFCL} ${SRC_SMP}

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
F77 = ifort
FORT_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
LOAD_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
INCLUDE_DIR = -I../.. -I../include 
LIB_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/libs_gmp/lib
INSTALLDLIB = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXP_PARENT = ..
LIBS = ${INSTALLDLIB}/share.a ${INSTALLDLIB}/common.a

${LIB} :  ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

${SRCS} :
	ln -fs ${EXP_PARENT}/fcst/$@ $@
  
depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	sed 's/module_trans/&_reduce/;s/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90

.F90.o :
	sed 's/module_trans/&_reduce/;s/ *use *[0-z_]*/&_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o *.a core 
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
# ----------------------------------------------------
