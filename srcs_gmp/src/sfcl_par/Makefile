LIB = sfcl_par.a
OBJ_SFCL = sfc_check_file.o   sfc_grib_file.o    sfc_read_file_driver.o    sfc_read_file.o   \
	sfc_cycle_driver.o    sfc_merge_file.o   sfc_interp_setup.o        sfc_interp_grid.o \
	sfc_quick_print.o     sfc_usgs_rough.o   sfc_usgs_albedo.o         sfc_land_type.o   \
	sfc_file_subst.o      sfc_tsfc_correc.o  dyn_lon_lat.o                               \
	sfc_grib_sum.o        sfc_get_latlon.o   sfc_max_unfwat.o

OBJS = ${OBJ_SFCL}

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

${LIB} : ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}

${SRCS} :
	if [ -f ../sfcl/$@ ]; then \
		ln -fs ../sfcl/$@ $@; \
	elif [ -f ../fcst/$@ ]; then \
		ln -fs ../fcst/$@ $@; \
	fi

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
	- ${RM} -f *.F90 *.f90 *.o *.a core ${SRCS}
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
