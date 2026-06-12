EXECUTABLE = chgr.x
OBJ_CHGR = \
    chgr_main.o chgr_horizontal_wavecut.o \
    sph_comp_index.o chgr_divor2wind.o chgr_vertical_interp.o \
    sph_fft2grid.o sph_coeff_lat.o \
    chgr_gaussian_interp.o chgr_sph_gaussian.o \
    chgr_sph_derivative.o chgr_read_sigma.o \
    chgr_write_sigma.o chgr_comp_rms.o chgr_make_header.o \
    chgr_div_coeff_lat.o chgr_get_sfcp.o \
    chgr_new_sigma.o chgr_sph_funct.o    \
    chgr_sph_operator.o \
    dyn_sigma_setup.o chgr_sph_transp.o \
    chgr_sph_transp_wind.o chgr_sum_coeff.o chgr_sum_coeff_wind.o \
    chgr_fourier_diff.o dyn_transpose_coeff.o chgr_cubic_interp.o \
    chgr_comp_maxmin.o chgr_linear_interp.o dfs_standard_temp.o \
    chgr_hybrid_coeff.o chgr_hybrid_setup.o
OBJ_DFS = \
	dfs_advection.o      dfs_divergence.o   dfs_wind_transform.o \
	dfs_fft_transform.o  dfs_global_mean.o  dfs_sincos_lat.o        dfs_cut_alias.o

OBJS =  ${OBJ_CHGR}

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
INSTALL_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXTRA_LIBS = 
EXTRA_LIBS_CHGR = 

LIBS = ${INSTALL_DIR}/common.a ${INSTALL_DIR}/share.a ${INSTALL_DIR}/sfcl.a ${LIB_DIR}/modelib.a ${LIB_DIR}/w3lib.a ${EXTRA_LIBS_CHGR}

EXP_PARENT = ..

${EXECUTABLE} : ${SRCS} ${OBJS} ${LIBS}
	if [ linux_intel = mac ] ; then \
		ranlib -s ${LIBS} ; \
	fi
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS} ${EXTRA_LIBS}
	${CP} $@ ${INSTALL_DIR}/chgr.x

${SRCS} :
	if [ -f ${EXP_PARENT}/mpi_dfs/$@ ]; then \
		ln -fs ${EXP_PARENT}/mpi_dfs/$@ $@; \
	elif [ -f ${EXP_PARENT}/fcst/$@ ]; then \
		ln -fs ${EXP_PARENT}/fcst/$@ $@; \
	fi

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
	- ${RM} -f *.f90 *.o  ${EXECUTABLE} \
		core ${INSTALL_DIR}/${EXECUTABLE}
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
