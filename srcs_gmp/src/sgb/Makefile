EXECUTABLE = sgb.x
OBJ_SGB =  \
	post_sgb_main.o             post_sgb_solver.o           \
	post_byte2char.o            post_wind_transform.o       \
	post_gauss2latlon.o         post_interp_spline.o        \
	post_pgb_3d.o               post_grad_scalar.o          \
	post_spectral_field.o       post_quick_print.o          \
	post_grib_binary.o          post_sph_lat.o              \
	post_pole_value.o           post_sph_fft_coeff.o        \
	post_legendre_poly.o        post_synth_coeff.o          \
	post_read.o                 post_interp_sigma.o         \
	post_trans_output_grid.o    post_trans_wave2grid.o      \
	post_dfs_lat.o

OBJ_DFS1 =  \
	dfs_fft_transform.o         dfs_advection.o             \
	dfs_wind_transform.o        dfs_global_mean.o           \
	dfs_divergence.o            dfs_semi_implicit.o         \
	dfs_cut_alias.o

OBJ_DFS2 = post_dfs_sfcp.o          post_dfs_latlon.o           post_dfs2sph_trans.o

OBJ_DFS=${OBJ_DFS1} ${OBJ_DFS2}

OBJ_ETC=file_make_grib.o

OBJS =  ${OBJ_SGB} ${OBJ_ETC}

SRC_DFS = ${OBJ_DFS:.o=.F90}
SRC_SGB = ${OBJ_SGB:.o=.F90}
SRC_ETC = ${OBJ_ETC:.o=.F90}
SRCS =  ${SRC_SGB} ${SRC_ETC}

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
INSTALLDLIB = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXTRA_LIBS = 
EXP_PARENT = ..

LIBS = ${INSTALLDLIB}/share.a ${INSTALLDLIB}/common.a ${LIB_DIR}/w3lib.a ${LIB_DIR}/modelib.a

${EXECUTABLE} : ${SRCS} ${OBJS} ${LIBS}
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} ${LIBS}
	${CP} $@ ${INSTALL_DIR}

${SRCS} :
	if [ -f ../pgb/$@ ]; then    \
		ln -fs ${EXP_PARENT}/pgb/$@ $@; \
	fi

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90

.F90.f90 :
	sed 's/ *module *[0-z_]*/&_ser/;s/ *use *[0-z_]*/&_ser/;s/_ser_ser/_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90

.F90.o :
	sed 's/ *module *[0-z_]*/&_ser/;s/ *use *[0-z_]*/&_ser/;s/_ser_ser/_ser/;s/use constant_ser/use constant/g;/define.h/{G;s/$$/#undef MP/;}' $*.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i TEMP.F90
	${F77} ${FORT_FLAGS} -c $*.f90

.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -f *.f90 *.o  ${EXECUTABLE} core
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
