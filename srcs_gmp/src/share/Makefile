LIB = share.a

OBJ_SHARE = paramodel.o   comshare.o         \
	file_write_byte.o file_name.o        file_name_second.o    \
	sph_comp_index.o  share_funct_init.o share_funct_comp.o    \
	ftdpq.o           ftdpxg.o           ftheq.o               \
	fpvsq.o           ftlcl.o            ftmaq.o               \
	ftmaxg.o          isrcheqx.o         isrchfltx.o           \
	isrchnex.o        ismin.o            ismax.o               \
	isrcheq.o         isrchflt.o                               \
	snodfix.o         snosfc.o           znlavs.o              \
	znlody.o          znlsfc.o           znlwgt.o              \
	count.o           mpiset_sub.o       rmpiset_sub.o         \
	print_maxmin_six.o                   print_maxmin_seven.o  \
  coupling.o
#soojin_couple_above line

OBJ_DCMIP = dcmip_initial_conditions_test_1_2_3_v5.o           \
            dcmip_initial_conditions_test_4_v3.o               \
            dcmip_initial_conditions_test_5_v1.o               \
            simple_physics_v5_1.o
            

OBJS=${OBJ_SHARE} ${OBJ_DCMIP} paramodel_ser.o

SRC_SHARE = ${OBJ_SHARE:.o=.F90}
SRCS      = ${SRC_SHARE}

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
EXTC_FLAGS = 
LOAD_FLAGS = -r8 -O2 -DLINUX -convert big_endian -module ../include -m64 -mcmodel medium -I/usr/local/netcdf/4.6.1_intel21/include -L/usr/local/netcdf/4.6.1_intel21/lib -lnetcdff -lnetcdf
INCLUDE_DIR = -I../.. -I../include 
LIB_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/libs_gmp/lib
INSTALLDLIB = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin
EXP_PARENT = ..

ALL : constant.o ${LIB}

share.a : ${SRCS} ${OBJS}
	${AR} -r $@ ${OBJS}
	${CP} $@ ${INSTALLDLIB}
	- mv *.mod ../include

${SRCS} :
	ln -fs ${EXP_PARENT}/fcst/$@ $@

depend :
	${LIB_DIR}/../etc/mkdep

.SUFFIXES : .o .f90 .F90 .il

.F90.f90 :
	${CPP} -P ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i

.F90.o :
	${CPP} -P ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

.F90.il :
	${CPP} -P ${INCLUDE_DIR} $*.F90 >$*.i
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90 ; rm $*.i
	${F77} ${EXTC_FLAGS} $*.f90 -o $*.il
	
.f90.o :
	${F77} ${FORT_FLAGS} -c $*.f90

constant.o:constant.F90
	sed 's/paramodel/&_ser/;/define.h/{G;s/$$/#undef MP/;}' constant.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i ; rm -f  TEMP.F90
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90; rm -f $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

paramodel_ser.o:paramodel.F90
	sed 's/paramodel/&_ser/;/define.h/{G;s/$$/#undef MP/;}' paramodel.F90 > TEMP.F90
	${CPP} -P ${INCLUDE_DIR} TEMP.F90 >$*.i ; rm -f  TEMP.F90
	sed '/^ *$$/d;s/\/\*/\!/' $*.i >$*.f90; rm -f $*.i
	${F77} ${FORT_FLAGS} -c $*.f90

clean :
	- ${RM} -rf *.f90 *.o *.x *.a core *.il ../include/constant.mod ../include/paramodel*.mod
	- find ./ -type l -exec rm -f {} \;

dist_clean :
	make clean
	- ${RM} -rf rii_files Makefile
