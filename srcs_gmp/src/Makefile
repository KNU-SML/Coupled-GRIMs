#
SHELL = /bin/sh
RM = rm
CD = cd
MAKE = make
INSTALL_DIR = /rhome/soo/NEMO/couple/couple_ice/gmp_v3.2/srcs_gmp/bin

PROGS = 			share 			common 			common_par 			mpi 			 		sfcl 		sfcl_par 		chgr 		chgdates 		co2 		mtn 		pgb 		fcst 		sgb 		fgb   		sfc0   		p2sig 	

SSI_PROGS = 
GOP_PROGS = 
GOM_PROGS = 
ROP_PROGS = 
GSI_PROGS = 

all : 
	for prog in ${PROGS} ${GOP_PROGS} ${GOM_PROGS} ${ROP_PROGS} ${GSI_PROGS}; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

ssi :
	for prog in ${SSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gop :
	for prog in ${GOP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gom :
	for prog in ${GOM_PROGS} ; do \
		if [ -s gom/$$prog ] ; then \
			${CD} gom/$$prog ; \
			echo now in gom/$$prog ; \
			${MAKE} || exit 8 ; \
			${CD} ../.. ; \
		fi ; \
	done

rop :
	for prog in ${ROP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

gsi :
	for prog in ${GSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} || exit 8 ; \
			${CD} .. ; \
		fi ; \
	done

clean :  
	for prog in ${PROGS} ${GOP_PROGS} ${ROP_PROGS} ${GSI_PROGS}; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done
	- ${RM} -f include/*.mod
	- ${RM} -f include/*.MOD

clean_ssi :  
	for prog in ${SSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

clean_rop :  
	for prog in ${ROP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

clean_gsi :  
	for prog in ${GSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done
	- ${RM} -f gsi_include/*.mod
	- ${RM} -f gsi_include/*.MOD

clean_gop :  
	for prog in ${GOP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} clean ; \
			${CD} .. ; \
		fi ; \
	done

dist_clean :  
	for prog in ${PROGS} ${GOP_PROGS} ${ROP_PROGS} ${GSI_PROGS}; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} dist_clean ; \
			${CD} .. ; \
		fi ; \
	done
	- ${RM} -rf ${INSTALL_DIR} Makefile
	- ${RM} -f include/*.mod
	- ${RM} -f include/*.MOD

dist_clean_ssi :  
	for prog in ${SSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} dist_clean ; \
			${CD} .. ; \
		fi ; \
	done
dist_clean_gop :  
	for prog in ${GOP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} dist_clean ; \
			${CD} .. ; \
		fi ; \
	done
dist_clean_rop :  
	for prog in ${ROP_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} dist_clean ; \
			${CD} .. ; \
		fi ; \
	done
dist_clean_gsi :  
	for prog in ${GSI_PROGS} ; do \
		if [ -s $$prog ] ; then \
			${CD} $$prog ; \
			echo now in $$prog ; \
			${MAKE} dist_clean ; \
			${CD} .. ; \
		fi ; \
	done
