#include <define.h>
march=_march_
sfcftyp=_sfcftyp_ 
lsoil=_lsoil_
lalbd=_lalbd_ 
mtnres=_mtnres_
nstype=_nstype_
nvtype=_nvtype_ 
mtnvar=_mtnvar_
#ifdef ROP
roms_levs=ROMS_LEVS
roms_startup=ROMS_STARTUP
roms_ini_theta_s=ROMS_INI_THETA_S
roms_ini_theta_b=ROMS_INI_THETA_B
roms_ini_tcline=ROMS_INI_TCLINE
roms_ini_tstart=ROMS_INI_TSTART
roms_ini_tend=ROMS_INI_TEND
roms_ini_smthinp=ROMS_INI_SMTHINP
roms_ini_nli=ROMS_INI_NLI
roms_bny_theta_s=ROMS_BNY_THETA_S
roms_bny_theta_b=ROMS_BNY_THETA_B
roms_bny_tcline=ROMS_BNY_TCLINE
roms_bny_tstart=ROMS_BNY_TSTART
roms_bny_tend=ROMS_BNY_TEND
roms_bny_smthinp=ROMS_BNY_SMTHINP
roms_for_nli=ROMS_FOR_NLI
roms_for_theta_s=ROMS_FOR_THETA_S
roms_for_theta_b=ROMS_FOR_THETA_B
roms_for_tcline=ROMS_FOR_TCLINE
roms_for_tstart=ROMS_FOR_TSTART
roms_for_tend=ROMS_FOR_TEND
roms_for_smthinp=ROMS_FOR_SMTHINP
roms_theta_s=ROMS_THETA_S
roms_theta_b=ROMS_THETA_B
roms_tcline=ROMS_TCLINE
roms_bny_nli=ROMS_BNY_NLI
#endif
#ifdef MP
mp=yes
#endif
#ifdef REDUCE_GRID
reduce_grid=yes
#endif

#ifdef R2_PHYSICS
r2=yes
#endif
#ifdef ICECHECK
sfc=yes
#endif
#ifdef SFCMRG
sfcmrg=yes
#endif
#ifdef OSULSM1
osulsm1=yes
#endif
#ifdef OSULSM2
osulsm2=yes
#endif
#ifdef NOALSM1
noalsm1=yes
#endif
#ifdef VICLSM1
viclsm1=yes
#endif
#ifdef USGS_SFC
usgs_sfc=yes
#endif
#ifdef STATSGO_SOIL
statsgo_soil=yes
#endif
#ifdef NCAR_EDIR
ncar_edir=yes
#endif
#ifdef KIMGWD
kimgwd=yes
#endif
#ifdef RAD_SMOOTH_CLOUD
rad_smooth_cloud=yes
#endif
#ifdef INSA_CLD
insa_cld=yes
#endif
#ifdef CLDADJ
cldadj=yes
#endif
#ifdef CLDSLINGO
cldslingo=yes
#endif
#ifdef ICE
ice=yes
#endif
#ifdef ICECLOUD
icecloud=yes
#endif
#ifdef CLD2
cld2=yes
#endif
#ifdef INTERACTIVE_STRATUS
interactive_stratus=yes
#endif
#ifdef LWRMDC
lwrmdc=yes
#endif
#ifdef SWRMDC
swrmdc=yes
#endif
#ifdef SW99UPDATE
sw99update=yes
#endif
#ifdef ALB4
alb4=yes
#endif
#ifdef SAS
sas=yes
#endif
#ifdef RAS
ras=yes
#endif
#ifdef RASV2
rasv2=yes
#endif
#ifdef CCMCNV
ccmcnv=yes
#endif
#ifdef CCMSNV
ccmsnv=yes
#endif
#ifdef DIFFQ
diffq=yes
#endif
#ifdef RDAMPUX
rdampux=yes
#endif
#ifdef HYDRO
hydro=yes
#endif
#ifdef NOAHYDRO
noahydro=yes
#endif
#ifdef PSPLIT
psplit=yes
#endif
#ifdef STDAMP
stdamp=yes
#endif
#ifdef STDAMP3
stdamp3=yes
#endif
#ifdef STDAMP_P
stdamp_p=yes
#endif
#ifdef DG
dg=yes
#endif
#ifdef MRG_POST
mrg_post=yes
#endif
#ifdef PGB_PARALLEL
pgb_parallel=yes
#endif
#ifdef PGB_SEMIPARALLEL
pgb_semiparallel=yes
#endif
#ifdef PGB_NOPARALLEL
pgb_noparallel=yes
#endif
#ifdef PGBGATHER
pgbgather=yes
#endif
#ifdef PGBMPI
pgbmpi=yes
#endif
#ifdef VERY_LARGE_DOMAIN
very_large_domain=yes
#endif
#ifdef DG3
dg3=yes
#endif
#ifdef CLR
clr=yes
#endif
#ifdef DGP
ken=yes
#endif
#ifdef DBG
dbg=yes
#endif
#ifdef LFM
lfm=yes
#endif
#ifdef LFC
lfc=yes
#endif
#ifdef DFI
dfi=yes
#endif
#ifdef ONELOOP
oneloop=yes
#endif
#ifdef SKIPSFCMRG
skipsfcmrg=yes
#endif
#ifdef NOPRINT
noprint=yes
#endif
#ifdef CHGR_SMPL
chgr_smpl=yes
#endif
#ifdef SSI
ssi=yes
#endif
#ifdef CRAY_THREAD
cray_thread=yes
#endif
#ifdef ORIGIN_THREAD
origin_thread=yes
#endif
#ifdef OPENMP
openmp=yes
#endif
#ifdef REAL4_W3LIB
real4_w3lib=yes
#endif
#ifdef DYNAMIC_ALLOC
dynamic_alloc=yes
#endif
#ifdef MINV
minv=yes
#endif
#ifdef RFFTMLT
rfftmlt=yes
#endif
#ifdef ASSIGN
assign=yes
#endif
#ifdef GETENV
getenv=yes
#endif
#ifdef PXFGETENV
pxfgetenv=yes
#endif
#ifdef SGEMVX1
sgemvx1=yes
#endif
#ifdef SGERX1
sgerx1=yes
#endif
#ifdef FL2I
fl2i=yes
#endif
#ifdef DCRFT
dcrft=yes
#endif
#ifdef RANF
ranf=yes
#endif
#ifdef CRAY_BUFRLIB
cray_bufrlib=yes
#endif
#ifdef GTOPO30
gtopo30=yes
#endif
#ifdef SETMTNZERO
setmtnzero=yes
#endif
#ifdef RMPVECTORIZE
rsmvectorize=yes
#endif
