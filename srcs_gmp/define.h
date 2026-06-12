#include <machine.h>
/* ----- grims-gmp global spherical harmonics (sph) setup -----  */

#undef RMP          /* rmp option off*/
#undef GDAS         /* GDAS option off */

/* dynamics options */
#define GRIMS        /* grims or wrf */
#define MPFINE        /* mpi finalize */
#define NIM_LAHEY
#define NIM_BITWISE
#define REDUCE_GRID  /* reduced grid for sph  for _jcap_ > 42 */
#define DIFFQ        /* horizontal diffusion of moisture */
#define SPH_GFS_DIFFUSION /* enhanced horizontal diffusion with height */
#define HYBRID       /* vertical cordiate with hybrid sigma-pressure */
#define MTN4MIN      /* use global 4 min mountain instead of 8 min */
#undef REIGH_FRICTION /* enhanced diffusion in the upper stratosphere */
#undef  NISLQ        /* semi-lagrangian advection for gases and water */
#undef STOCH         /* stochastic tendency perturbation in dyn and phy */
#undef GTOPO30       /* use original gtopo 30 data for mountain in rmp */
#undef SETMTNZERO    /* set orography to zero over ocean */
#undef PSPLIT        /* process splitting time scheme for both gmp and rmp */
#undef LFM           /* low-frequency model setup */
#undef LFC           /* lsm forcing for smp */
#undef DFI           /* digital filter initilization */
#undef ONELOOP       /* single loop for dynamics and physics */
#undef SKIPSFCMRG    /* surface merge logical flag */
#undef CHGR_SMPL     /* change resolution option for gmp */

/*soojin_couple*/
#define OASIS3_mct   /* OASIS3 mct for coupling */
#define AOMG         /* Atmosphere-Ocean coupled Model Global version */

#undef NONHYD        /* nonhydrostatic opton */
#define CORIOLIS     /* Coriolis terms */
#define RDAMPUX      /* rmp speed dependent diffusion to allow longer timestep */
#define STDAMP       /* spectral tendency damping scheme for rsm */
#define STDAMP3      /* area average log ps correction to be used with stdamp */
#define STDAMP_VER   /* vertical weighted damping coefficient with stdamp */
#undef STDAMP_P      /* spectral damping of perturbation tendency */
#define LTB_RELAX    /* lateral boundary relaxation */
#define TOP_RELAX    /* top-layer boundary relaxation */
#undef VDIFF         /* gridpoint second-order vertical diffusion */
#define _kvdif_ 75.  /* vertical diffusion coefficient (m2/s) */
#undef RMPVECTORIZE  /* vectorize mpi rsm physics */
#undef RCHGR_SMPL    /* change resolution option for rmp */
#undef VERY_LARGE_DOMAIN /* rmp for large arrays model set up */
#undef GMPDAMP       /* gmp damping with spectral nudging for wave > 2000 km */
#define _gdamp_crlat_ 90 /* gmp damping critical latitude  < 90 */

#undef DBG           /* debugging option */
#undef NOPRINT       /* suppress print out */
#undef MRG_POST      /* merge postprocessor into forecast code */
#undef DG3           /* 3d diagnostic grib file */
#undef DGP           /* diagnostic point binary output */

/* physics options */
#undef R2_PHYSICS    /* simulate reanalysis-2 physics */
#define  LWRMDC      /* ming-dah chou long wave radiation */
#undef  SWRMDC       /* ming-dah chou short wave radiation */
#undef  RRTMGLW      /* RRTMG long wave radiation */      
#undef  RRTMGSW      /* RRTMG short wave radiation */  
#undef  CLDADJ       /* empirical cloudiness adjustment */
#define vvadj        /* vertical velocity adjustment for low clouds */

#define  Z0T         /* thermal roughness length */
#define  VSGD        /* subgrid surface wind velocity scale */
#define SFCMRG       /* surface merge flag of guess and climatology */
#undef OSULSM1       /* osu lsm with homogeneous   vegatation and soil type */
#undef OSULSM2       /* osu lsm with heterogeneous vegatation and soil type */
#define  NOALSM1     /* noah land scheme */
#define  KIMHONG2010 /* kim and hong ocean surface layer model */
#undef  HYDRO        /* residual free hydrology with osu2 */
#undef  NOBSFLW      /* remove base flow from osu2 */
#undef  NOAHYDRO     /* residual free hydrology with noah */
#undef  USGS_SFC     /* usgs surface and vegetation data */
#ifdef RMP
#define  USGS_SFC
#endif
#undef VIC
#undef  VICLSM1      /* vic single tile scheme */

/*soojin_couple start*/
#undef OMLWRF       /* ocean mixed layer coupling */
/*soojin_couple end*/
#define RIVER        /* river flow model : half deg resolution */

#undef MRFPBL        /* ysu pbl scheme or mrf scheme */
#define YSUPBL       /* ysu pbl scheme or mrf scheme */
#define MAR2012_YSU  /* revised pr and mixing length scale as in wrf v3.4 */
#define JUN2012_YSU  /* bug fix in ws in stable pbl as in wrf 3.4.1 */
#undef YSUTKE        /* 3d tke for ysu pbl */

#define _ngases_ 1   /* number of gase 1:o3 2:tke */
#undef O3CHEM        /* kim jhoon's ozone chemistry */

#define KSAS
#undef SAS          /* simplified arakawa-chubert cps */
#undef UV_CMT        /* convective momentum transport in sas */
#undef NCEP2010      /* han and pan sas changes in july 2010 */
#undef DX_FACTOR_NSAS /* critical omega is optimized with respect to dx */
#undef SAS_CCN       /* ccn number concentration in sas */
#undef  RAS          /* relaxed arakawa-schubert ras cps */
#undef  RASV2        /* ras with downdraft */
#undef  CCMCNV       /* zhang-mcfarlane ccm cps */
#undef  KF2          /* kain-fritsch kf2 cps */
#undef  KUO          /* kuo cps */
#undef  CPS_ENS      /* cps ensemble */

#ifdef KSAS
#undef SAS
#define CPS_KSAS_DIUR   /* diurnal cycle with cldwrk in pbl in sas cps */
#define EXPLICIT_CLOUDINESS
#define CPS_KSAS_SIGMA  /* with sigma value */
#define DET_HYDRO
#endif

#undef  MRFSCV       /* tiedke shallow convection */
#define GRIMSCV      /* explicitly coupled shallow convection and pbl */
#define FEB2013_GRIMSCV /* revisions in cloud top and kz profile */
#undef  CCMSCV       /* hack shallow convection */
#undef  HANPANSCV    /* han and pan shallow convection */

/*soojin_couple*/
#ifdef  AOMG
#undef WSM1
#undef WSM2
#define WSM3
#undef WSM5
#undef WSM6
#undef WDM5
#undef WDM6
#else /*AOMG*/
#define WSM1         /* prognostic qv mps */
#undef  WSM2         /* prognostic qv qci mps */
#undef  WSM3         /* prognostic qv qci qrs mps */
#undef  WSM5         /* prognostic qv qc qr qi qs mps */
#undef  WSM6         /* prognostic qv qc qr qi qs qg mps */
#undef  WDM5         /* prognostic qv qc qr qi qs nccn nc nr mps */
#undef  WDM6         /* prognostic qv qc qr qi qs qg nccn nc nr mps */
#endif /*AOMG*/

#undef   NOGWD       /* no gravity wave drag */
#define  KAGWD       /* kim-arakawa orography gravity wave drag */
#define  GWDC	     /* Chun Baik convection gravity wave drag */

#undef  DRY_SOIL     /* initial soil moisture at wilting point */
#undef  WET_SOIL     /* initial soil moisture at saturation point */
#undef  ANL_SOIL     /* soil moisture from a analysis */
#undef  CLM_SOIL     /* soil moisture from a climatology */
#undef  NO_DIURNAL   /* no diurnal cycle for radiation */
#define _march_ mpi
/* dependencies */

#ifdef GRIMS
#undef WRF
#endif

#ifdef DFS
#define REIGH_FRICTION
#undef REDUCE_GRID
#endif

/* soojin_couple*/
#ifdef AOMG
#undef REDUCE_GRID
#define NCEP2010
#endif

#ifdef RMP
#undef MRG_POST
#define VCI
#undef CLDADJ
#define LWRMDC
#undef SFCMRG
#undef REDUCE_GRID
#define HD6ORDER
#undef HYBRID
#ifdef NONHYD
#define NONHYD_HYD
#endif
#ifdef IDEAL
#define NO_PHYSICS
#define NOGWD
#undef  CORIOLIS      /* Coriolis terms */
#undef  RDAMPUX
#undef  STDAMP
#undef  STDAMP3
#undef  STDAMP_VER
#undef  STDAMP_P
#if defined(COLD) || defined(WARM)
#undef  LTB_RELAX     /* lateral boundary relaxation */
#undef  TOP_RELAX     /* top-layer boundary relaxation */
#define POST_SEC      /* time interval in seconds to write output file */
#endif
#ifdef COLD
#define VDIFF         /* gridpoint second-order vertical diffusion */
#define _kvdif_ 75.   /* vertical diffusion coefficient (m2/s) */
#endif
#else
#undef COLD
#undef WARM
#undef HILL
#endif
#endif

#ifdef SMP
#undef MRG_POST
#undef REDUCE_GRID
#undef HYBRID
#define SKIPSFCMRG
#undef  DG3
#undef  DG
#ifndef OSULSM1
#define SMP_RA2SFC
#endif
#undef PSPLIT
#undef SFCMRG
#undef  R_CNVCLD
#undef  NEW_AUTO
#undef RIVER
#undef GWDC
#endif

#ifdef CHEM
#undef  REDUCE_GRID
#undef  LWRMDC        /* ming-dah chou long wave radiation */
#undef  SWRMDC        /* ming-dah chou short wave radiation */
#define RRTMGLW       /* RRTMG long wave radiation */      
#define RRTMGSW       /* RRTMG short wave radiation */  
#define NISLQ         /* Chemistry module needs semi-Lagrangian */
#define NCEP2010      /* Chemistry module needs NCEP2010 */
#undef  DST1BIN       /* Simplified dust simulation with 1 size bin */
#undef  VOLCOFF       /* Turn off daily volcano SO2 emission */
#endif

#ifdef NIM
#define NONHYD      /* NIM is NONHYD */
#undef NONHYD_HYD
#undef GRIMS          /* rmp option off */
#undef RMP            /* rmp option off */
#undef GDAS           /* GDAS option off */
#undef DFS            /* dfs */
#undef REDUCE_GRID
#undef RIVER
#ifdef GWDC
#define AUG2012_GWDC
#endif
#endif /* NIM end */

#ifdef KIM
#undef _ngases_
#define _ngases_ 0   /* should be 0 for prescribed o3 */
#undef NIM_LAHEY
#undef NIM_BITWISE
#endif

#ifdef NIMAQUA
#define AQUA_PLANET
#undef  O3SYM         /* ozone symmetric between north and south hem */
#undef  O3SPR         /* fixed ozone as spring distribution in nor hem */
#define YSUPBL
#define NCEP2010
#define Z0T
#define VSGD
#endif

#ifdef RIVER
#define _io2_ 720    /* i-grid number for river model */
#define _jo2_ 361    /* j-grid number for river model */
#define _delxo2_ 0.5 /* river model resolution in degree */
#endif

#ifdef OSULSM1
#define _lsoil_ 2
#define _nsoil_ 0
#define _msub_ 0
#define _sfcftyp_ osu1
#define _lalbd_ 1
#endif

#ifdef OSULSM2
#define _lsoil_ 2
#define _nsoil_ 0
#define _msub_ 0
#define _sfcftyp_ osu2
#define _lalbd_ 4
#endif

#ifdef NOALSM1
#define _lsoil_ 4
#define _nsoil_ 0
#define _msub_ 0
#define _sfcftyp_ noa1
#define _lalbd_ 4
#endif

#ifdef VICLSM1
#define _lsoil_ 3
#define _nsoil_ 5
#define _msub_ 2
#define _sfcftyp_ vic1
#define _lalbd_ 4
#endif

#ifdef VICLSM2
#define _lsoil_ 3
#define _nsoil_ 5
#define _msub_ 12
#define _sfcftyp_ vic1
#define _lalbd_ 4
#endif

#ifndef OSULSM1
#undef R2_PHYSICS
#endif

#ifdef R2_PHYSICS
#undef RAD_SMOOTH_CLOUD
#undef CLDADJ
#undef LWRMDC
#undef SWRMDC
#define SAS
#undef  RAS
#undef RASV2
#undef CCMCNV
#undef CCMSCV
#define SFC
#define SFCMRG
#endif

#ifdef MRG_POST
#ifndef MP
#undef MPIGRIB
#else
#ifdef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#undef PGB_NOPARALLEL
#undef MPIGRIB
#undef PGBGATHER
#endif
#ifdef PGB_SEMIPARALLEL
#undef PGB_PARALLEL
#undef PGB_NOPARALLEL
#undef MPIGRIB
#define PGBGATHER
#endif
#ifdef PGB_NOPARALLEL
#undef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#define MPIGRIB
#undef PGBGATHER
#endif
#endif
#else
#undef PGB_NOPARALLEL
#undef PGB_PARALLEL
#undef PGB_SEMIPARALLEL
#undef MPIGRIB
#undef PGBGATHER
#endif

#ifdef NOALSM1
#undef HYDRO
#endif

#ifdef OSULSM1
#undef OSULSM2
#undef NOALSM1
#undef VICLSM1
#undef MRGLSM
#undef USGS_SFC
#undef NOAHYDRO
#endif

#ifdef OSULSM2
#undef OSULSM1
#undef NOALSM1
#undef VICLSM1
#undef MRGLSM
#undef NOAHYDRO
#endif

#ifdef NOALSM1
#undef OSULSM2
#define MRGLSM
#undef USGS_SFC
#undef OSULSM1
#undef VICLSM1
#endif

#undef OSU
#ifdef OSULSM1
#define OSU
#endif
#ifdef OSULSM2
#define OSU
#endif

#undef NOA
#ifdef NOALSM1
#define NOA
#endif
#ifdef NOALSM2
#define NOA
#endif

#undef VIC
#ifdef VICLSM1
#define VIC
#endif
#ifdef VICLSM2
#define VIC
#endif

#ifdef VICLSM1
#define MRGLSM
#undef OSULSM1
#undef OSULSM2
#undef NOALSM1
#undef USGS_SFC
#undef NOAHYDRO
#endif

#ifdef USGS_SFC
#define _nstype_ 16
#define _nvtype_ 12
#define USGS
#define OSULSM2
#define STATSGO_SOIL
#define NCAR_EDIR
#define SM_SPINUP
#else
#define _nstype_ 9
#define _nvtype_ 13
#undef USGS
#undef NCAR_EDIR
#undef STATSGO_SOIL
#undef SM_SPINUP
#endif

#ifdef WSM1
#define _nwmass_ 1
#endif
#ifdef WSM2
#define _nwmass_ 2
#endif
#ifdef WSM3
#define _nwmass_ 3
#endif
#if defined (WSM5) || defined (WDM5)
#define _nwmass_ 5
#endif
#if defined (WSM6) || defined (WDM6)
#define _nwmass_ 6
#endif

#if defined (WDM5) || defined (WDM6)
#define _nwsize_ 3
#else
#define _nwsize_ 0
#endif

#if (_nwmass_ == 1)
#undef ICE
#undef ICECLOUD
#else
#define ICE
#define ICECLOUD
#ifdef NCEP2010
#define DET_HYDRO
#endif
#endif

#define INTERACTIVE_STRATUS
#if (_nwmass_ == 1)
#undef INTERACTIVE_STRATUS
#endif

#ifndef RMP
#if (_nwmass_ > 1) || (_ngases_ > 1)
#define NISLQ
#endif
#endif

#ifdef NOGWD
#undef KAGWD
#undef GWDC
#else
#ifdef GTOPO30
#define KAGWD
#endif
#endif

#ifdef FLOW_BLOCKING
#define _mtnvar_ 11
#define _rmtnvar_ 11
#else
#ifdef KAGWD
#define _mtnvar_ 10
#define _rmtnvar_ 10
#else
#define _mtnvar_ 1
#define _rmtnvar_ 1
#endif
#endif

#ifndef RMP
#undef RMPVECTORIZE
#else
#ifndef MP
#undef RMPVECTORIZE
#endif
#endif

#if ( _igrd_ > 256 )
#define VERY_LARGE_DOMAIN
#endif

#ifdef VERY_LARGE_DOMAIN
#define RINPG_SMPL
#endif

#ifdef DG3
#define DG
#define CLR
#endif

#ifdef DGP
#undef CLR
#endif

#ifdef SAS
#if defined(SAS2005)
#undef  SAS_DIAG
#define  NEW_TRIGGER
#define  UV_CMT
#define  MUL_CLDTOP
#define _ncldtop_ 1
#define  CLM_CWF
#define NEW_AUTO
#elif  defined (NCEP2010)
#undef  SAS_DIAG
#undef  NEW_TRIGGER
#define  UV_CMT
#undef  MUL_CLDTOP
#define _ncldtop_ 1
#undef  CLM_CWF
#undef NEW_AUTO
#else
#undef  SAS_DIAG
#undef  NEW_TRIGGER
/* #undef  UV_CMT */
#undef  MUL_CLDTOP
#define _ncldtop_ 1
#undef  CLM_CWF
#undef  NEW_AUTO
#endif
#ifdef SAS_CCN
#define _iccn_ 100  /* initial CCN number concentration. 100, 1000, or 7000 */
#endif
#endif
!
#ifdef RAS
#undef RASV2
#undef SAS
#undef CCMCNV
#undef CCMSCV
#undef KUO
#endif
#ifdef RASV2
#undef RAS
#undef SAS
#undef CCMCNV
#undef CCMSCV
#undef KUO
#endif
#ifdef KF2
#undef RAS
#undef RASV2
#undef SAS
#undef CCMCNV
#undef CCMSCV
#undef KUO
#endif
#ifdef CPS_ENS
#define SAS
#define RAS
#undef  RASV2
#define CCMCNV
#undef  CCMSCV
#define KF2
#undef  KUO
#endif
#if defined (DRY_SOIL) || defined(WET_SOIL) || defined(ANL_SOIL) || defined(CLM_SOIL)
#define SM_SPINUP
#endif

#if defined (WDM5) || defined (WDM6)
#define WDM
#define _iccn_ 100000000  /* initial CCN number concentration. 1.e8, 10.e8, or 70.e8 */
#endif

#if defined SAS || defined KSAS
#define GWDC
#else 
#undef GWDC
#endif
 
#ifdef GTOPO30
#define _mtnres_ 0.5
#define _rmtnres_ 0.5
#else
#ifdef MTN4MIN
#define _mtnres_ 4
#define _rmtnres_ 4
#else
#define _mtnres_ 8
#define _rmtnres_ 4
#endif
#endif

#ifdef NISLQ
#undef NISLQ_MASS     /* semi-Lagrangian with mass conservation   */
#define NISLQ_MONO    /* semi-Lagrangian with only monotonicity   */
#define NISLQ_PHYS    /* YSU: upward borrowing, SAS: scattering   */
#define NISLQ_GRIB     /* write values by nislq */
#undef NISLQ_CHECK    /* check negative values in physics */
#undef NISLQ_MONOMASS /* mass adjustment for global conservation  */
#undef TRCON          /* check global conservation */
#endif

#ifdef STOCH
#define STOCH_TIME    /* time dependency of interval */
#define STOCH_LEVEL   /* layer dependency of interval */
#define STOCH_INT     /* temporal correlation (sampling time interval) */
#undef STOCH_SC      /* spatial correlation */
#undef STOCH_PRINT    /* print out */
#undef STOCH_DBG      /* for debug */

#define DYN_U         /* for dynamics */
#define DYN_V
#define DYN_T
#define DYN_Q
#define PHY_U         /* for physics */
#define PHY_V
#define PHY_T
#define PHY_Q
#undef NLT_U          /* for nonlinear tendency in dynamcis */
#undef NLT_V
#undef NLT_T
#undef NLT_Q
#endif

#ifdef AQUA_PLANET
#define SKIPSFCMRG
#undef REDUCE_GRID
#undef OMLWRF
#endif

#ifdef DCMIP
#undef _ngases_
#define _ngases_ 6
#define NO_PHYSICS
#undef REDUCE_GRID
#undef RIVER
#ifdef NISLQ
#define NISLQ_WRITE
#define NCEP_GFS_DIFFUSION
#endif
#endif /* DCMIP end */
#include <define_mpi.h>
#undef CRAY_THREAD
#undef ORIGIN_THREAD
#undef OPENMP
#undef REAL4_W3LIB
#undef DYNAMIC_ALLOC
#undef MINV
#undef RFFTMLT
#undef ASSIGN
#undef GETENV
#undef PXFGETENV
#undef SGEMVX1
#undef SGERX1
#undef FL2I
#undef DCRFT
#undef RANF
#undef CRAY_BUFRLIB
#undef ASLES
#undef FASTBAREAD

#ifdef CRA
#define DYNAMIC_ALLOC
#define CRAY_THREAD
#define ASSIGN
#define RFFTMLT
#define MINV
#define GETENV
#define SGEMVX1
#define SGERX1
#define FL2I
#define RANF
#define CRAY_BUFRLIB
#endif

#ifdef T90
#define DYNAMIC_ALLOC
#define CRAY_THREAD
#define ASSIGN
#define PXFGETENV
#define RANF
#define CRAY_BUFRLIB
#endif

#ifdef T3E
#define DYNAMIC_ALLOC
#define ASSIGN
#define RANF
#define CRAY_BUFRLIB
#endif

#ifdef IBMSP
#define REAL4_W3LIB
#define DYNAMIC_ALLOC
#define OPENMP
#undef  DCRFT
#define FASTBAREAD
#define flush flush_
#ifdef RMP
#undef DCRFT
#endif
#endif

#ifdef ORIGIN
#define REAL4_W3LIB
#define DYNAMIC_ALLOC
#define ORIGIN_THREAD
#endif

#ifdef X1E
#define OPENMP
#endif

#ifdef SGI
#define REAL4_W3LIB
#endif

#ifdef DEC
#define REAL4_W3LIB
#define FASTBAREAD
#define OPENMP
#endif

#ifdef LINUX_PGI
#define REAL4_W3LIB
#define OPENMP
#define FASTBAREAD
#endif

#ifdef LINUX_INTEL
#define REAL4_W3LIB
#define OPENMP
#define FASTBAREAD
#endif

#ifdef LINUX_GNU
#define REAL4_W3LIB
#define OPENMP
#define FASTBAREAD
#endif

#ifdef MAC
#define REAL4_W3LIB
#define DYNAMIC_ALLOC
#endif

#ifdef NEC
#define DYNAMIC_ALLOC
#define RFFTMLT
#define FASTBAREAD
#endif

#ifdef SX6
#define DYNAMIC_ALLOC
#define RFFTMLT
#define FASTBAREAD
#endif

#ifdef ES
#define DYNAMIC_ALLOC
#define ASLES
#define FASTBAREAD
#ifdef MP
#endif
#endif

#ifdef X1E
#define REAL4_W3LIB
#undef FASTBAREAD
#undef ASSIGN
#undef RANF
#endif
