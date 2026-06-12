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
#undef NCEP2010
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
