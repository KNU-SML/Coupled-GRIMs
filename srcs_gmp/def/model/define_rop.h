#include <machine.h>
/* ----- grims-rop  setup ----- */

#define ROP        /* define rop */
#define RMP        /* define rmp */

#define SQK        /* spectral solver quick but more memory */
#undef SPT         /* spectral solver slow  but less memory */
#define G2R        /* rmp with global   large scale forcing */
#undef C2R         /* rmp with regional large scale forcing */

/**---------------------------------------------------------*/
/** Options for Model Configuration                         */
/**---------------------------------------------------------*/
/* roms ocean model vertical number of levels */
#define ROMS_LEVS 30
/* ocean model bathimetry file */
#define BATH_FILE topo_R2.nc    /* low resolution bathimetry.  High resolution is etopo2.nc */
#define SODA
#undef  LEVITUS
#define R2FRC
#undef  CaRD10V2
#ifdef  SODA
#define ROMS_STARTUP "SODA"
#define ROMS_INI_NLI 40
#elif defined LEVITUS
#define ROMS_STARTUP "LEVITUS"
#define ROMS_INI_NLI 33
#else
#define ROMS_STARTUP "GROMS"
#define ROMS_INI_NLI 20
#endif
#ifdef  SODA
#define ROMS_BNY_NLI 40
#else
#define ROMS_BNY_NLI 20
#endif
#define ROMS_INI_THETA_S 4
#define ROMS_INI_THETA_B 0.4
#define ROMS_INI_TCLINE  75
#define ROMS_INI_TSTART  15
#define ROMS_INI_TEND    15
#define ROMS_INI_SMTHINP 1
#define ROMS_BNY_THETA_S 4
#define ROMS_BNY_THETA_B 0.4
#define ROMS_BNY_TCLINE  75
#define ROMS_BNY_TSTART  15
#define ROMS_BNY_TEND    15
#define ROMS_BNY_SMTHINP 1
#define ROMS_FOR_NLI 20
#define ROMS_FOR_THETA_S 4
#define ROMS_FOR_THETA_B 0.4
#define ROMS_FOR_TCLINE  75
#define ROMS_FOR_TSTART  15
#define ROMS_FOR_TEND    15
#define ROMS_FOR_SMTHINP 1
#define ROMS_THETA_S 4
#define ROMS_THETA_B 0.4
#define ROMS_TCLINE  75

/**----------------------------------------------------------*/
/** Options for Physics and Output Configuration             */
/**----------------------------------------------------------*/
/* output options */
#undef  DIAGNOSTICS_UV      /* define if writing out momentum diagnostics */
#undef  DIAGNOSTICS_TS      /* define if writing out tracer diagnostics */
#define OUT_DOUBLE
#define AVERAGES
#define AVERAGES_FLUXES
/* advection and coriolis */
#define UV_ADV
#define UV_COR
#define TS_U3HADVECTION
#define TS_C4VADVECTION
 
/* pressure gradient */
#define DJ_GRADPS       /* Splines density  Jacobian (Shchepetkin, 2000) */
/* TIDAL FORCING */
#undef TIDAL
#ifdef TIDAL
#define SSH_TIDES       /* turn on computation of tidal elevation */
#define UV_TIDES        /* turn on computation of tidal currents */
#define ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#define ADD_M2OBC       /* Add tidal currents  to processed OBC data */
#endif

/* diffusion and viscosities schemes */
#define UV_VIS2
#define TS_DIF2
#define DIFF_GRID       /* Diffusivity coefficient scaled by grid size */
#define VISC_GRID       /* viscosity coefficient scaled by grid size */
#define MIX_GEO_TS
#define MIX_S_UV
 
#define UV_QDRAG        /* turn ON or OFF quadratic bottom friction */
 
/* vertical mixing schemes */
# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
 
/* equation of state - relaxation - salinity - coordinates - etc */
#define NONLIN_EOS
#define SOLAR_SOURCE

#define SALINITY

#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES
/* boundary conditions */
#define BRYFILE
#ifndef BRYFILE
#define M2CLIMATOLOGY
#define M3CLIMATOLOGY
#define TCLIMATOLOGY
#undef  ZCLIMATOLOGY
#define M2CLM_NUDGING
#define M3CLM_NUDGING
#define TCLM_NUDGING
#endif
#ifdef TIDAL
#define SOUTH_M2FLATHER
#define SOUTH_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_FSCHAPMAN
#define NORTH_M2FLATHER
#define NORTH_FSCHAPMAN
#undef  WEST_VOLCONS
#undef  SOUTH_VOLCONS
#endif

#define RADIATION_2D
#define SPONGE
/* Climats processing and relaxations */
#define QCORRECTION
#define SCORRECTION
#undef  DIURNAL_SRFLUX  /* impose shortwave radiation local diurnal cycle */
#define ANA_BSFLUX
#define ANA_BTFLUX

/*
/**-------------------------------------------------------------
/** Options for Model Opening Boundary Conditions.
/**-------------------------------------------------------------
/*/
#undef  CLOSED_OBC
#ifdef  CLOSED_OBC
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#else
#undef  SOUTHERN_WALL
#define  WESTERN_WALL
#undef  NORTHERN_WALL
#undef  EASTERN_WALL
#endif
#ifndef SOUTHERN_WALL
#define SOUTH_TNUDGING
#define SOUTH_M3NUDGING
#define SOUTH_M2NUDGING
#define SOUTH_M2FLATHER
#define SOUTH_FSCHAPMAN
#define SOUTH_TRADIATION
#define SOUTH_M3RADIATION
#endif
#ifndef WESTERN_WALL
#define WEST_TNUDGING
#define WEST_M3NUDGING
#define WEST_M2NUDGING
#define WEST_M2FLATHER
#define WEST_FSCHAPMAN
#define WEST_TRADIATION
#define WEST_M3RADIATION
#endif
#ifndef NORTHERN_WALL
#define NORTH_TNUDGING
#define NORTH_M3NUDGING
#define NORTH_M2NUDGING
#define NORTH_M2FLATHER
#define NORTH_FSCHAPMAN
#define NORTH_TRADIATION
#define NORTH_M3RADIATION
#endif
#ifndef EASTERN_WALL
#define EAST_TNUDGING
#define EAST_M3NUDGING
#define EAST_M2NUDGING
#define EAST_M2FLATHER
#define EAST_FSCHAPMAN
#define EAST_TRADIATION
#define EAST_M3RADIATION
#endif

#define MPI
#define COUPLE_ROP
#define MPI_COMM_ROMS mpi_comm_roms
#define MPI_COMM_RMP mpi_comm_rmp

