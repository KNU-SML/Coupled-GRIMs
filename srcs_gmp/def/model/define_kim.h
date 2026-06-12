#include <machine.h>
/* ----- grims-kim non-hydrostatic icosahedral model (kim) setup -----  */
/* KIM options */
#define KIM         /* for KIM */
/* NIM options */
#define NIM         /* for NIM (KIM in NIM) */
#undef NIM_OMP      /* turn on WRF-style "tile" loops for OpenMP in GRIMs */
#undef NIM_DIAG     /* diagnostics dbg */
#undef NIM_DBG      /* NIM debugging option */
#undef NIM_DBGM     /* NIM debugging option for mpi */
#undef NIM_BIN      /* NIM binary plotting */
#define NIMAQUA     /* for aqua simulation */
#undef NIM_AEROSOL  /* read aerosol data directly for nim prep */
#define GFS_SFC     /* read gfs sfc data */
#undef LANDSEA      /* read gfs sfc data */
