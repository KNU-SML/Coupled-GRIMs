#include <machine.h>
/* ----- grims-nim non-hydrostatic icosahedral model (nim) setup -----  */
/* NIM options */
#define NIM          /* for NIM */
#undef  NIM_DIAG     /* diagnostics dbg */
#undef NIM_DBG       /* NIM debugging option */
#undef NIM_DBGM      /* NIM debugging option for mpi */
#define NIMAQUA      /* for aqua simulation */
#undef GFS_SFC         /* read gfs sfc data */
#undef LANDSEA       /* read gfs sfc data */
