#include <machine.h>
/* ----- grims-smp single column model setup ----- */

#undef GDAS        /* GDAS option off */
#define SMP        /* single column model program option */

#define CHOOSE_OBS /* observation data options */
#undef EXTRA1_OBS  /* observation data options */
#undef EXTRA2_OBS  /* observation data options */
#undef EXTRA3_OBS  /* observation data options */
#define EXPLICIT_CLOUDINESS /* qci qrs from convective clouds output */
#undef SMP_NUDGING_SFC      /* surface variable nudging option */
#undef SMP_NUDGING_FLX      /* surface flux nudging option */
#undef SMP_NUDGING_RAD      /* radiation flux nudging option */
#undef SMP_RA2SFC  /* surface data from ra2 sfc cycling when no osu1 lsm */
 
#ifdef ARM         /* arm case */
#define SGPyy      /* sgp year (yy) */
#undef  REV_FRC    /* relaxation parameter for winds */
#undef  HOR_FRC    /* relaxation parameter for horizontal advection */ 
#define RLX_FRC    /* relaxation to observations */
#ifdef SGPyy
#define _gxlon_ 262.51 /* longitude point */
#define _gylat_ 36.61  /* latitude point */
#define _lsmsk_ 1      /* land sea mask */
#define _lobsl_ 0      /* ceop option */
#endif
#endif

#ifdef GATE            /* gate case */
#define _gxlon_ 336.5  /* longitude point */
#define _gylat_ 9.0    /* latitude point */
#define _lsmsk_ 0      /* land sea mask */
#define _lobsl_ 0      /* ceop option */
#define RLX_FRC        /* relaxation to observations */
#endif

#ifdef TOGA
#define _gxlon_ 156.0
#define _gylat_ -2.0
#define _lsmsk_ 0
#define _lobsl_ 0
#undef  REV_FRC
#undef  HOR_FRC
#define RLX_FRC
#endif

#ifdef TRMM
#define _gxlon_ 167.4
#define _gylat_ 8.6
#define _lsmsk_ 0
#define _lobsl_ 0
#undef  REV_FRC
#undef  HOR_FRC
#define RLX_FRC
#endif

