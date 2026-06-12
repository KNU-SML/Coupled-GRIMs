#include <machine.h>
/* ----- grims-rmp  setup ----- */

#undef GDAS        /* GDAS option off */
#define RMP        /* define rmp */

#define SQK        /* spectral solver quick but more memory */
#undef SPT         /* spectral solver slow  but less memory */
#define G2R        /* rmp with global   large scale forcing */
#undef C2R         /* rmp with regional large scale forcing */
#undef IDEAL       /* idealized simulations -> undef g2r, define c2r */
#undef COLD        /* idealized 2D cold bubble case */
#undef WARM        /* idealized 2D warm bubble case */
#undef HILL        /* idealized 2D flow case over a bell-shaped hill */
