#if defined(SCM) || !defined(MP)
#define MPABORT abort
#else
#ifdef RMP
#define MPABORT rmpabort
#else
#define MPABORT mpabort
#endif
#endif
