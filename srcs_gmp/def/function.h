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
