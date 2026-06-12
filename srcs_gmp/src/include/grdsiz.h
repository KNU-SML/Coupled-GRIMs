#ifndef RMP
#define LONF2F lonf2_
#define LATG2F latg2_
#ifdef MP
#define LONF2S lonf2p_
#define LONF22S lonf22p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LONF22S lonf22_
#define LATG2S latg2_
#endif
#else
#define LONF2F igrd12_
#define LATG2F jgrd12_
#ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif

