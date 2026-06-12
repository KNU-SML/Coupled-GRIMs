/* version 1 of grib headers  w. ebisuzaki */
/* this version is incomplete */

/* I did it in this manner because I was worried that the cray
would not handle structures gracefully */


struct PDS {
	unsigned char bytes[28];
};

#define PDSLen1		bytes[0]
#define PDSLen2		bytes[1]
#define PDSLen3		bytes[2]
#define PDSVsn		bytes[3]
#define Center		bytes[4]
#define Model		bytes[5]
#define Grid		bytes[6]
#define GDS_BMS		bytes[7]
#define PARAM		bytes[8]
#define L_TYPE		bytes[9]
#define LEVEL1		bytes[10]
#define LEVEL2		bytes[11]

#define KPDS5(pds)	(pds->bytes[8])
#define KPDS6(pds)	(pds->bytes[9])
#define KPDS7(pds)	(pds->bytes[10]*256 + pds->bytes[11])

#define Year		bytes[12]
#define Month		bytes[13]
#define Day		bytes[14]
#define Hour		bytes[15]
#define Minute		bytes[16]
#define ForecastTimeUnit	bytes[17]
#define P1		bytes[18]
#define P2		bytes[19]
#define TimeRange	bytes[20]
#define Century		bytes[24]
#define Subcenter	bytes[25]
