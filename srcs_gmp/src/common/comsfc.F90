#include "define.h"
   module comsfc
!-------------------------------------------------------------------------------
!
!  common for variables in the surface file
!
!  important:  do not change the order of the variables
!              do not insert any variable, simply add them
!              at the end
!
!-------------------------------------------------------------------------------
   character(len=4)  ::  sfcftyp
#ifdef OSULSM1
   data sfcftyp/'osu1'/
#endif
#ifdef OSULSM2
   data sfcftyp/'osu2'/
#endif
#ifdef NOALSM1
   data sfcftyp/'noa1'/
#endif
#ifdef VICLSM1
   data sfcftyp/'vic1'/
#endif

   real,allocatable,target,dimension(:,:,:)     :: sfcv
   real,pointer,dimension(:,:)      ::  sfcfcs
   real,pointer,dimension(:,:)      ::  z0cmt,tml,t0ml,hml,h0ml,huml,hvml,dtw1
   real,pointer,dimension(:,:)      ::  hml0,tmoml
   real,pointer,dimension(:,:)      ::  kprfi
   real,pointer,dimension(:,:,:)    ::  paer,denni,idxci,cmixi
#ifdef OSULSM1
   real,pointer,dimension(:,:)      ::  tsea,snoweq,tg3,z0cm,cv,cvb,cvt,       &
                                        slmsk,plantr,canopy,f10m
   real,pointer, dimension(:,:,:)   ::  smc,stc,albedo
#endif
#ifdef OSULSM2
   real,pointer, dimension(:,:)     ::  tsea,snoweq,tg3,z0cm,cv,cvb,cvt,       &
                                        slmsk,vfrac,canopy,f10m,vtype,stype,   &
                                        uustar,ffmm,ffhh
   real,pointer, dimension(:,:,:)   ::  smc,stc,albedo,facalf
#endif
#ifdef NOALSM1
   real,pointer, dimension(:,:)     ::  tsea,snoweq,tg3,z0cm,cv,cvb,cvt,       &
                                        slmsk,vfrac,canopy,f10m,vtype,stype,   &
                                        uustar,ffmm,ffhh,prcp,srflag,snwdph,   &
                                        shdmin,shdmax,slope,snoalb
   real,pointer, dimension(:,:,:)   ::  smc,stc,albedo,facalf,slc
#endif
#ifdef VICLSM1
   real,pointer, dimension(:,:)     ::  tsea,snoweq,tg3,z0cm,cv,cvb,cvt,       &
                                        slmsk,vfrac,canopy,f10m,vtype,         &
                                        uustar,ffmm,ffhh,prcp,srflag,binf,     &
                                        ds,dsm,ws,cef,flai,silz,snwz,csno,     &
                                        rsno,tsf,tpk,sfw,pkw,lstsn
   real,pointer, dimension(:,:,:)   ::  smc,stc,albedo,vroot,facalf,expt,      &
                                        kest,dph,bub,qrt,bkd,sld,wcr,wpw,smr,  &
                                        smx,dphn,smxn,expn,bubn,alpn,betn,     &
                                        gamn,sic
#endif
#ifdef VICLSM2
   real,pointer, dimension(:,:)     ::  tsea,snoweq,tg3,z0cm,cv,cvb,cvt,       &
                                        slmsk,vfrac,canopy,f10m,vtype,         &
                                        uustar,ffmm,ffhh,prcp,srflag,binf,     &
                                        ds,dsm,ws,cef,silz,snwz,nveg
   real,pointer, dimension(:,:,:)   ::  smc,stc,albedo,vroot,facalf,expt,      &
                                        kest,dph,bub,qrt,bkd,sld,wcr,wpw,smr,  &
                                        smx,dphn,smxn,expn,bubn,alpn,betn,gamn,&
                                        mvfr,mcnp,mvty,flai,msno,msmc,msic,    &
                                        mstc,csno,rsno,tsf,tpk,sfw,pkw,lstsn
#endif
!
   contains
!-------------------------------------------------------------------------------
   subroutine init_comsfc(LONF2S,LATG2S)
!-------------------------------------------------------------------------------
   use varsfc, only : lsoil=>lsoil_,nsoil=>nsoil_,lalbd=>lalbd_,msub=>msub_,   &
                      kslmb=>kslmb_,nslmb=>nslmb_,numsfcs,                     &
                      naer=>naer_,nden=>nden_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                                      ::  LONF2S, LATG2S
   integer                                      ::  ns
!
! note. numsfcs is smaller than allocated size by 8 because of z0cmt
!       additional 7 variables are added for ocean mixed layer, 
!       tml, t0ml, hml, h0ml, huml, hvml   (hml0 is cycling)
!
   allocate(sfcv(LONF2S,LATG2S,numsfcs+9))
!
#ifdef OSULSM1 /* numsfcs = 11+lsoil_*2+lalbd_ */
   tsea  =>sfcv(1:LONF2S,1:LATG2S,1)
   smc   =>sfcv(1:LONF2S,1:LATG2S,2:1+lsoil)
   snoweq=>sfcv(1:LONF2S,1:LATG2S,2+  lsoil)
   stc   =>sfcv(1:LONF2S,1:LATG2S,3+  lsoil:2+2*lsoil)
   tg3   =>sfcv(1:LONF2S,1:LATG2S,3+2*lsoil)
   z0cm  =>sfcv(1:LONF2S,1:LATG2S,4+2*lsoil)
   cv    =>sfcv(1:LONF2S,1:LATG2S,5+2*lsoil)
   cvb   =>sfcv(1:LONF2S,1:LATG2S,6+2*lsoil)
   cvt   =>sfcv(1:LONF2S,1:LATG2S,7+2*lsoil)
   albedo=>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil:7+2*lsoil+lalbd)
   slmsk =>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil+lalbd)
   plantr=>sfcv(1:LONF2S,1:LATG2S,9+2*lsoil+lalbd)
   canopy=>sfcv(1:LONF2S,1:LATG2S,10+2*lsoil+lalbd)
   f10m  =>sfcv(1:LONF2S,1:LATG2S,11+2*lsoil+lalbd)
!  hml0  =>sfcv(1:LONF2S,1:LATG2S,12+2*lsoil+lalbd)
#endif
#ifdef OSULSM2 /* numsfcs = 20+lsoil_*2+lalbd_ */
   tsea  =>sfcv(1:LONF2S,1:LATG2S,1)
   smc   =>sfcv(1:LONF2S,1:LATG2S,2:1+lsoil)
   snoweq=>sfcv(1:LONF2S,1:LATG2S,2+  lsoil)
   stc   =>sfcv(1:LONF2S,1:LATG2S,3+  lsoil:2+2*lsoil)
   tg3   =>sfcv(1:LONF2S,1:LATG2S,3+2*lsoil)
   z0cm  =>sfcv(1:LONF2S,1:LATG2S,4+2*lsoil)
   cv    =>sfcv(1:LONF2S,1:LATG2S,5+2*lsoil)
   cvb   =>sfcv(1:LONF2S,1:LATG2S,6+2*lsoil)
   cvt   =>sfcv(1:LONF2S,1:LATG2S,7+2*lsoil)
   albedo=>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil:7+2*lsoil+lalbd)
   slmsk =>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil+lalbd)
   vfrac =>sfcv(1:LONF2S,1:LATG2S,9+2*lsoil+lalbd)
   canopy=>sfcv(1:LONF2S,1:LATG2S,10+2*lsoil+lalbd)
   f10m  =>sfcv(1:LONF2S,1:LATG2S,11+2*lsoil+lalbd)
   vtype =>sfcv(1:LONF2S,1:LATG2S,12+2*lsoil+lalbd)
   stype =>sfcv(1:LONF2S,1:LATG2S,13+2*lsoil+lalbd)
   facalf=>sfcv(1:LONF2S,1:LATG2S,14+2*lsoil+lalbd:15+2*lsoil+lalbd)
   uustar=>sfcv(1:LONF2S,1:LATG2S,16+2*lsoil+lalbd)
   ffmm  =>sfcv(1:LONF2S,1:LATG2S,17+2*lsoil+lalbd)
   ffhh  =>sfcv(1:LONF2S,1:LATG2S,18+2*lsoil+lalbd)
   hml0  =>sfcv(1:LONF2S,1:LATG2S,19+2*lsoil+lalbd)
   paer  =>sfcv(1:LONF2S,1:LATG2S,20+2*lsoil+lalbd:19+2*lsoil+lalbd+naer)
   kprfi =>sfcv(1:LONF2S,1:LATG2S,20+2*lsoil+lalbd+naer)
   denni =>sfcv(1:LONF2S,1:LATG2S,21+2*lsoil+lalbd+naer:20+2*lsoil+lalbd+naer+nden)
   ns = 20+2*lsoil+lalbd+naer+nden
   idxci =>sfcv(1:LONF2S,1:LATG2S,1+ns:ns+naer)
   cmixi =>sfcv(1:LONF2S,1:LATG2S,1+ns+naer:ns+2*naer)
#endif
#ifdef NOALSM1 /* numsfcs = 27+lsoil_*3+lalbd_ */
   tsea  =>sfcv(1:LONF2S,1:LATG2S,1)
   smc   =>sfcv(1:LONF2S,1:LATG2S,2:1+lsoil)
   snoweq=>sfcv(1:LONF2S,1:LATG2S,2+  lsoil)
   stc   =>sfcv(1:LONF2S,1:LATG2S,3+  lsoil:2+2*lsoil)
   tg3   =>sfcv(1:LONF2S,1:LATG2S,3+2*lsoil)
   z0cm  =>sfcv(1:LONF2S,1:LATG2S,4+2*lsoil)
   cv    =>sfcv(1:LONF2S,1:LATG2S,5+2*lsoil)
   cvb   =>sfcv(1:LONF2S,1:LATG2S,6+2*lsoil)
   cvt   =>sfcv(1:LONF2S,1:LATG2S,7+2*lsoil)
   albedo=>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil:7+2*lsoil+lalbd)
   slmsk =>sfcv(1:LONF2S,1:LATG2S,8+2*lsoil+lalbd)
   vfrac =>sfcv(1:LONF2S,1:LATG2S,9+2*lsoil+lalbd)
   canopy=>sfcv(1:LONF2S,1:LATG2S,10+2*lsoil+lalbd)
   f10m  =>sfcv(1:LONF2S,1:LATG2S,11+2*lsoil+lalbd)
   vtype =>sfcv(1:LONF2S,1:LATG2S,12+2*lsoil+lalbd)
   stype =>sfcv(1:LONF2S,1:LATG2S,13+2*lsoil+lalbd)
   facalf=>sfcv(1:LONF2S,1:LATG2S,14+2*lsoil+lalbd:15+2*lsoil+lalbd)
   uustar=>sfcv(1:LONF2S,1:LATG2S,16+2*lsoil+lalbd)
   ffmm  =>sfcv(1:LONF2S,1:LATG2S,17+2*lsoil+lalbd)
   ffhh  =>sfcv(1:LONF2S,1:LATG2S,18+2*lsoil+lalbd)
   prcp  =>sfcv(1:LONF2S,1:LATG2S,19+2*lsoil+lalbd)
   srflag=>sfcv(1:LONF2S,1:LATG2S,20+2*lsoil+lalbd)
   snwdph=>sfcv(1:LONF2S,1:LATG2S,21+2*lsoil+lalbd)
   slc   =>sfcv(1:LONF2S,1:LATG2S,22+2*lsoil+lalbd:21+3*lsoil+lalbd)
   shdmin=>sfcv(1:LONF2S,1:LATG2S,22+3*lsoil+lalbd)
   shdmax=>sfcv(1:LONF2S,1:LATG2S,23+3*lsoil+lalbd)
   slope =>sfcv(1:LONF2S,1:LATG2S,24+3*lsoil+lalbd)
   snoalb=>sfcv(1:LONF2S,1:LATG2S,25+3*lsoil+lalbd)
   hml0  =>sfcv(1:LONF2S,1:LATG2S,26+3*lsoil+lalbd)
   paer  =>sfcv(1:LONF2S,1:LATG2S,27+3*lsoil+lalbd:26+3*lsoil+lalbd+naer)
   kprfi =>sfcv(1:LONF2S,1:LATG2S,27+3*lsoil+lalbd+naer)
   denni =>sfcv(1:LONF2S,1:LATG2S,28+3*lsoil+lalbd+naer:27+3*lsoil+lalbd+naer+nden)
   ns = 27+3*lsoil+lalbd+naer+nden
   idxci =>sfcv(1:LONF2S,1:LATG2S,1+ns:ns+naer)
   cmixi =>sfcv(1:LONF2S,1:LATG2S,1+ns+naer:ns+2*naer)
#endif
#ifdef VICLSM1  /* numsfcs=36+lsoil_*14+nsoil_*8+lalbd_ */
   tsea  =>sfcv(1:LONF2S,1:LATG2S,1)
   smc   =>sfcv(1:LONF2S,1:LATG2S,2:1+lsoil)
   snoweq=>sfcv(1:LONF2S,1:LATG2S,2+  lsoil)
   stc   =>sfcv(1:LONF2S,1:LATG2S,3+  lsoil:2+  lsoil+  nsoil)
   ns=3+lsoil+  nsoil
   tg3   =>sfcv(1:LONF2S,1:LATG2S,ns)
   z0cm  =>sfcv(1:LONF2S,1:LATG2S,1+ns)
   cv    =>sfcv(1:LONF2S,1:LATG2S,2+ns)
   cvb   =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   cvt   =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   albedo=>sfcv(1:LONF2S,1:LATG2S,5+ns:4+ns+lalbd)
   ns=5+ns+lalbd
   slmsk =>sfcv(1:LONF2S,1:LATG2S,ns)
   vfrac =>sfcv(1:LONF2S,1:LATG2S,1+ns)
   canopy=>sfcv(1:LONF2S,1:LATG2S,2+ns)
   f10m  =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   vtype =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   vroot =>sfcv(1:LONF2S,1:LATG2S,5+ns:4+ns+lsoil)
   ns=5+ns+lsoil
   facalf=>sfcv(1:LONF2S,1:LATG2S,ns:1+ns)
   uustar=>sfcv(1:LONF2S,1:LATG2S,2+ns)
   ffmm  =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   ffhh  =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   prcp  =>sfcv(1:LONF2S,1:LATG2S,5+ns)
   srflag=>sfcv(1:LONF2S,1:LATG2S,6+ns)
   binf  =>sfcv(1:LONF2S,1:LATG2S,7+ns)
   ds    =>sfcv(1:LONF2S,1:LATG2S,8+ns)
   dsm   =>sfcv(1:LONF2S,1:LATG2S,9+ns)
   ws    =>sfcv(1:LONF2S,1:LATG2S,10+ns)
   cef   =>sfcv(1:LONF2S,1:LATG2S,11+ns)
   expt  =>sfcv(1:LONF2S,1:LATG2S,12+ns:11+ns+lsoil)
   ns=12+ns+lsoil
   kest  =>sfcv(1:LONF2S,1:LATG2S,ns:ns-1+lsoil)
   dph   =>sfcv(1:LONF2S,1:LATG2S,ns+lsoil:ns-1+2*lsoil)
   bub   =>sfcv(1:LONF2S,1:LATG2S,ns+2*lsoil:ns+3*lsoil-1)
   qrt   =>sfcv(1:LONF2S,1:LATG2S,ns+3*lsoil:ns+4*lsoil-1)
   bkd   =>sfcv(1:LONF2S,1:LATG2S,ns+4*lsoil:ns+5*lsoil-1)
   sld   =>sfcv(1:LONF2S,1:LATG2S,ns+5*lsoil:ns+6*lsoil-1)
   wcr   =>sfcv(1:LONF2S,1:LATG2S,ns+6*lsoil:ns+7*lsoil-1)
   wpw   =>sfcv(1:LONF2S,1:LATG2S,ns+7*lsoil:ns+8*lsoil-1)
   smr   =>sfcv(1:LONF2S,1:LATG2S,ns+8*lsoil:ns+9*lsoil-1)
   smx   =>sfcv(1:LONF2S,1:LATG2S,ns+9*lsoil:ns+10*lsoil-1)
   dphn  =>sfcv(1:LONF2S,1:LATG2S,ns+10*lsoil:ns+10*lsoil+nsoil-1)
   ns=ns+10*lsoil+nsoil
   smxn  =>sfcv(1:LONF2S,1:LATG2S,ns        :ns+  nsoil-1)
   expn  =>sfcv(1:LONF2S,1:LATG2S,ns+  nsoil:ns+2*nsoil-1)
   bubn  =>sfcv(1:LONF2S,1:LATG2S,ns+2*nsoil:ns+3*nsoil-1)
   alpn  =>sfcv(1:LONF2S,1:LATG2S,ns+3*nsoil:ns+4*nsoil-1)
   betn  =>sfcv(1:LONF2S,1:LATG2S,ns+4*nsoil:ns+5*nsoil-1)
   gamn  =>sfcv(1:LONF2S,1:LATG2S,ns+5*nsoil:ns+6*nsoil-1)
   ns=ns+6*nsoil
   flai  =>sfcv(1:LONF2S,1:LATG2S,ns)
   silz  =>sfcv(1:LONF2S,1:LATG2S,ns+1)
   snwz  =>sfcv(1:LONF2S,1:LATG2S,ns+2)
   sic   =>sfcv(1:LONF2S,1:LATG2S,ns+3:ns+3+lsoil-1)
   ns=ns+3+lsoil
   csno  =>sfcv(1:LONF2S,1:LATG2S,ns  )
   rsno  =>sfcv(1:LONF2S,1:LATG2S,ns+1)
   tsf   =>sfcv(1:LONF2S,1:LATG2S,ns+2)
   tpk   =>sfcv(1:LONF2S,1:LATG2S,ns+3)
   sfw   =>sfcv(1:LONF2S,1:LATG2S,ns+4)
   pkw   =>sfcv(1:LONF2S,1:LATG2S,ns+5)
   lstsn =>sfcv(1:LONF2S,1:LATG2S,ns+6)
   hml0  =>sfcv(1:LONF2S,1:LATG2S,ns+7)
   paer  =>sfcv(1:LONF2S,1:LATG2S,ns+8:ns+8+naer-1)
   ns=ns+8+naer
   kprfi =>sfcv(1:LONF2S,1:LATG2S,ns+1)
   denni =>sfcv(1:LONF2S,1:LATG2S,ns+2:ns+2+nden-1)
   ns=ns+2+nden
   idxci =>sfcv(1:LONF2S,1:LATG2S,1+ns:ns+naer)
   cmixi =>sfcv(1:LONF2S,1:LATG2S,1+ns+naer:ns+2*naer)
#endif
#ifdef VICLSM2 /* numsfcs = 29+lsoil*12+nsoil*8+msub*11+kslmb*4+nslmb+lalbd */
   tsea  =>sfcv(1:LONF2S,1:LATG2S,1)
   smc   =>sfcv(1:LONF2S,1:LATG2S,2:1+lsoil)
   snoweq=>sfcv(1:LONF2S,1:LATG2S,2+  lsoil)
   stc   =>sfcv(1:LONF2S,1:LATG2S,3+  lsoil:2+  lsoil+  nsoil)
   ns=3+lsoil+  nsoil
   tg3   =>sfcv(1:LONF2S,1:LATG2S,ns)
   z0cm  =>sfcv(1:LONF2S,1:LATG2S,1+ns)
   cv    =>sfcv(1:LONF2S,1:LATG2S,2+ns)
   cvb   =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   cvt   =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   albedo=>sfcv(1:LONF2S,1:LATG2S,5+ns:4+ns+lalbd)
   ns=5+ns+lalbd
   slmsk =>sfcv(1:LONF2S,1:LATG2S,ns)
   vfrac =>sfcv(1:LONF2S,1:LATG2S,1+ns)
   canopy=>sfcv(1:LONF2S,1:LATG2S,2+ns)
   f10m  =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   vtype =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   vroot =>sfcv(1:LONF2S,1:LATG2S,5+ns:4+ns+kslmb)
   ns=5+ns+kslmb
   facalf=>sfcv(1:LONF2S,1:LATG2S,ns:1+ns)
   uustar=>sfcv(1:LONF2S,1:LATG2S,2+ns)
   ffmm  =>sfcv(1:LONF2S,1:LATG2S,3+ns)
   ffhh  =>sfcv(1:LONF2S,1:LATG2S,4+ns)
   prcp  =>sfcv(1:LONF2S,1:LATG2S,5+ns)
   srflag=>sfcv(1:LONF2S,1:LATG2S,6+ns)
   binf  =>sfcv(1:LONF2S,1:LATG2S,7+ns)
   ds    =>sfcv(1:LONF2S,1:LATG2S,8+ns)
   dsm   =>sfcv(1:LONF2S,1:LATG2S,9+ns)
   ws    =>sfcv(1:LONF2S,1:LATG2S,10+ns)
   cef   =>sfcv(1:LONF2S,1:LATG2S,11+ns)
   expt  =>sfcv(1:LONF2S,1:LATG2S,12+ns:11+ns+lsoil)
   ns=12+ns+lsoil
   kest  =>sfcv(1:LONF2S,1:LATG2S,ns:ns-1+lsoil)
   dph   =>sfcv(1:LONF2S,1:LATG2S,ns+lsoil:ns-1+2*lsoil)
   bub   =>sfcv(1:LONF2S,1:LATG2S,ns+2*lsoil:ns+3*lsoil-1)
   qrt   =>sfcv(1:LONF2S,1:LATG2S,ns+3*lsoil:ns+4*lsoil-1)
   bkd   =>sfcv(1:LONF2S,1:LATG2S,ns+4*lsoil:ns+5*lsoil-1)
   sld   =>sfcv(1:LONF2S,1:LATG2S,ns+5*lsoil:ns+6*lsoil-1)
   wcr   =>sfcv(1:LONF2S,1:LATG2S,ns+6*lsoil:ns+7*lsoil-1)
   wpw   =>sfcv(1:LONF2S,1:LATG2S,ns+7*lsoil:ns+8*lsoil-1)
   smr   =>sfcv(1:LONF2S,1:LATG2S,ns+8*lsoil:ns+9*lsoil-1)
   smx   =>sfcv(1:LONF2S,1:LATG2S,ns+9*lsoil:ns+10*lsoil-1)
   dphn  =>sfcv(1:LONF2S,1:LATG2S,ns+10*lsoil:ns+10*lsoil+nsoil-1)
   ns=ns+10*lsoil+nsoil
   smxn  =>sfcv(1:LONF2S,1:LATG2S,ns        :ns+  nsoil-1)
   expn  =>sfcv(1:LONF2S,1:LATG2S,ns+  nsoil:ns+2*nsoil-1)
   bubn  =>sfcv(1:LONF2S,1:LATG2S,ns+2*nsoil:ns+3*nsoil-1)
   alpn  =>sfcv(1:LONF2S,1:LATG2S,ns+3*nsoil:ns+4*nsoil-1)
   betn  =>sfcv(1:LONF2S,1:LATG2S,ns+4*nsoil:ns+5*nsoil-1)
   gamn  =>sfcv(1:LONF2S,1:LATG2S,ns+5*nsoil:ns+6*nsoil-1)
   ns=ns+6*nsoil
   silz  =>sfcv(1:LONF2S,1:LATG2S,ns  )
   snwz  =>sfcv(1:LONF2S,1:LATG2S,ns+1)
   nveg  =>sfcv(1:LONF2S,1:LATG2S,ns+2)
   mvfr  =>sfcv(1:LONF2S,1:LATG2S,ns+3:ns+3+msub-1)
   ns=ns+3+msub
   mcnp  =>sfcv(1:LONF2S,1:LATG2S,ns       :ns+  msub-1)
   mvty  =>sfcv(1:LONF2S,1:LATG2S,ns+  msub:ns+2*msub-1)
   flai  =>sfcv(1:LONF2S,1:LATG2S,ns+2*msub:ns+3*msub-1)
   msno  =>sfcv(1:LONF2S,1:LATG2S,ns+3*msub:ns+4*msub-1)
   msmc  =>sfcv(1:LONF2S,1:LATG2S,ns+4*msub:ns+4*msub+kslmb-1)
   ns=ns+4*msub+kslmb
   msic  =>sfcv(1:LONF2S,1:LATG2S,ns:ns+kslmb-1)
   mstc  =>sfcv(1:LONF2S,1:LATG2S,ns+kslmb:ns+kslmb+nslmb-1)
   ns=ns+kslmb+nslmb
   csno  =>sfcv(1:LONF2S,1:LATG2S,ns       :ns+  msub-1)
   rsno  =>sfcv(1:LONF2S,1:LATG2S,ns+  msub:ns+2*msub-1)
   tsf   =>sfcv(1:LONF2S,1:LATG2S,ns+2*msub:ns+3*msub-1)
   tpk   =>sfcv(1:LONF2S,1:LATG2S,ns+3*msub:ns+4*msub-1)
   sfw   =>sfcv(1:LONF2S,1:LATG2S,ns+4*msub:ns+5*msub-1)
   pkw   =>sfcv(1:LONF2S,1:LATG2S,ns+5*msub:ns+6*msub-1)
   lstsn =>sfcv(1:LONF2S,1:LATG2S,ns+6*msub:ns+7*msub-1)
   hml0  =>sfcv(1:LONF2S,1:LATG2S,ns+7*msub:ns+7*msub+1)
   paer  =>sfcv(1:LONF2S,1:LATG2S,ns+8*msub:ns+8*msub+1+naer-1)
   ns=ns+8*msub+1+naer
   kprfi =>sfcv(1:LONF2S,1:LATG2S,ns:ns+1)
   denni =>sfcv(1:LONF2S,1:LATG2S,2+ns:2+ns+nden-1)
   idxci =>sfcv(1:LONF2S,1:LATG2S,2+ns+nden:1+ns+nden+naer)
   cmixi =>sfcv(1:LONF2S,1:LATG2S,2+ns+nden+naer:1+ns+nden+2*naer)
#endif
!
! thermal roughness lengh and ocean mixed layer vairables
! for inside model and output, not cycling
!
   z0cmt =>sfcv(1:LONF2S,1:LATG2S,numsfcs+1)
   tml   =>sfcv(1:LONF2S,1:LATG2S,numsfcs+2)
   t0ml  =>sfcv(1:LONF2S,1:LATG2S,numsfcs+3)
   hml   =>sfcv(1:LONF2S,1:LATG2S,numsfcs+4)
   h0ml  =>sfcv(1:LONF2S,1:LATG2S,numsfcs+5)
   huml  =>sfcv(1:LONF2S,1:LATG2S,numsfcs+6)
   hvml  =>sfcv(1:LONF2S,1:LATG2S,numsfcs+7)
   dtw1  =>sfcv(1:LONF2S,1:LATG2S,numsfcs+8)
   tmoml =>sfcv(1:LONF2S,1:LATG2S,numsfcs+9)
!
   sfcfcs=>sfcv(1:LONF2S*LATG2S,1,1:numsfcs)
!
   return
   end subroutine init_comsfc
!-------------------------------------------------------------------------------
   end module comsfc
