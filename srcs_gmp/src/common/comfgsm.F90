!
   module comfgsm
!-------------------------------------------------------------------------------
   integer  ::  lnts2,lnoffset,lons2,lats2,n1                                 ,&
                kpfix,ksfcx,komlx,ksig,ksfc,kpost,krestart,klfm,krsm 
#ifdef RMP
   integer  ::  nrsmi1,nrsmi2,nrflip                                          ,&
                nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken
#endif
#ifdef LFM
   integer  ::  nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp
#endif
   real     ::  solsec,dthr,hdthr
#ifdef LFM
   real     ::  weight,filtwin
#endif
!
   end module comfgsm
