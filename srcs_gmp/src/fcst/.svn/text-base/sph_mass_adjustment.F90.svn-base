#include <define.h>
   subroutine sph_mass_adjustment(q,avprs0)                                               
!-------------------------------------------------------------------------------
!
! subroutine: sph_mass_adjustment
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
!                                                                               
!  surface pressure correction for long integrations                            
!                                                                               
   use paramodel, only :  JCAP1S,LNT22S,LLN22S,LCAP22S,                        &
                          lnt2_,lonf22_,latg2_,                                &
#ifdef MP
                          lonf_,jcap1_,lcap22p_
#else
                          lonf_,jcap1_
#endif
   use comio
   use comfgrid
   use comgpln
#ifdef REDUCE_GRID
   use comreduce
#endif
#ifdef MP
   use commpi
   use paramodel, only   :  lln22p_,latg2p_,lnt22p_
   real                 ::  qa(lln22p_),syf(lonf22_,latg2p_),avprsjp(latg2p_)
#endif
!                                                                               
   real,intent(inout)   ::  q(LNT22S)
   real,intent(inout)   ::  avprs0
!
!  local
!
   real                 ::  prs(LCAP22S,latg2_)
   real                 ::  qnew(LLN22S,latg2_)
   real                 ::  flp(2,JCAP1S,latg2_),flm(2,JCAP1S,latg2_)
   real                 ::  avprsj(latg2_)
!=============================================================
#ifdef MP
   llstr=lwvstr(mype)
   llens=lwvlen(mype)
!                                                                               
!  loop to find average pressure                                                
!                                                                               
   avprs=0.                                                                  
   sumwgt=0.                                                                 
!
   call mpnn2n (q,lnt22p_,qa,lln22p_,1)
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(qtt,colrad,qa,prs,llstr,llens,lwvdef,lcapd,lcapdp,mype)
!$& local(lat,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,llensd)
#endif
   do lat = 1,latg2_
#ifdef REDUCE_GRID
     llensd=lcapdp(lat,mype)
#else
     llensd=llens
#endif
     call sph_sum_coeff(qa,prs(1,lat),qtt(1,lat),llstr,llensd,lwvdef,1)
   enddo
!
   call mpnl2ny(prs,lcap22p_,latg2_,syf,lonf22_,latg2p_,1,1,1)
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(syf,latlen,latstr,mype,avprsjp,wgb,lcapd,lonfd,latdef)
!$&        local(lan,lat,i,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lan,lat,i,lcapf,lonff)
#endif
   do lan = 1,latlen(mype)
     lat=latstr(mype)+lan-1
#ifdef REDUCE_GRID
     lcapf=lcapd(latdef(lat))
     lonff=lonfd(latdef(lat))
#else
     lcapf=jcap1_
     lonff=lonf_
#endif
     call sph_wave2grid(syf(1,lan),syf(1,lan),2,lcapf,lonff,latdef(lat),1)
     avprsjp(lan)=0.
     do i = 1,lonff*2
       avprsjp(lan)=avprsjp(lan)+syf(i,lan)
     enddo                                                                   
     avprsjp(lan)=avprsjp(lan)/float(lonff*2)
   enddo
!
   call mplatall(avprsjp,latg2p_,avprsj,latg2_,1)
!
! ******** no multithread or multi-node for next loop
!          otherwise non-reproduciable ****************!
!
   do lat = 1,latg2_
     sumwgt=sumwgt+wgt(lat)
     avprs=avprs+avprsj(lat)*wgt(lat)
   enddo                                                                     
   avprs=avprs/sumwgt                                                        
#ifndef NOPRINT
   if(iope) then
     write(6,*) ' sumwgt=',sumwgt,' avprs=',avprs,' avprs0=',avprs0 
   endif
#endif
!                                                                               
!  pressure correction loop                                                     
!                                                                               
   if(avprs0.gt.0) then 
#ifdef ORIGIN_THREAD
!$doacross share(qa) local(i)
#endif
#ifdef OPENMP
!$omp parallel do private(i)
#endif
     do i = 1,lln22p_
       qa(i)=0.0
     enddo
#ifdef ORIGIN_THREAD
!$doacross share(syf,avprs,avprs0,latlen,mype,lcapd,lonfd,latdef)                         
!$& local(lan,lat,i,lcapf,lonff) 
#endif
#ifdef OPENMP
!$omp parallel do private(lan,lat,i,lcapf,lonff)
#endif
     do lan = 1,latlen(mype)
       lat=latstr(mype)+lan-1
#ifdef REDUCE_GRID
       lcapf=lcapd(latdef(lat))
       lonff=lonfd(latdef(lat))
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       do i = 1,lonff*2
         syf(i,lan)=syf(i,lan)-avprs+avprs0
       enddo                                                                 
       call sph_wave2grid(syf(1,lan),syf(1,lan),2,lcapf,lonff,latdef(lat),-1)
     enddo
     call mpny2nl(syf,lonf22_,latg2p_,prs,lcap22p_,latg2_,1,1,1)
!
#ifdef ORIGIN_THREAD
!$doacross share(prs,flp,flm,qa,llens,
!$&              llstr,llens,lwvdef,lcapdp,mype)
!$& local(lat,i,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,llensd)
#endif
     do lat = 1,latg2_
#ifdef REDUCE_GRID
       llensd=lcapdp(lat,mype)
#else
       llensd=llens
#endif
       call sph_grid2wave(flp(1,1,lat),flm(1,1,lat),prs(1,lat),llensd,1)
     enddo
!
! ******** no multithread for the next loop
!          otherwise non-reproduciable ****************!
!
     do lat = 1,latg2_
#ifdef REDUCE_GRID
       llensd=lcapdp(lat,mype)
#else
       llensd=llens
#endif
       call sph_comp_coeff2(flp(1,1,lat),flm(1,1,lat),qa,qww(1,lat),           &
                  llstr,llensd,lwvdef,1)
     enddo
!
     call mpn2nn(qa,lln22p_,q,lnt22p_,1)
!
   else
     avprs0=avprs
!                                                                               
   endif
!                                                                               
#else
! ================================================================
! ============= following is non MP ==============================
! ================================================================
   llstr=0
   llens=jcap1_
!                                                                               
!  loop to find average pressure                                                
!                                                                               
   avprs=0.                                                                  
   sumwgt=0.                                                                 
#ifdef ORIGIN_THREAD
!$doacross share(qtt,colrad,q,prs,avprsj,wgt,
!$&              llstr,llens,lwvdef,lcapd,lonfd),
!$&        local(lat,i,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,lcapf,lonff)
#endif
   do lat = 1,latg2_
#ifdef REDUCE_GRID
     lcapf=lcapd(lat)
     lonff=lonfd(lat)
#else
     lcapf=jcap1_
     lonff=lonf_
#endif
#ifndef SMP
     call sph_sum_coeff(q,prs(1,lat),qtt(1,lat),llstr,lcapf,lwvdef,1)
     call sph_wave2grid(prs(1,lat),prs(1,lat),2,lcapf,lonff,lat,1)
#endif
     avprsj(lat)=0. 
     do i = 1,lonff*2
#ifdef SMP
       ij = (lat-1)*latg2_ + i
       prs(i,lat) = q(ij)
#endif
       avprsj(lat)=avprsj(lat)+prs(i,lat)
     enddo
     avprsj(lat)=avprsj(lat)/float(lonff*2)
   enddo
!
! ******** no multithread or multi-node computational
!          otherwise non-reproduciable ****************!
!
   do lat = 1,latg2_
     sumwgt=sumwgt+wgt(lat)
     avprs=avprs+avprsj(lat)*wgt(lat)
   enddo
   avprs=avprs/sumwgt
#ifndef NOPRINT
   if (iope)                                                                   &
     write(6,*) 'sumwgt=',sumwgt,' avprs=',avprs,' avprs0=',avprs0
#endif
!
!  pressure correction loop
!
   if(avprs0.gt.0) then
#ifdef ORIGIN_THREAD
!$doacross share(prs,avprs,avprs0,flp,flm,qnew,qww,
!$&              llstr,llens,lwvdef,lcapd,lonfd),
!$&        local(lat,i,lcapf,lonff)  
#endif
#ifdef OPENMP
!$omp parallel do private(lat,i,lcapf,lonff)
#endif
     do lat = 1,latg2_ 
#ifdef REDUCE_GRID
       lcapf=lcapd(lat)
       lonff=lonfd(lat)
#else
       lcapf=jcap1_
       lonff=lonf_
#endif
       do i = 1,lnt2_
         qnew(i,lat)=0.e0
       enddo
       do i = 1,lonff*2
           prs(i,lat)=prs(i,lat)-avprs+avprs0
#ifdef SMP
           qnew(i,lat) = prs(i,lat)
#endif
       enddo
#ifndef SMP
       call sph_wave2grid(prs(1,lat),prs(1,lat),2,lcapf,lonff,lat,-1)
       call sph_grid2wave(flp(1,1,lat),flm(1,1,lat),prs(1,lat),lcapf,1)
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
       call sph_grid2wave(flp(1,1,lat),flm(1,1,lat),qnew(1,lat),qww(1,lat),    &
                 llstr,lcapf,lwvdef,1)         
#endif
#ifdef DEFAULT
       call sph_comp_coeff2(flp(1,1,lat),flm(1,1,lat),qnew(1,lat),qww(1,lat),  &
                  llstr,lcapf,lwvdef,1)        
#endif
#endif				/* not SMP */
     enddo
!                                                                               
#ifdef ORIGIN_THREAD
!$doacross share(q,qnew) local(i,lat)                                           
#endif
#ifdef OPENMP
!$omp parallel do private(i,lat)
#endif
     do i = 1,lnt2_ 
       q(i)=0.0
       do lat=1,latg2_
         q(i)=q(i)+qnew(i,lat)
       enddo
     enddo
   else
     avprs0=avprs
   endif
!
!===================================================================
!
#endif
   return
   end subroutine sph_mass_adjustment
!
