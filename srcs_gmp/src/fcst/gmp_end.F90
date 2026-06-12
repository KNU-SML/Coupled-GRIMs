#include "define.h"
   subroutine gmp_end
!-------------------------------------------------------------------------------
!
! subroutine: gmp_end         
!
! abstract: finish global forecast 
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
! references : 
!   hong et al. (2013, apjas): global/regional integrated model system (grims)
!   park et al. (2013, mwr): dfs dynamical core
!   byun and hong (2007, j. climate): single-column model (SMP)
!   kanamitsu et al. (2002, bams): ncep dynamical seasonal forecast system 2000
!
! input files:                                                                  
!   unit   11    sigma file (analysis or at time t-dt)                          
!   unit   12    sigma file (at time t if not analysis)                         
!   unit   14    surface file                                                   
!   unit   15    co2 constants (dependent on vertical resolution)               
!   unit   24    mountain variance (dependent on horizontal resolution)         
!   unit   43    cloud tuning                                                   
!                                                                               
! output files:                                                                 
!   unit   51    sigma file (at time t-dt)                                      
!   unit   52    sigma file (at time t)                                         
!   unit   53    surface file                                                   
!   unit   61    initial zonal diagnostics                                      
!   unit   63    flux diagnostics                                               
!   unit   64    final zonal diagnostics                                        
!   unit   67    grid point diagnostics                                         
!
!-------------------------------------------------------------------------------
   use paramodel, only : lnt22_,levs_
#ifdef DFS
   use dfsvar, only : psl1,t1,d1,v1,q1,sf1,tavexy,psl2,t2,d2,v2,q2            ,&
                      sl=>sigma,si=>sigmafull
#else
#ifdef MP
   use commpi
#endif
#endif /* DFS end */
   use comfibm
   use comcon
   use comgpd
   use comfgsm
#ifdef LFM
   use comlfm
#endif
#ifdef SMP
   use comfcst, only : dtbdy,curtime,vvel,hour1,hour2,hours,houre,ftim1,ftim2
#ifdef CLM_CWF
   use comfcst, only : wdiv,hadq
#endif
#endif /* SMP end */
!
#ifndef HYBRID
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
#endif
!
! write sigitdt and sigit at the end of inchour
!
!#if defined(SMP)
!   call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm ,sl,si,gz,z00,vvel,1)
!   call file_write_sigma(n1,thour,idate,q ,te ,di ,ze ,rq ,sl,si,gz,z00,vvel,2)
!#elif defined(DFS)
!#ifndef HYBRID
!   call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,sl,si,sf1,tavexy,1)
!   call file_write_sigma(n1,thour,idate,psl2,t2,d2,v2,q2,sl,si,sf1,tavexy,2)
!#else
!   call file_write_sigma(n1,thour,idate,psl1,t1,d1,v1,q1,ak5,bk5,sf1,tavexy,1)
!   call file_write_sigma(n1,thour,idate,psl2,t2,d2,v2,q2,ak5,bk5,sf1,tavexy,2)
!#endif
!#else /* SPH */
!#ifndef HYBRID
!   call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm ,sl,si,gz,z00,1)
!   call file_write_sigma(n1,thour,idate,q ,te ,di ,ze ,rq ,sl,si,gz,z00,2)
!#else
!   call file_write_sigma(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,1)
!   call file_write_sigma(n1,thour,idate,q ,te ,di ,ze ,rq,ak5,bk5,gz,z00,2)
!#endif
!#endif
!
#ifdef CHEM
   call end_chem
#endif
!
   return
   end subroutine gmp_end
!-------------------------------------------------------------------------------
