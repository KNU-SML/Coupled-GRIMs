#include <define.h>
   subroutine phys_scv_mrf(klev,dt,delprsi,prsi,prsik,prsl,prslk,zl,icps,q,t,  &
                       cp,g,hvap,rd,rv,                                        &
                       ids,ide, jds,jde, kds,kde,                              &
                       ims,ime, jms,jme, kms,kme,                              &
                       its,ite, jts,jte, kts,kte)
!-------------------------------------------------------------------------------
! 
! subroutine:   phys_scv_mrf
!
! abstract: 
! - computes shallow convective heating and moisture
! - sub-grid-scale shallow convective cloud parameterization.
!   this routine computes the effects of shallow convection
!   based on tiedtke (1984), ecmwf workshop on convection in
!   large-scale numerical models.
!
! program history log:
!   1991-03-19  peter caplan
!   1991-03-19  hua-lu pan
!   1991-05-07  iredell     arguments changed, scv_tri_diagonal split off
!   1999-03-01  hong        add profile function with newshal option
!
! usage:    call phys_scv_mrf(im,km,dt,del,si,sl,slk,icps,ps,q,t)
!
!   input argument list:
!     im       - integer number of points
!     km       - integer number of levels
!     dt       - real time step in seconds
!     del      - real (km) sigma layer thickness
!     sl       - real (km) sigma values
!     slk      - real (km) sigma values to the kappa
!     ps       - real (im) surface pressure in kilopascals (cb)
!     q        - real (im,km) current specific humidity in kg/kg
!     t        - real (im,km) current temperature in kelvin
!
!   output argument list:
!     q        - real (im,km) adjusted specific humidity in kg/kg
!     t        - real (im,km) adjusted temperature in kelvin
!
! subprograms called:
!   phys_moist_adiabat  - computes moist adiabat and returns cloud values
!   scv_tri_diagonal    - solves tridiagonal matrix problem
!
! remarks: nonstandard automatic arrays are used.
!
!  ::: structure :::
!
!    [phys_scv_mrf] --- [scv_tri_diagonal]
! 
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  klev,                                              &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
   integer              ::  icps(ims:ime)
   real                 ::  g,rd,rv,cp,hvap,akapa
   real                 ::  delprsi(ims:ime,kms:kme)
   real                 ::  prsi(ims:ime,kms:kme+1),prsik(ims:ime,kms:kme+1)
   real                 ::  prsl(ims:ime,kms:kme),prslk(ims:ime,kms:kme)
   real                 ::  zl(ims:ime,kms:kme)
   real                 ::  q(ims:ime,kms:kme),t(ims:ime,kms:kme)
!
!  bounds of parcel origin
!
   integer,parameter    ::  kliftl=2,kliftu=2
!
!  local variables and arrays
!
   logical              ::  lshc(its:ite)
   integer              ::  i,ik,ik1,iku,k,k1,k2,kt,n2
   integer              ::  index2(its:ite)
   integer              ::  klcl(its:ite),kbot(its:ite)
   integer              ::  ktop(its:ite)
   real                 ::  eps,epsm1
   real                 ::  eldq,ck,cpdt,rtdls
   real                 ::  dmse,dtodsu,dtodsl,dt,dsig,dsdz1,dsdz2 
   real                 ::  q2(its:ite*kte)
   real                 ::  t2(its:ite*kte)
   real                 ::  al(its:ite*(kte-1))
   real                 ::  ad(its:ite*kte)
   real                 ::  au(its:ite*(kte-1))
   real                 ::  delprsi2(its:ite*kte)
   real                 ::  prsi2(its:ite*kte),prsik2(its:ite*kte)
   real                 ::  prsl2(its:ite*kte),prslk2(its:ite*kte)
   real                 ::  prslimsi(its:ite,kts:kte)
!-------------------------------------------------------------------------------
!
! constants
!
!
!  compress fields to points with no deep convection
!  and moist static instability.
!
   do i = its,ite
     lshc(i)=.false.
   enddo
!
   do k = kts,klev-1
     do i = its,ite
       if(icps(i).eq.0) then
         eldq=hvap*(q(i,k)-q(i,k+1))
         cpdt=cp*(t(i,k)-t(i,k+1))
         rtdls=(prsl(i,k)-prsl(i,k+1))/prsi(i,k+1)*rd*0.5*(t(i,k)+t(i,k+1))
         dmse=eldq+cpdt-rtdls
         lshc(i)=lshc(i).or.dmse.gt.0.
       endif
     enddo
   enddo
!
   n2=0
   do i = its,ite
     if(lshc(i)) then
       n2=n2+1
       index2(n2)=i
     endif
   enddo
#if defined(DBG) && !defined(MP)
       print * , 'number of SCV = ' ,n2 
#endif
!
   if(n2.eq.0) return
!
   do k = kts,klev
     do i = 1,n2
       ik=(k-1)*n2+i
       q2(ik)=q(index2(i),k)
       t2(ik)=t(index2(i),k)
       delprsi2(ik)=delprsi(index2(i),k)
       prsi2(ik)=prsi(index2(i),k)
       prsl2(ik)=prsl(index2(i),k)
       prsik2(ik)=prsik(index2(i),k)
       prslk2(ik)=prslk(index2(i),k)
     enddo
   enddo
!
!  compute moist adiabat and determine limits of shallow convection.
!  check for moist static instability again within cloud.
!
   call phys_moist_adiabat(n2,klev-1,kliftl,kliftu,                             &
                           prsl2,prsik2,prslk2,t2,q2,                           &
                           klcl,kbot,ktop,al,au,rd,rv,                          &
                            ids,ide, jds,jde, kds,kde,                          &
                            ims,ime, jms,jme, kms,kme,                          &
                            its,ite, jts,jte, kts,kte)
!
   do i = 1,n2
     kbot(i)=klcl(i)-1
     ktop(i)=ktop(i)+1
     lshc(i)=.false.
   enddo
!
   do k = kts,klev-1
     do i = 1,n2
       if(k.ge.kbot(i).and.k.lt.ktop(i)) then
         ik=(k-1)*n2+i
         iku=k*n2+i
         eldq=hvap*(q2(ik)-q2(iku))
         cpdt=cp*(t2(ik)-t2(iku))
         rtdls=(prsl2(ik)-prsl2(iku))/prsi2(iku)*rd*0.5*(t2(ik)+t2(iku))
         dmse=eldq+cpdt-rtdls
         lshc(i)=lshc(i).or.dmse.gt.0.
         au(ik)=g/rtdls
       endif
     enddo
   enddo
!
   k1=klev+1
   k2=0
   do i = 1,n2
     if(.not.lshc(i)) then
       kbot(i)=klev+1
       ktop(i)=0
     endif
     k1=min(k1,kbot(i))
     k2=max(k2,ktop(i))
   enddo
   kt=k2-k1+1
   if(kt.lt.2) return
!
!  set eddy viscosity coefficient cku at sigma interfaces.
!  compute diagonals and rhs for tridiagonal matrix solver.
!  expand final fields.
!
   do i = 1,n2
     ik=(k1-1)*n2+i
     ad(ik)=1.
   enddo
!
   do k = k1,k2-1
     do i = 1,n2
       ik=(k-1)*n2+i
       iku=k*n2+i
       dtodsl=2.*dt/delprsi2(ik)*prsi2(i)
       dtodsu=2.*dt/delprsi2(iku)*prsi2(i)
       dsig=(prsl2(ik)-prsl2(iku))/prsi2(i)
       if(k.eq.kbot(i)) then
         ck=1.5
       elseif(k.eq.ktop(i)-1) then
         ck=1.
       elseif(k.eq.ktop(i)-2) then
         ck=3.
       elseif(k.gt.kbot(i).and.k.lt.ktop(i)-2) then
         ck=5.
       else
         ck=0.
       endif
       dsdz1=ck*dsig*au(ik)*g/cp
       dsdz2=ck*dsig*au(ik)*au(ik)
       au(ik)=-dtodsl*dsdz2
       al(ik)=-dtodsu*dsdz2
       ad(ik)=ad(ik)-au(ik)
       ad(iku)=1.-al(ik)
       t2(ik)=t2(ik)+dtodsl*dsdz1
       t2(iku)=t2(iku)-dtodsu*dsdz1
     enddo
   enddo
!
   ik1=(k1-1)*n2+1
   call scv_tri_diagonal(n2,n2,kt,al(ik1),ad(ik1),au(ik1),q2(ik1),t2(ik1),     &
                                                  au(ik1),q2(ik1),t2(ik1))
   do k = k1,k2
     do i = 1,n2
       ik=(k-1)*n2+i
       q(index2(i),k)=q2(ik)
       t(index2(i),k)=t2(ik)
     enddo
   enddo
!
   return
   end subroutine phys_scv_mrf
!-----------------------------------------------------------------------
   subroutine scv_tri_diagonal(lons2,l,n,cl,cm,cu,r1,r2,au,a1,a2)
!-----------------------------------------------------------------------        
!
! subprogram:  scv_tri_diagonal    
!                                                                               
! abstract: this routine solves multiple tridiagonal matrix problems            
!   with 2 right-hand-side and solution vectors for every matrix.               
!   the solutions are found by eliminating off-diagonal coefficients,           
!   marching first foreward then backward along the matrix diagonal.            
!   the computations are vectorized around the number of matrices.              
!   no checks are made for zeroes on the diagonal or singularity.               
!                                                                               
! program history log:                                                          
!   1991-05-07  iredell                                                           
!                                                                               
! usage:    call scv_tri_diagonal(l,n,cl,cm,cu,r1,r2,au,a1,a2)      
!                                                                               
!   input argument list:                                                        
!     l        - integer number of tridiagonal matrices                         
!     n        - integer order of the matrices                                  
!     cl       - real (l,2:n) lower diagonal matrix elements                    
!     cm       - real (l,n) main diagonal matrix elements                       
!     cu       - real (l,n-1) upper diagonal matrix elements                    
!                (may be equivalent to au if no longer needed)                  
!     r1       - real (l,n) 1st right-hand-side vector elements                 
!                (may be equivalent to a1 if no longer needed)                  
!     r2       - real (l,n) 2nd right-hand-side vector elements                 
!                (may be equivalent to a2 if no longer needed)                  
!                                                                               
!   output argument list:                                                       
!     au       - real (l,n-1) work array                                        
!     a1       - real (l,n) 1st solution vector elements                        
!     a2       - real (l,n) 2nd solution vector elements                        
!                                                                               
! remarks: this routine can be easily modified to solve a different             
!   number of right-hand-sides and solutions per matrix besides 2.              
!                                                                               
!-------------------------------------------------------------------------------
   real                 ::  cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n),       &
                            au(l,n-1),a1(l,n),a2(l,n)                                       
!-------------------------------------------------------------------------------
   do i = 1,lons2
     fk=1./cm(i,1)                                                           
     au(i,1)=fk*cu(i,1)                                                      
     a1(i,1)=fk*r1(i,1)                                                      
     a2(i,1)=fk*r2(i,1)                                                      
   enddo                                                                     
   do k = 2,n-1                                                                
     do i = 1,lons2
       fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))                                     
       au(i,k)=fk*cu(i,k)                                                    
       a1(i,k)=fk*(r1(i,k)-cl(i,k)*a1(i,k-1))                                
       a2(i,k)=fk*(r2(i,k)-cl(i,k)*a2(i,k-1))                                
     enddo                                                                   
   enddo                                                                     
   do i = 1,lons2
     fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))                                       
     a1(i,n)=fk*(r1(i,n)-cl(i,n)*a1(i,n-1))                                  
     a2(i,n)=fk*(r2(i,n)-cl(i,n)*a2(i,n-1))                                  
   enddo                                                                     
   do k = n-1,1,-1                                                             
     do i=1,lons2
       a1(i,k)=a1(i,k)-au(i,k)*a1(i,k+1)                                     
       a2(i,k)=a2(i,k)-au(i,k)*a2(i,k+1)                                     
     enddo                                                                   
   enddo                                                                     
!-----------------------------------------------------------------------        
   return                                                                    
   end subroutine scv_tri_diagonal
