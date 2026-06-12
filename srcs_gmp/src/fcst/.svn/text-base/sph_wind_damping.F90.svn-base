#include <define.h>
   subroutine sph_wind_damping(div,vor,tem,rmx,deltim,fdamp,arn,spdmax,        &
                     snnp1,lnts2,lnoffset)                
!-------------------------------------------------------------------------------
!
! program history log:
!   1988-04-25  joseph sela
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatable with namelist input
!
!-------------------------------------------------------------------------------
   use constant, only    :  rerth_
   use paramodel, only   : LNT22S,levs_,levh_,lnt2_,ntotal_
   real                 ::  fdamp(LNT22S),arn(LNT22S),spdmax(levs_)
   real                 ::  rmx(LNT22S,levh_),snnp1(lnt2_)
   real                 ::  vor(LNT22S,levs_), div(LNT22S,levs_)
   real                 ::  tem(LNT22S,levs_)
!-------------------------------------------------------------------------------
   do j = 1,lnts2                                                          
     jj=lnoffset+j
     arn(j)=snnp1(jj)+0.25e0                                                   
     arn(j)= sqrt(arn(j))                                                      
     arn(j)=arn(j)-0.5e0                                                      
   enddo
   alfa=2.5e0                                                               
   beta=rerth_*1.009e0/deltim                                               
   alfadt=alfa*deltim/rerth_                                                 
!
   do k = 1,levs_                                                           
     rncrit=beta/spdmax(k)                                                     
     coef=alfadt*spdmax(k)                                                     
     do j = 1,lnts2                                                          
       if (arn(j).gt.rncrit) then                                                
         div(j,k) =div(j,k)/(1.+(arn(j)-rncrit)*coef)                              
         vor(j,k) =vor(j,k)/(1.+(arn(j)-rncrit)*coef)                              
         tem(j,k) =tem(j,k)/(1.+(arn(j)-rncrit)*coef)                              
       end if                                                                    
     enddo
#ifndef NISLQ
      do it=1,ntotal_
        is=(it-1)*levs_                                                           
        l=is+k                                                                    
        do j = 1,lnts2                                                          
          if (arn(j).gt.rncrit) then                                                
            rmx(j,l) =rmx(j,l)/(1.+(arn(j)-rncrit)*coef)                              
          end if                                                                    
        enddo
      enddo
#endif
   enddo
!
   return                                                                    
   end subroutine sph_wind_damping  
