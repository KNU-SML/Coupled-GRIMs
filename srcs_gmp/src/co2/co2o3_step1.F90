!
   subroutine co2o3_step1(sgtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)             
!-------------------------------------------------------------------------------
   use paramodel, only : nl=>levs_,nlp=>levp1_,nlp2=>levp2_
!-------------------------------------------------------------------------------
!                                                                               
! **         this program calculates temperatures ,h2o mixing ratios            
! **         and o3 mixing ratios by using an analytical                        
! **         function which approximates                                        
! **         the us standard (1976).  this is                                   
! **         calculated in function 'co2_temp', which is called by the            
! **         main program.  the form of the analytical function was             
! **         suggested to me in 1971 by richard s. lindzen.                     
! ******************************************************************            
!         code to save stemp,gtemp on data set,bracketed by co222  **           
!             ....k. campana march 88,october 88                                
!
   real  ::  sgtemp(nlp,2),t41(nlp2,2),t42(nlp),                               & 
             t43(nlp2,2),t44(nlp)                                            
   real  ::  sglvnu(nlp),siglnu(nl)                                          
!
! *****this version is only usable for 1976 us std atm and obtains               
!     quantities for co2 interpolation and insertion into opera-                
!     tional radiation codes                                                    
!                                                                               
   character(len=20) profil                                                       
   real  ::  press(nlp),temp(nlp),alt(nlp),wmix(nlp),o3mix(nlp)              
   real  ::  wmxint(nlp,4),wmxout(nlp2),omxint(nlp,4),omxout(nlp2)           
   real  ::  pd(nlp2),gtemp(nlp)                                             
   real  ::  prs(nlp),temps(nlp),prsint(nlp),tmpint(nlp,4),a(nlp,4)          
   real  ::  prout(nlp2),tmpout(nlp2),tmpflx(nlp2),tmpmid(nlp2)              
!                                                                               
   data profil/                                                                & 
        'us standard 1976'/                                                    
   data psmax/1013.250/                                                      
!                                                                               
! **         ntype is an integer variable which has the following               
! **        values:    0 =sigma levels are used;   1= skyhi l40 levels          
! **        are used;   2 = skyhi l80 levels are used. default: 0               
!                                                                               
   ntype=0                                                                   
!o222 read (*,*) ntype                                                          
!    5 nlev=nl                                                                   
   nlev=nl                                                                   
   delzap=0.5                                                                
   r=8.31432                                                                 
   g0=9.80665                                                                
   zmass=28.9644                                                             
   aa=6356.766                                                               
   alt(1)=0.0                                                             
   temp(1)=co2_temp(6,0.0)                                                  
!
!*******determine the pressures (press)                                         
!                                                                               
   pstar=psmax                                                               
!                                                                               
!o222 if (ntype.eq.1) call skyp(pstar,pd,gtemp)                                 
!o222 if (ntype.eq.2) call sky80p(pstar,pd,gtemp)                               
!o222 if (ntype.eq.0) call co2_sigma_interp(pstar,pd,gtemp)        
!cc----      call co2_sigma_interp(pstar,pd,gtemp)                
!                                                                               
   call co2_sigma_interp(pstar,pd,gtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)
   pd(nlp2)=pstar                                                            
!
   do n = 1,nlp                                                             
     prsint(n)=pd(nlp2+1-n)                                                    
   enddo
!
!    *** calculate temps for several pressures to do quadrature                 
!
   do nq = 1,4                                                             
     do n = 2,nlp                                                            
       press(n)=prsint(n)+0.25*(nq-1)*(prsint(n-1)-prsint(n))                    
     enddo
     press(1)=prsint(1)                                                        
!
     do n = 1,nlev                                                           
!                                                                               
! **         establish computatational levels between user levels at            
! **         intervals of approximately 'delzap' km.                            
!                                                                               
       dlogp=7.0*alog(press(n)/press(n+1))                                       
       nint=dlogp/delzap                                                         
       nint=nint+1                                                               
       znint=nint                                                                
       g=g0                                                                      
       dz=r*dlogp/(7.0*zmass*g*znint)                                            
       ht=alt(n)                                                                 
!                                                                               
! **         calculate height at next user level by means of                    
! **                   runge-kutta integration.                                 
!                                                                               
       do m = 1,nint                                                           
         rk1=co2_temp(6,ht)*dz                                                       
         rk2=co2_temp(6,ht+0.5*rk1)*dz                                               
         rk3=co2_temp(6,ht+0.5*rk2)*dz                                               
         rk4=co2_temp(6,ht+rk3)*dz                                                   
         ht=ht+0.16666667*(rk1+rk2+rk2+rk3+rk3+rk4)                                
       enddo
       alt(n+1)=ht                                                               
       temp(n+1)=co2_temp(6,ht)                                                    
     enddo
     do n = 1,nlp                                                            
       tmpint(n,nq)=temp(n)                                                      
       a(n,nq)=alt(n)                                                            
     enddo
     do n = 1,nlp                                                            
       call co2_mixing_ratio(6,press(n),tmpint(n,nq),wmxint(n,nq))     
       wmxint(n,nq)=wmxint(n,nq)*1.0e-6                                          
       call co2_ozone_data(6,press(n),tmpint(n,nq),omxint(n,nq))
     enddo
   enddo
!
!    ***apply simpson s rule                                                    
!
   do n = 2,nlp                                                            
     temp(n)=1./12.*(tmpint(n-1,1)+tmpint(n,1)+4.*tmpint(n,2)+                 & 
           2.*tmpint(n,3)+4.*tmpint(n,4))                                           
     wmix(n)=1./12.*(wmxint(n-1,1)+wmxint(n,1)+4.*wmxint(n,2)+                 & 
           2.*wmxint(n,3)+4.*wmxint(n,4))                                           
     o3mix(n)=1./12.*(omxint(n-1,1)+omxint(n,1)+4.*omxint(n,2)+                & 
           2.*omxint(n,3)+4.*omxint(n,4))                                           
   enddo
!
!***output for line-by-line calcs                                               
!
   tmpout(1)=tmpint(nlp,1)                                                   
   wmxout(1)=wmxint(nlp,1)                                                   
   omxout(1)=omxint(nlp,1)                                                   
   tmpmid(1)=tmpint(nlp,1)                                                   
   do i = 2,nlp2                                                           
     tmpout(i)=temp(nlp2+1-i)                                                  
     wmxout(i)=wmix(nlp2+1-i)                                                  
     omxout(i)=o3mix(nlp2+1-i)                                                 
     tmpmid(i)=tmpint(nlp2+1-i,3)                                              
   enddo
!
   do i = 1,nlp2                                                           
     prout(i)=pd(i)                                                            
   enddo
!
   do i = 2,nlp2                                                          
     tmpflx(i)=tmpint(nlp2+1-i,1)                                              
   enddo
!
   tmpflx(1)=tmpflx(2)                                                       
!                                                                               
! **        calculate water mixing ratio using luther program                   
!                                                                               
   do n = 1,nlp                                                            
     call co2_mixing_ratio(6,prsint(n),tmpint(n,1),wmix(n))           
     wmix(n)=wmix(n)*1.0e-06                                                   
   enddo
!                                                                               
!   calculate ozone mixing ratio using schwarzkopf program                      
!                                                                               
   do n = 1,nlp                                                            
     call co2_ozone_data(6,prsint(n),tmpint(n,1),o3mix(n))
   enddo
!                                                                               
   write (6,101) profil                                                      
   101 format (1x,a20)                                                           
   write (6,201)                                                             
   201 format(5x,'   height    temperature     pressure     r(h2o)     ',      & 
            ' r(o3)')                                                                 
   write (6,202)  a(1,1),tmpint(1,1),prsint(1),wmix(1),o3mix(1)              
   do n = 2,nlp                                                            
     write (6,203)  tmpout(nlp2+1-n),wmxout(nlp2+1-n),omxout(nlp2+1-n)         
     write (6,202)  a(n,1),tmpint(n,1),prsint(n),wmix(n),o3mix(n)              
   enddo
!
   203 format (1x,14x,f14.6,14x,2e14.6)                                          
   202 format(1x,2f14.6,e14.6,2e14.6)                                            
!
!   determine what formats to use***                                            
!
   nrep=nlp/5                                                                
   nrem=nlp-5*nrep                                                           
   if (nrem.eq.0) then                                                       
     nrep=nrep-1                                                             
     nrem=5                                                                  
   endif                                                                     
   nrep2=nl/5                                                                
   nrem2=nl-5*nrep2                                                          
   if (nrem2.eq.0) then                                                      
     nrep2=nrep2-1                                                           
     nrem2=5                                                                 
   endif                                                                     
   nrep3=nlp/4                                                               
   nrem3=nlp-4*nrep3                                                         
   if (nrem3.eq.0) then                                                      
     nrep3=nrep3-1                                                           
     nrem3=4                                                                 
   endif                                                                     
!
!o222   *****************************************************                   
!c    rewind 66                                                                 
!o222   *****************************************************                   
!***output temperatures                                                         
!
   do iout = 1,2                                                           
     if (iout.eq.1)  then                                                   
       write (16,701)                                                       
     else                                                                   
       write (16,702)                                                       
     endif                                                                  
     701   format (6x,'data dtemp /')                                                
     702   format (6x,'data stemp /')                                                
     nf=0                                                                      
!
!o222   *****************************************************                   
!         save stemp                                                            
!c      if(iout.eq.2) write(66) (tmpint(nlp2-n,1),n=1,nlp)                      
     do n = 1,nlp                                                            
       sgtemp(n,1) = tmpint(nlp2-n,1)                                          
     enddo
!
!o222   *****************************************************                   
!
     if (nrep.ne.0) then                                                       
       do nr = 1,nrep                                                       
         ns=5*(nr-1)+1                                                        
         nf=ns+4                                                              
         write (16,656) (tmpint(nlp2-n,1),n=ns,nf)                            
       enddo
     endif                                                                     
     nf=nf+1                                                                   
     if (nrem.eq.1) then                                                       
       write (16,616) (tmpint(nlp2-n,1),n=nf,nlp)                            
     endif                                                                     
     if (nrem.eq.2) then                                                       
       write (16,626) (tmpint(nlp2-n,1),n=nf,nlp)                            
     endif                                                                     
     if (nrem.eq.3) then                                                       
       write (16,636) (tmpint(nlp2-n,1),n=nf,nlp)                            
     endif                                                                     
     if (nrem.eq.4) then                                                       
       write (16,646) (tmpint(nlp2-n,1),n=nf,nlp)                            
     endif                                                                     
     if (nrem.eq.5) then                                                       
       write (16,606) (tmpint(nlp2-n,1),n=nf,nlp)                            
     endif                                                                     
     616   format (5x,'*',1x,f12.6,'/')                                              
       626   format (5x,'*',1x,f12.6,',',f12.6,'/')                                    
     636   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,'/')                          
     646   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                    & 
                  f12.6,'/')                                              
     656   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                    & 
                  f12.6,',',f12.6,',')                                    
     606   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                    & 
                  f12.6,',',f12.6,'/')                                    
   enddo
!
!***output gtemp                                                                
!o222   *****************************************************                   
!         save gtemp                                                            
!c             write(66) (gtemp(n),n=1,nlp)                                     
!
   do n = 1,nlp                                                            
     sgtemp(n,2) = gtemp(n)                                                  
   enddo
!
!o222   *****************************************************                   
!
   write (16,706)                                                            
   706   format (6x,'data gtemp /')                                                
   nf=0                                                                      
!
   if (nrep3.ne.0) then                                                      
     do nr = 1,nrep3                                                      
       ns=4*(nr-1)+1                                                        
       nf=ns+3                                                              
       write (16,648) (gtemp(n),n=ns,nf)                                    
     enddo
   endif                                                                     
!
   nf=nf+1                                                                   
   if (nrem3.eq.1) then                                                      
     write (16,618) (gtemp(n),n=nf,nlp)                                    
   endif                                                                     
!
   if (nrem3.eq.2) then                                                      
     write (16,628) (gtemp(n),n=nf,nlp)                                    
   endif                                                                     
!
   if (nrem3.eq.3) then                                                      
     write (16,638) (gtemp(n),n=nf,nlp)                                    
   endif                                                                     
!
   if (nrem3.eq.4) then                                                      
     write (16,608) (gtemp(n),n=nf,nlp)                                    
   endif                                                                     
!
!***output wmix in 10**2 gm/gm                                                  
!
   do n = 2,nlp                                                            
     wmix(n)=wmix(n)*1.0e2                                                   
   enddo
!
   nf=0                                                                      
   write (16,703)                                                            
   703  format (6x,'data rr /')                                                   
!
   if (nrep2.ne.0) then                                                      
     do nr = 1,nrep2                                                      
       ns=5*(nr-1)+1                                                        
       nf=ns+4                                                              
       write (16,657) (wmix(nlp2-n),n=ns,nf)                                
     enddo
   endif                                                                     
!
   nf=nf+1                                                                   
   if (nrem2.eq.1) then                                                      
     write (16,617) (wmix(nlp2-n),n=nf,nl)                                 
   endif                                                                     
!
   if (nrem2.eq.2) then                                                      
     write (16,627) (wmix(nlp2-n),n=nf,nl)                                 
   endif                                                                     
!
   if (nrem2.eq.3) then                                                      
     write (16,637) (wmix(nlp2-n),n=nf,nl)                                 
   endif                                                                     
!
   if (nrem2.eq.4) then                                                      
     write (16,647) (wmix(nlp2-n),n=nf,nl)                                 
   endif                                                                     
!
   if (nrem2.eq.5) then                                                      
     write (16,607) (wmix(nlp2-n),n=nf,nl)                                 
   endif                                                                     
!
   617   format (5x,'*',1x,e12.6,'/')                                              
   627   format (5x,'*',1x,e12.6,',',e12.6,'/')                                    
   637   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,'/')                          
   647   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                      &
               e12.6,'/')                                              
   657   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                      &
               e12.6,',',e12.6,',')                                    
   607   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                      &
               e12.6,',',e12.6,'/')                                    
!
!***output o3mix in 10**3 gm/gm                                                 
!
   do n = 2,nlp                                                            
     o3mix(n)=o3mix(n)*1.0e3                                                
   enddo
!
   nf=0                                                                      
   write (16,704)                                                            
   704   format (6x,'data qqo3 /')                                                 
!
   if (nrep2.ne.0) then                                                      
     do nr = 1,nrep2                                                      
       ns=5*(nr-1)+1                                                        
       nf=ns+4                                                              
       write (16,657) (o3mix(nlp2-n),n=ns,nf)                               
     enddo
   endif                                                                     
!
   nf=nf+1                                                                   
   if (nrem2.eq.1) then                                                      
     write (16,617) (o3mix(nlp2-n),n=nf,nl)                                
   endif                                                                     
!
   if (nrem2.eq.2) then                                                      
     write (16,627) (o3mix(nlp2-n),n=nf,nl)                                
   endif                                                                     
!
   if (nrem2.eq.3) then                                                      
     write (16,637) (o3mix(nlp2-n),n=nf,nl)                                
   endif                                                                     
!
   if (nrem2.eq.4) then                                                      
     write (16,647) (o3mix(nlp2-n),n=nf,nl)                                
   endif                                                                     
!
   if (nrem2.eq.5) then                                                      
     write (16,607) (o3mix(nlp2-n),n=nf,nl)                                
   endif                                                                     
!
!***output prout in cgs                                                         
!
   do n = 2,nlp2                                                           
     prout(n)=prout(n)*1.0e3                                                 
   enddo
!
   write (16,705)                                                            
   705   format (6x,'data ppress /')                                               
   nf=0                                                                      
!
   if (nrep3.ne.0) then                                                      
     do nr = 1,nrep3                                                      
       ns=4*(nr-1)+2                                                        
       nf=ns+3                                                              
       write (16,648) (prout(n),n=ns,nf)                                    
     enddo
   endif                                                                     
!
   nf=nf+1                                                                   
   if (nrem3.eq.1) then                                                      
    write (16,618) (prout(n),n=nf,nlp2)                                   
   endif                                                                     
!
   if (nrem3.eq.2) then                                                      
    write (16,628) (prout(n),n=nf,nlp2)                                   
   endif                                                                     
!
   if (nrem3.eq.3) then                                                      
    write (16,638) (prout(n),n=nf,nlp2)                                   
   endif                                                                     
!
   if (nrem3.eq.4) then                                                      
    write (16,608) (prout(n),n=nf,nlp2)                                   
   endif                                                                     
!
   618   format (5x,'*',1x,e15.9,'/')                                              
   628   format (5x,'*',1x,e15.9,',',e15.9,'/')                                    
   638   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,'/')                          
   648   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,',',e15.9,',')                
   608   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,',',e15.9,'/')                
!
   return                                                                    
   end subroutine co2o3_step1
!-------------------------------------------------------------------------------
