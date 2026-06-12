!
   subroutine znlsfc(slimsk,idim,jdim,ijdim,lsoil,nrcznl,                      &
                     tsfanl,wetanl,snoanl,stcanl,tg3anl,                       &
                     zoranl,plranl,cvanl,cvbanl,cvtanl,                        &
                     albanl,aisanl,smcanl,cnpanl,                              &
                     znlsl,nzl1,nzl2,lsmsk,weis,gaul)                        
!-------------------------------------------------------------------------------
   use paramodel, only : latg_,lonf_
!-------------------------------------------------------------------------------
   real   ::  tsfanl(idim,jdim),wetanl(idim,jdim),snoanl(idim,jdim),           &
              stcanl(idim,jdim,lsoil),tg3anl(idim,jdim),                       &
              zoranl(idim,jdim),plranl(idim,jdim),                             &
              cvanl (idim,jdim),cvbanl(idim,jdim),cvtanl(idim,jdim),           &
              albanl(idim,jdim),aisanl(idim,jdim),                             &
              smcanl(idim,jdim,lsoil),cnpanl(idim,jdim)                        
!                                                                               
   logical lsmsk(idim,2,6)                                                   
   real   ::   znlsl(6,6,nrcznl),weis(6,6)                                     
!                                                                               
   real   ::   gaul(jdim)                                                      
!                                                                               
   real   ::   wgt(latg_)                                                      
   real   ::   work1(lonf_,2), work2(lonf_,2)                                  
!                                                                               
   real   ::   slimsk(idim,jdim)                                               
!                                                                               
   idimt=idim*2                                                              
   jdimhf=jdim/2                                                             
!                                                                               
   do j = 1,jdim                                                               
     wgt(j)=cos((90.-gaul(j))*0.01745329)                                    
   enddo                                                                     
!                                                                               
!    lat loop                                                                  
!                                                                               
   do lat = 1,jdimhf                                                    
     latco = jdim + 1 - lat                                                    
!                                                                               
!   zonal average monitoring                                                    
!                                                                               
     do i = 1,idim                                                         
       work1(i,1) = slimsk(i,lat  )                                              
       work1(i,2) = slimsk(i,latco)                                              
       work2(i,1) = snoanl(i,lat )                                               
       work2(i,2) = snoanl(i,latco)                                              
     enddo
     call znlwgt(work1,work2,wgt(lat),lat,                                     &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
!                                                                               
     do i = 1,idim                                                         
       work1(i,1) = tsfanl(i,lat  )                                              
       work1(i,2) = tsfanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat), 1,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = smcanl(i,lat,1)                                              
       work1(i,2) = smcanl(i,latco,1)                                            
     enddo
     call znlavs(work1  ,lat,wgt(lat), 2,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = smcanl(i,lat,lsoil)                                          
       work1(i,2) = smcanl(i,latco,lsoil)                                        
     enddo
     call znlavs(work1  ,lat,wgt(lat), 3,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = snoanl(i,lat       )                                         
       work1(i,2) = snoanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat), 4,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = stcanl(i,lat       ,1)                                       
       work1(i,2) = stcanl(i,latco,1)                                            
     enddo
     call znlavs(work1  ,lat,wgt(lat), 5,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = stcanl(i,lat       ,lsoil)                                   
       work1(i,2) = stcanl(i,latco,lsoil)                                        
     enddo
     call znlavs(work1  ,lat,wgt(lat), 6,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = zoranl(i,lat       )                                         
       work1(i,2) = zoranl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat), 7,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) =   cvanl(i,lat  )                                             
       work1(i,2) =   cvanl(i,latco)                                             
     enddo
     call znlavs(work1  ,lat,wgt(lat), 8,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                         
       work1(i,1) = cvbanl(i,lat  )                                              
       work1(i,2) = cvbanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat), 9,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                        
       work1(i,1) = cvtanl(i,lat  )                                              
       work1(i,2) = cvtanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat),10,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                        
       work1(i,1) = albanl(i,lat  )                                              
       work1(i,2) = albanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat),11,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                        
       work1(i,1) = plranl(i,lat  )                                              
       work1(i,2) = plranl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat),12,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
     do i = 1,idim                                                        
       work1(i,1) = cnpanl(i,lat  )                                              
       work1(i,2) = cnpanl(i,latco)                                              
     enddo
     call znlavs(work1  ,lat,wgt(lat),13,                                      &
              znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)            
!                                                                               
   enddo
!                                                                               
   return                                                                    
   end subroutine znlsfc                                                                       
!-------------------------------------------------------------------------------
