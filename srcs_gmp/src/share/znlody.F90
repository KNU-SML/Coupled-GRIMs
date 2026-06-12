!
   subroutine znlody(fh,                                                       &
                     znlsl,nzl1,nzl2,lsmsk,weis,idim,jdim,idimt,nrcznl)        
!-------------------------------------------------------------------------------
   logical lsmsk(idimt,6)                                                    
   real                 ::   znlsl(6,6,nrcznl),weis(6,6)                                     
!                                                                               
!  this routine prints out zonal averages                                       
!                                                                               
   character(len=21)    ::   ifmt0                                                        
   character(len=31)    ::   ifmt1                                                        
   character(len=37)    ::   ifmt2                                                        
   character(len=7)     ::   nlat(6)                                                      
!                                                                               
   character(len=8)     ::   lttlsl(100)                                                  
   character(len=8)     ::   lblnk                                                        
   real                 ::   ipwrsl(30)                                                     
!                                                                               
   data lblnk/'        '/                                                    
!                                                                               
   data nlat/'90n-90s','90n-60n','60n-30n','30n-30s','30s-60s',                &
             '60s-90s'/                                                      
!                                                                               
   lttlsl( 1)='tsfc    '                                                     
   lttlsl( 2)='soilm1  '                                                     
   lttlsl( 3)='soilm2  '                                                     
   lttlsl( 4)='snow    '                                                     
   lttlsl( 5)='tg1     '                                                     
   lttlsl( 6)='tg2     '                                                     
   lttlsl( 7)='z0cm    '                                                     
   lttlsl( 8)='cv      '                                                     
   lttlsl( 9)='cvb     '                                                     
   lttlsl(10)='cvt     '                                                     
   lttlsl(11)='alb     '                                                     
   lttlsl(12)='plantr  '                                                     
   lttlsl(13)='cnpwc   '                                                     
   lttlsl(14)='slimsk  '                                                     
   data ipwrsl/0,0,0,1,0,0,0,0,0,0,-2,0,0,0,16*0/                            
!                                                                               
   jdimhf=jdim/2                                                             
!                                                                               
   print *,'@@@@ start of zonal diagnostic print @@@'                        
!                                                                               
!  single level field                                                           
!                                                                               
   k1=1                                                                      
   k2=6                                                                      
   kxxx=k2-k1+1                                                              
!
   do 300 itm = 1,nrcznl                                                       
     if(lttlsl(itm).eq.lblnk) go to 300                                        
!
     write(6,100) lttlsl(itm),ipwrsl(itm),fh                                   
     write(6,110)                                                              
     kyyy=-ipwrsl(itm)                                                         
     write(ifmt0,130) kyyy,kxxx                                                
     write(6,ifmt0)(nlat(j),(znlsl(j,k,itm),k=k1,k2),j=1,6)                    
   300 continue                                                                  
!                                                                               
   print *,'@@@@ end of zonal diagnostic print @@@'                          
!                                                                               
   100 format(1x,a8,5h(10**,i3,1h),' fh=',f7.1)                                  
   110 format(2x,5h lat ,7x,'mean',6x,'lnd',3x,'sn-lnd',6x,'ice',              &
                    3x,'sn-ice',6x,'sea')                                   
   130 format(10h(1x,a7,1x,,i2,1hp,i2,5hf9.2))                                   
!
   return                                                                    
   end subroutine znlody                                                                       
!-------------------------------------------------------------------------------
