!
   subroutine co2o3_step2(itape,t15a,t15b,t22,ratio,ir,nmethd)                    
!-------------------------------------------------------------------------------
   use paramodel, only : nlevls=>levs_, nlp1=>levp1_, nlp2=>levp2_
   use comco2
!-------------------------------------------------------------------------------
!nov89                                                                          
!     *********************************************************                 
!       changes to data read  and format see co222     ***                      
!          ..... k.campana   march 1988,october 1988                            
!       changes to pass itape,and if ir=4,read 1 co2 rec..kac nov89             
!     *********************************************************                 
!       co2o3_step2 interpolates carbon dioxide transmission functions               
!  from the 109 level grid,for which the transmission functions                 
!  have been pre-calculated, to the grid structure specified by the             
!  user.                                                                        
!                                                                               
!        method:                                                                
!                                                                               
!      co2o3_step2 is employable for two purposes: 1) to obtain transmis-            
!  sivities between any 2 of an array of user-defined pressures; and            
!  2) to obtain layer-mean transmissivities between any 2 of an array           
!  of user-defined pressure layers.to clarify these two purposes,see            
!  the diagram and discussion below.                                            
!      co2o3_step2 may be used to execute only one purpose at one time.              
!                                                                               
!     let p be an array of user-defined pressures                               
!     and pd be user-defined pressure layers.                                   
!                                                                               
!       - - - - - - - - -   pd(i-1) ---                                         
!                                     ^                                         
!       -----------------   p(i)      ^  pressure layer i  (plm(i))             
!                                     ^                                         
!       - - - - - - - - -   pd(i)  ---                                          
!                                     ^                                         
!       -----------------   p(i+1)    ^  pressure layer i+1 (plm(i+1))          
!                                     ^                                         
!       - - - - - - - - -   pd(i+1)---                                          
!            ...                          (the notation used is                 
!            ...                          consistent with the code)             
!            ...                                                                
!      - - - - - - - - -    pd(j-1)                                             
!                                                                               
!      -----------------    p(j)                                                
!                                                                               
!      - - - - - - - - -    pd(j)                                               
!                                                                               
!      purpose 1:   the transmissivity between specific pressures               
!      p(i) and p(j) ,tau(p(i),p(j))  is computed by this program.              
!      in this mode,there is no reference to layer pressures pd                 
!      (pd,plm are not inputted).                                               
!                                                                               
!      purpose 2:   the layer-mean transmissivity between a layer-              
!      mean pressure plm(j) and pressure layer i is given by                    
!         taulm(plm(i),plm(j)). it is computed by the integral                  
!                                                                               
!                           pd(i)                                               
!                           ----                                                
!             1             ^                                                   
!        -------------  *   ^   tau ( p',plm(j) )  dp'                          
!        pd(i)-pd(i-1)      ^                                                   
!                        ----                                                   
!                        pd(i-1)                                                
!                                                                               
!           the layer-mean pressure plm(i) is specified by the user.            
!        for many purposes,plm will be chosen to be the average                 
!        pressure in the layer-ie,plm(i)=0.5*(pd(i-1)+pd(i)).                   
!           for layer-mean transmissivities,the user thus inputs                
!        a pressure array (pd) defining the pressure layers and an              
!        array (plm) defining the layer-mean pressures.the calculation          
!        does not depend on the p array used for purpose 1 (p is not            
!        inputted).                                                             
!                                                                               
!            the following paragraphs depict the utilization of this            
!       code when used to compute transmissivities between specific             
!       pressures. later paragraphs describe additional features needed         
!       for layer-mean transmissivities.                                        
!                                                                               
!          for a given co2 mixing ratio and standard temperature                
!      profile,a table of transmission functions for a fixed grid               
!     of atmospheric pressures has been pre-calculated.                         
!      the standard temperature profile is computed from the us                 
!     standard atmosphere (1977) table.additionally, the                        
!     same transmission functions have been pre-calculated for a                
!     temperature profile increased and decreased (at all levels)               
!     by 25 degrees.                                                            
!         this program reads in the prespecified transmission functions         
!     and a user-supplied pressure grid (p(i)) and calculates trans-            
!     mission functions ,tau(p(i),p(j)), for all p(i) s and p(j) s.             
!     a logarithmic interpolation scheme is used.                               
!         this method is repeated for the three temperature profiles            
!     given above .therefore outputs from the program are three tables          
!     of transmission functions for the user-supplied pressure grid.            
!     the existence of the three tables permits subsequent interpo-             
!     lation to a user-supplied temperature profile using the method            
!     described in the reference.see limitations section if the                 
!     user desires to obtain only 1 table of transmissivities.                  
!                                                                               
!     modifications for layer-mean transmissivities:                            
!          the pressures inputted are the layer-mean pressures,pd,              
!     and the layer-mean pressures ,plm. a series of transmissivities           
!     (tau(p'',plm(j)) are computed and the integral given in the               
!     discussion of purpose 2 is computed.for plm(i) not equal to               
!     plm(j) simpson s rule is used with 5 points. if plm(i)=plm(j)             
!     (the -nearby layer- case) a 49-point quadrature is used for               
!     greater accuracy.the output is in taulm(plm(i),plm(j)).                   
!        note:                                                                  
!     taulm is not a symmetrical matrix. for the array element                  
!     taulm(plm(i),plm(j)),the inner(first,most rapidly varying)                
!     dimension is the varying layer-mean pressure,plm(i);the outer             
!     (second) dimension is the fixed layer-mean pressure plm(j).               
!     thus the element taulm(2,3) is the transmission function between          
!     the fixed pressure plm(3)  and the pressure layer having an averag        
!     pressure of plm(2).                                                       
!         also note that no quadrature is performed over the layer              
!     between the smallest nonzero pressure and zero pressure;                  
!     taulm is taulm(0,plm(j)) in this case,and taulm(0,0)=1.                   
!                                                                               
!                                                                               
!             reference:                                                        
!         s.b.fels and m.d.schwarzkopf,-an efficient accurate                   
!     algorithm for calculating co2 15 um band cooling rates-,journal           
!     of geophysical research,vol.86,no. c2, pp.1205-1232,1981.                 
!        modifications to the algorithm have been made by the authors;          
!     contact s.b.f.or m.d.s. for further details.a note to j.g.r.              
!     is planned to document these changes.                                     
!                                                                               
!            author:    m.daniel schwarzkopf                                    
!                                                                               
!            date:      14 july 1983                                            
!                                                                               
!            address:                                                           
!                                                                               
!                      g.f.d.l.                                                 
!                      p.o.box 308                                              
!                      princeton,n.j.08540                                      
!                      u.s.a.                                                   
!            telephone:  (609) 452-6521                                         
!                                                                               
!            information on tape: this source is the first file                 
!        on this tape.the six files that follow are co2 trans-                  
!        missivities for the 500-850 cm-1 interval for co2                      
!        concentrations of 330 ppmv (1x) ,660 ppmv (2x), and                    
!        1320 ppmv (4x). the files are arranged as follows:                     
!          file 2   1x,consolidated using b(250) weighting fctn.                
!          file 3   1x,consolidated with no weighting fctn.                     
!          file 4   2x,consolidated using b(250) weighting fctn.                
!          file 5   2x,consolidated with no weighting fctn.                     
!          file 6   4x,consolidated using b(250) weighting fctn.                
!          file 7   4x,consolidated with no weighting fctn.                     
!            files 2,4,6 are recommended for use in obtaining                   
!        transmission functions for use in heating rate                         
!        computations;they correspond to the transmissivities                   
!        discussed in the 1980 paper.files 3,5,7 are provided                   
!        to facilitate comparison with observation and with other               
!        calculations.                                                          
!                                                                               
!            program language: fortran 1977,including parameter                 
!        and program statements.the program is written on a                     
!        cyber 170-730.see the section on limitations for                       
!        adaptations to other machines.                                         
!                                                                               
!          input units,formats and format statement nos:                        
!                                                                               
!   unit no    variables       format      statement no.    type                
!      5        p (purpose 1)  (5e16.9)        201         cards                
!      5        pd (purpose 2) (5e16.9)        201         cards                
!      5        plm(purpose 2) (5e16.9)        201         cards                
!      5        nmethd         (i3)            202         cards                
!      20       transa         (4f20.14)       102          tape                
!nov89                                                                          
!      itape    transa         (4f20.14)       102          tape                
!nov89                                                                          
!                                                                               
!         output units,formats and format statement nos:                        
!                                                                               
!   unit no    variables       format     statement no.                         
!      6         trnfct        (1x,8f15.8)     301         print                
!      22        trnfct        (4f20.14)       102          tape                
!                                                                               
!            parameter inputs:                                                  
!     a) nlevls    : nlevls is an (integer) parameter denoting                  
!        the number of nonzero pressure levels for purpose 1                    
!        or the number of nonzero layer pressures needed to                     
!        specify the pressure layers(purpose 2) in the output                   
!        grid. for example,in purpose 1,if p=0,100,1000,nlevls=2.               
!        if,in purpose 2,pd=0,100,500,1000,the number of nonzero                
!        pressure layers=2,so nlevls=2                                          
!           in the code as written,nlevls=40; the user should                   
!        change this value to a user-specified value.                           
!     b) nlp1,nlp2 : integer parameters defined as: nlp1=nlevls+1;              
!        nlp2=nlevls+2.                                                         
!           see limitations for code modifications if parameter                 
!        statements are not allowed on your machine.                            
!                                                                               
!            inputs:                                                            
!                                                                               
!     a) transa    : the 109x109 grid of transmission functions                 
!            transa is a  double precision real array.                          
!                                                                               
!           transa  is read from file 20. this file contains 3                  
!     records,as follows:                                                       
!        1)   transa, standard temperature profile                              
!        3)   transa, standard temperatures + 25 deg                            
!        5)   transa, standard temperatures - 25 deg                            
!                                                                               
!     b)   nmethd: an integer whose value is either 1 (if co2o3_step2 is             
!       to be used for purpose 1) or 2 (if co2o3_step2 is to be used for             
!       purpose 2).                                                             
!                                                                               
!     c)     p,pd,plm :                                                         
!          p is a real array (length nlp1) specifying the pressure              
!       grid at which transmission functions are to be computed for             
!       purpose 1.the dimension  of p is  in millibars.the                      
!       following limitations will be explained more                            
!       in the section on limitations: p(1) must be zero; p(nlp1),the           
!       largest pressure, must not exceed 1165 millibars.                       
!         pd is a real array (length nlp2) specifying the pressure              
!       layers for which layer-averaged transmission functions are              
!       to be computed.the dimension of pd is millibars.the limitations         
!       for pd are the same as for p,and are given in the section on            
!       limitations.                                                            
!         plm is a real array (length nlp2) specifying the layer-mean           
!       pressures. the dimension of plm is millibars. the limitations           
!       for plm are the same as for p,and are given in the section on           
!       limitations.pd is read in before plm.                                   
!                                                                               
!          note: again,we note that the user will input either p (for           
!       purpose 1) or pd and plm(for purpose 2) but not both.                   
!                                                                               
!                                                                               
!                                                                               
!                                                                               
!           limitations:                                                        
!     1)       p(1)=0.,pd(1)=0.,plm(1)=0. the top pressure level                
!       must be zero,or the top pressure layer must be bounded by zero.         
!       the top layer-mean pressure (plm(1)) must be zero; no                   
!       quadrature is done on the top pressure layer.even if one is             
!       not interested in the transmission function between 0 and p(j),         
!       one must include such a level.                                          
!     2)      pd(nlp2)=p(nlp1) is less than or equal to 1165 mb.                
!       extrapolation to higher pressures is not possible.                      
!     3)      if program is not permitted on your compiler,                     
!       simply delete the line.                                                 
!     4)      if parameter is not permitted,do the following:                   
!            1) delete all parameter statements in co2o3_step2                       
!            2) at the point where nmethod is read in,add:                      
!                read (5,202) nlevls                                            
!                nlp1=nlevls+1                                                  
!                nlp2=nlevls+2                                                  
!            3) change dimension and/or common statements defining              
!              arrays trns,delta,p,pd,trnfct,ps,pds,plm in co2o3_step2.              
!              the numerical value of (nlevls+1) should be inserted             
!              in dimension or common statements for trns,delta,                
!              p,trnfct,ps,plm; the numerical value of (nlevls+2)               
!              in dimension or common statements for pd,pds.                    
!      5)    parameter (nlevls=40) and the other parameter                      
!       statements are written in cdc fortran; on other machines the            
!       same statement may be written differently,for example as                
!       parameter   nlevls=40                                                   
!      6) -double precision- is used instead of -real*8- ,due to                
!       requirements of cdc fortan.                                             
!      7) the statement -do 400 kkk=1,3- controls the number of                 
!       transmissivity output matrices porduced by the program.to               
!       produce 1 output matrix,delete this statement.                          
!                                                                               
!     output:                                                                   
!         a) trnfct is an (nlp1,nlp1) real array of the transmission            
!     functions appropriate to your array. it is to be saved on file 22.        
!     the procedure for saving may be modified; as given here,the               
!     output is in card image form with a format of (4f20.14).                  
!                                                                               
!         b)  printed  output is a listing of trnfct on unit 6, in              
!     the format (1x,8f15.8) (format statement 301). the user may               
!     modify or eliminate this at will.                                         
!                                                                               
!      ************   function interpolater routine  *****************          
!                                                                               
!                                                                               
!     ******   the following parameter gives the number of     *******          
!     ******           data levels in the model                *******          
!     ****************************************************************          
!-------------------------------------------------------------------------------
      real     ::  ps(nlp1),pds(nlp2),plm(nlp1)
      real     ::  t15a(nlp2,2),t15b(nlp1)
      real     ::  t22(nlp1,nlp1,3)
      integer  ::  nrtab(3)                                                        
      data nrtab/1,2,4/                                                         
!
!   the following are the input formats                                         
!
   100   format (4f20.14)                                                          
   201   format (5e16.9)                                                           
   202   format (i3)                                                               
!o222   203   format (f12.6,i2)                                                 
   203   format (f12.6)                                                            
!
!    the following are the output formats                                       
!
   102   format (4f20.14)                                                          
   301   format (1x,8f15.8)                                                        
!                                                                               
!cc   rewind 15                                                                 
!cc   rewind 20                                                                 
!nov89                                                                          
   rewind itape                                                              
!nov89                                                                          
!cc   rewind 22                                                                 
!                                                                               
!     calculation of pa -the -table- of 109 grid pressures                      
!     note-this code must not be changed by the user
!
   pa(1)=0.                                                                  
   fact15=10.**(1./15.)                                                      
   fact30=10.**(1./30.)                                                      
   pa(2)=1.0e-3                                                              
!
   do i = 2,76                                                             
     pa(i+1)=pa(i)*fact15                                                      
   enddo
!
   do i = 77,108                                                           
     pa(i+1)=pa(i)*fact30                                                      
   enddo
!                                                                               
   n=25                                                                      
   nlv=nlevls                                                                
   nlp1v=nlp1                                                                
   nlp2v=nlp2                                                                
!
!     read in the co2 mixing ratio(in units of 330 ppmv),and an index           
!     giving the frequency range of the lbl data                                
!o222    read (5,203) ratio,ir                                                  
!cc         ir = 1                                                              
!cc         read (5,203) ratio                                                  
!o222   ***********************************                                     
!***values for ir*****                                                          
!          ir=1     consol. lbl trans. =490-850                                 
!          ir=2     consol. lbl trans. =490-670                                 
!          ir=3     consol. lbl trans. =670-850                                 
!          ir=4     consol. lbl trans. =2270-2380                               
!*** ir must be 1,2,3 or 4 for the pgm. to work                                 
!     also read in the method no.(1 or 2)                                       
!cc         read (5,202) nmethd                                                 
!oct92if (ratio.eq.1.0) go to 621                                               
!oct92if (ratio.ge.2.0) go to 622                                               
!oct92if (ratio.eq.4.0) go to 623                                               
!oct92if (ratio.gt.1.0.and.ratio.lt.2.0) go to 624                              
!oct92if (ratio.gt.2.0.and.ratio.lt.4.0) go to 625                              
!oct92stop 8746                                                                 
!nov89  621   itap1=20                                                          
!621   itap1=itape                                                               
!nov89                                                                          
!
   ntap=1                                                                    
   go to 630                                                                 
   622   itap1=21                                                                  
   ntap=1                                                                    
   go to 630                                                                 
   623   itap1=22                                                                  
   ntap=1                                                                    
   go to 630                                                                 
!nov89   624   itap1=20                                                         
   624   itap1=itape                                                               
!nov89                                                                          
   ntap=2                                                                    
   ratstd=2.0                                                                
   ratsm=1.0                                                                 
   go to 630                                                                 
   625   itap1=21                                                                  
   ntap=2                                                                    
   ratstd=4.0                                                                
   ratsm=2.0                                                                 
   630   continue                                                                  
   if (nmethd.eq.2) go to 502                                                
!
!   *****cards for purpose 1(nmethd=1)                                          
!cc         read (15,201) (p(i),i=1,nlp1)                                       
!
   do i = 1,nlp1                                                           
     p(i)=t15b(i)                                                            
   enddo
!
   do i = 1,nlp1                                                           
     ps(i)=p(i)                                                                
   enddo
!
   go to 503                                                                 
   502   continue                                                                  
!
!  *****cards for purpose 2(nmethd=2)                                           
!cc         read (15,201) (pd(i),i=1,nlp2)                                      
!cc         read (15,201) (plm(i),i=1,nlp1)                                     
!
   do i = 1,nlp2                                                           
      pd(i)=t15a(i,1)                                                         
   enddo
!
   do i = 1,nlp1                                                           
      plm(i)=t15a(i,2)                                                        
   enddo
!
   do i = 1,nlp1                                                           
      pds(i)=pd(i+1)                                                            
      ps(i)=plm(i)                                                              
   enddo
!                                                                               
   503   continue                                                                  
!
!  *****do loop controlling number of output matrices                           
!nov89                                                                          
!nov89    do 400 kkk = 1,3                                                        
!
   icloop = 3                                                                
   if (ir.eq.4) icloop = 1                                                   
!
   do kkk = 1,icloop                                                       
!nov89                                                                          
     if (ntap.eq.2) call co2_trans_funct(itap1,ratstd,ratsm,ratio) 
!
!  **********************                                                       
!
     if (nmethd.eq.2) go to 505                                                
!
!   *****cards for purpose 1(nmethd=1)                                          
!
     do i = 1,nlp1                                                           
       p(i)=ps(i)                                                                
     enddo
     go to 506                                                                 
     505  continue                                                                  
!
!  *****cards for purpose 2(nmethd=2)                                           
!
     do i = 1,nlp1                                                           
       pd(i)=pds(i)                                                              
       p(i)=ps(i)                                                                
     enddo
!                                                                               
     506   continue                                                                  
     ia=108                                                                    
     iap=ia+1                                                                  
!nov89   if (ntap.eq.1) read (20,100) ((transa(i,j),i=1,109),j=1,109)           
     if (ntap.eq.1) read (itape,100) ((transa(i,j),i=1,109),j=1,109)           
!nov89                                                                          
!
     do i = 1,iap                                                              
       transa(i,i)=1.0                                                           
     enddo
!
     call co2_tau_interp(ratio,ir)   
!
     do i = 1,nlp1                                                           
       do j = 1,nlp1                                                           
         trns(j,i)=1.00                                                            
       enddo
     enddo
!
     do i = 1,nlp1                                                            
       do 20 j=1,i                                                               
         if (i.eq.j) go to 20                                                      
         p1=p(j)                                                                   
         p2=p(i)                                                                   
         call co2_interp    
         trns(j,i)=trnslo                                                          
       20  continue                                                                  
     enddo
!
     do i = 1,nlp1                                                            
       do j = i,nlp1                                                            
         trns(j,i)=trns(i,j)                                                       
       enddo
     enddo
!
!  *****this is the end of purpose 1 calculations                               
!
     if (nmethd.eq.1) go to 2872                                               
!                                                                               
     do j = 1,nlp1                                                            
       do i = 2,nlp1                                                            
         ia=i                                                                      
         ja=j                                                                      
         n=25                                                                      
         if (i.ne.j) n=3                                                           
         call co2_quardratic(nlv,nlp1v,nlp2v,p,pd,trns) 
       enddo
     enddo
!
!  *****this is the end of purpose 2 calculations                               
!
      2872  continue                                                                  
!                                                                               
      write (6,301) ((trns(i,j),i=1,nlp1),j=1,nlp1)                             
!cc         write (22,102) ((trns(i,j),i=1,nlp1),j=1,nlp1)                      
!
     do j = 1,nlp1                                                           
       do i = 1,nlp1                                                          
         t22(i,j,kkk) = trns(i,j)                                                
       enddo
     enddo
   enddo
!
   return                                                                    
   end subroutine co2o3_step2
!-------------------------------------------------------------------------------
