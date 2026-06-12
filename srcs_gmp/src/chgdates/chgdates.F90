#include <define.h>
   program chgdates
!-------------------------------------------------------------------------------
!
! main program:  chgdates    change date inside file  
!                                                                               
! abstract: change date inside the file
!
! program history log:
!   1981-01-01  sela                   initial mrf
!   2000-01-01  hann-ming henry juang  mpi 
!   2000-01-01  song-you hong          physcis options 
!   2001-01-01  masao kanamitsu        seasonal forecast system
!   2005-01-01  young-hwa byun         single column model
!   2006-04-01  hoon park              double-fourier spectral (dfs)
!   2009-10-01  jung-eun kim           f90 format with standard physics modules
!   2010-07-01  myung-seo koo          dimension allocatble with namelist input
!
! namelists:                                                                    
!   namin:      parameters determining new date
!  
! input files:                                                                  
!   unit   11  sigma file(s)                                                  
!                                                                               
! output files:                                                                 
!   unit   51  sigma file                                         
!
!-------------------------------------------------------------------------------
   use paramodel, only : para_init,lonf_,latg_,jcap_,levs_,ntotal_,igen_
   use comsfc, only : sfcftyp
!-------------------------------------------------------------------------------
   character(LEN=8)    ::  on85lab(4)
   integer, parameter  ::  kdum2=21,kens=2
   integer             ::  idate(4)
!
   integer             ::  ijdim
   integer             ::  nwave
   integer             ::  kdum
   real, allocatable   ::  si(:),sl(:)
   real, allocatable   ::  grid(:,:)
   real, allocatable   ::  wave(:)
   real, allocatable   ::  dummy(:)
   real                ::  dummy2(kdum2),ensemble(kens)
!
   character(len=8), allocatable    ::  gvar(:)
   integer                          ::  nrecs,maxlev
   integer, allocatable             ::  lev(:)
!
   data newyr,newmo,newdy,newhr,fhnew/-1,-1,-1,-1,-1./
   namelist/namin/newyr,newmo,newdy,newhr,fhnew
!-------------------------------------------------------------------------------
   call para_init
#ifdef NCO_TAG
   call w3tagb('clim_chgdates',2001,0000,0000,'np51')
!
#endif
   ijdim=lonf_*latg_
   nwave=(jcap_+1)*(jcap_+2)
   kdum=201-levs_-1-levs_
   allocate(si(levs_+1),sl(levs_))
   allocate(wave(nwave))
   allocate(dummy(kdum))
!
   read(5,namin)
   print*, newyr,newmo,newdy,newhr,fhnew
!
!  sigma file
!
   read(11,end=901,err=901) on85lab
   write(51) on85lab
!
!  read(11) fhour,idate,(si(k),k=1,levs_+1),(sl(k),k=1,levs_)
!
   read(11,err=201) fhour,idate,(si(k),k=1,levs_+1),(sl(k),k=1,levs_)          &
   ,dummy,waves,xlayers,trun,order,realform,gencode                            &
   ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                                  &
   ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                       &
   ,pdryini,dummy2,gases
   if(water.eq.0.) then
     water=1.
     print *,'water reset to 1.'
   endif
   iwater=nint(water)
   igases=nint(gases)
   print *,'tread unit fhour,idate=',11,fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
   goto 202
   201   continue
   rewind 11
   read(11) on85lab
   read(11,err=201) fhour,idate,(si(k),k=1,levs_+1),(sl(k),k=1,levs_)
!
   do i = 1,kdum
     dummy(i)=0.
   enddo
   waves=jcap_
   xlayers=levs_
   trun=1.
   order=2.
   realform=1.
   gencode=igen_
   rlond=lonf_
   rlatd=latg_
   rlonp=lonf_
   rlatp=latg_
   rlonr=lonf_
   rlatr=latg_
   gases=0.
   water=1.
   pdryini=0.
   subcen=0.
   ensemble(1)=0.
   ensemble(2)=0.
   ppid=0.
   slid=0.
   vcid=0.
   vmid=0.
   vtid=0.
   do k = 1,kdum2
     dummy2(k)=0.
   enddo
!
   igases=nint(gases)
   iwater=nint(water)
   print *,'tread old format fhour,idate=',fhour,idate
   print *,' number of water input = ',iwater
   print *,' number of gases input = ',igases
   202   continue
!
   if(newyr.lt.0 ) newyr=idate(4)
   if(newmo.lt.0 ) newmo=idate(2)
   if(newdy.lt.0 ) newdy=idate(3)
   if(newhr.lt.0 ) newhr=idate(1)
   if(fhnew.lt.0.) fhnew=fhour
! 
!     idate(4)=year, (2)=month, (3)=day, (1)=hour
!
   write(6,*) 'year,month,day,hour,fhour of sigma file' 
   write(6,*) idate(4),idate(2),idate(3),idate(1),fhour
   write(6,*) 'replaced by'
   write(6,*) newyr,newmo,newdy,newhr,fhnew
!
   idate(4)=newyr
   idate(2)=newmo
   idate(3)=newdy
   idate(1)=newhr
   fhour=fhnew
!
   write(51) fhour,idate,(si(k),k=1,levs_+1),(sl(k),k=1,levs_)                 &
   ,dummy,waves,xlayers,trun,order,realform,gencode                            &
   ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,water                                  &
   ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid                       &
   ,pdryini,dummy2,gases
!
   do k = 1,2+3*levs_+ntotal_*levs_
     read(11,end=900) (wave(i),i=1,nwave)
     write(51) (wave(i),i=1,nwave)
     write(6,*) ' write k=',k,' wave=',wave(1)
   enddo
   go to 900
!
   901 continue
   print *,'<WARNING> sig file empty'
   900 continue
!
!  sfc file
!
   write(6,*) ' '
   read(12,end=910,err=910) on85lab
   write(52) on85lab
!
   read(12) fhour,idate
!
   write(6,*) 'year,month,day,hour,fhour of surface file' 
   write(6,*) idate(4),idate(2),idate(3),idate(1),fhour
   write(6,*) 'replaced by'
   write(6,*) newyr,newmo,newdy,newhr,fhnew
!
   idate(4)=newyr
   idate(2)=newmo
   idate(3)=newdy
   idate(1)=newhr
   fhour=fhnew
!
   write(52) fhour,idate
!
!  get sfc file properties
!
   call sfcfld(sfcftyp,0,nrecs,lev,gvar,maxlev)
   allocate (lev(nrecs),gvar(nrecs))
   call sfcfld(sfcftyp,1,nrecs,lev,gvar,maxlev)
!
   kend=nrecs
   do k = 1,kend
     allocate(grid(ijdim,lev(k)))
     read(12,end=909) ((grid(ij,l),ij=1,ijdim),l=1,lev(k))
     write(52) ((grid(ij,l),ij=1,ijdim),l=1,lev(k))
     deallocate(grid)
!
     write(6,*) ' write k=',k
   enddo
!
   909   continue
   stop
   910   continue
   print *,'<WARNING> sfc file empty'

#ifdef NCO_TAG
   CALL w3tage('clim_chgdates')
#endif
   deallocate(si,sl,wave,dummy)
   stop
   end
!-------------------------------------------------------------------------------
