!
   subroutine fgb1(nflx,nfgb,fhour)
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_, latg_
!-------------------------------------------------------------------------------
   save
   real                        ::  work(lonf_,latg_)
   integer, parameter          ::  lenpds=28,lengds=32
   integer, parameter          ::  kchunk=1
   integer, parameter          ::  nchi=3, nxchi=3
!
   character             ::  grib(30+lenpds+lengds+lonf_*latg_*(32+1)/8,kchunk)
   character(len=5), parameter ::  fxno='fgb'
   character(len=5), parameter ::  fxni='flx'
   character(len=80)           ::  fno
!
   logical                     ::  lbm(lonf_*latg_)
   integer                     ::  lgrib(kchunk)
   integer                     ::  iens(5)
!-------------------------------------------------------------------------------
   call file_name(fxno,nchi,fhour,fno,ncho)
   open(unit=nfgb,file=fno(1:ncho),form='unformatted',err=900)
!
   go to 901
!
900 continue
   print *,' error in opening file ',fno(1:ncho)
   call abort
901 continue 
   print *,' opening file ',fno(1:ncho)
!
! tabl   kpds ipds name    kpds ipds name      kpds ipd names
!                1            7   11 il2k       17   21 ina
!                2            8   12 iyr        20   22 inm
!           1    3 icen       9   13 imo        21   23 icnt
!           2    4 igen      10   14 idy        23   24 icen2
!           3    5 igrid     11   15 ihr        22   25 idsk
!           4    6 iflg      12   16 imin
!                7 ibms      13   17 iftu
!           5    8 ipuk      14   18 ip1
!           6    9 itlk      15   19 ip2
!           7   10 il1k      16   20 itr
!
! readux output from model and gribit
!
   call file_name(fxni,nxchi,fhour,fno,ncho)
   open(unit=nflx,file=fno(1:ncho),form='unformatted',err=910)
!
   go to 911
!
910 continue
   print *,' error in opening file ',fno(1:ncho)
   call abort
911 continue 
   print *,' opening file ',fno(1:ncho)
   icheck=0
   kf=lonf_*latg_
#define WORK work
   do n = 1,100
     print *,'processing n=',n,' r_flx record'
     read(nflx,end=922,err=921)                                                &
                  WORK,lbm,idrt,io,jo,mxbit,colat1,                            &
                  lnpds,iptv,icn,igen,ibms,                                    &
                  ipuk,itlk,il1k,il2k,                                         &
                  iyr,imo,idy,ihr,iftu,ip1,ip2,itr,                            &
                  ina,inm,icen2,idsk,iens,                                     &
                  rlat1,rlon1,rlat2,rlon2,                                     &
                  delx,dely,ortru,proj
     call print_maxmin_six(work,lonf_*latg_,1,1,1,'flx file')
     call file_make_grib(WORK,lbm,idrt,io,jo,mxbit,colat1,                     &
                lnpds,iptv,icen,igen,ibms,                                     &
                ipuk,itlk,il1k,il2k,                                           &
                iyr,imo,idy,ihr,iftu,ip1,ip2,itr,                              &
                ina,inm,icen2,idsk,iens,                                       &
                rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj,truth,cotru,      &
                grib(1,1),lgrib(1),ierr)
     if(lgrib(1).gt.0) then
       call wryte(nfgb,lgrib(1),grib(1,1))
     endif
!
     go to 920
!
921  print *,' error reading flx file... but continue '
     call exit(10)
920  continue
   enddo
922 print *,' end of read flx file '
!
   close(nfgb)
   close(nflx)
!
   return
   end
!-------------------------------------------------------------------------------
