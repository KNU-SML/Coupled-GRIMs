!
   module comznl
!-------------------------------------------------------------------------------
   use paramodel, only   :  levs_,LATG2S
!-------------------------------------------------------------------------------
   private              ::  levs_,LATG2S
!
   integer,parameter    ::  nrm=23,nrs=32
   integer,parameter    ::  nlb=6,nst=6
   integer,parameter    ::  nmu=1,nmv=2,nmtv=3,nmq=4,nmvot2=5
   integer,parameter    ::  nmdiv2=6,nmomega=7,nmt=8,nmrh=9,nmke=10
   integer,parameter    ::  nmtconv=11,nmtlarg=12,nmtshal=13,nmtvrdf=14
   integer,parameter    ::  nmqshal=16,nmqvrdf=17,nmuvrdf=18,nmvvrdf=19
   integer,parameter    ::  nmqconv=15,nmthsw=20,nstskin=10
   integer,parameter    ::  nmthlw=21,nmtcld=22,nmtccv=23
   integer,parameter    ::  nsrain=1,nsrainc=2,nstsfc=3,nsqsfc=4
   integer,parameter    ::  nsusfc=5,nsvsfc=6,nsrcov=7,nsrcovc=8,nsps=9
   integer,parameter    ::  nswet=11,nssnow=12,nstg1=13,nstg2=14,nstg3=15
   integer,parameter    ::  nssfcsw=16,nssfclw=17,nsrhs=18,nstvs=19,nsts=20
   integer,parameter    ::  nsqs=21,nsz0cm=22,nsslmsk=23,nsugwd=24,nsvgwd=25
   integer,parameter    ::  nsuasfc=26,nsuagwd=27,nsuamtn=28,nsua=29,nsuap=30
   integer,parameter    ::  nsep=31,nscldwrk=32
   real                 ::  zhm(nrm)
   real, allocatable    ::  zdm(:,:,:,:),zwm(:,:)
   real, allocatable    ::  zds(:,:,:,:),zws(:,:,:)
!
   contains
!-------------------------------------------------------------------------------
   subroutine comznl_init
!-------------------------------------------------------------------------------
   allocate(zdm(2,levs_,nrm,LATG2S),zwm(2,LATG2S))
   allocate(zds(2,nst,nrs,LATG2S),zws(2,nst,LATG2S))
!
   end subroutine comznl_init
!-------------------------------------------------------------------------------
   end module comznl
