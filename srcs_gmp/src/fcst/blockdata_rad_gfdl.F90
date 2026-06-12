!
   module blockdata_rad_gfdl
!-------------------------------------------------------------------------------
!
!     block data intializes quantities needed by the gfdl codes.
!     bd2,bd3,bd4,bd5,blckfs all combined into 1 blockdata for frontend.
!
!     block data bd1 gives input data (temps,pressures,mixing ratios,
!     cloud amts and heights) for testing the radiation code as a
!     standalone model.
!
!     unused data cleaned out - nov 86 and mar 89 ..k.a.campana....
!
!-------------------------------------------------------------------------------
   use hcon
   use rdparm, only  :  nbly,nb
   use rnddta
   use co2dta, only  :  b0,b1,b2,b3
!-------------------------------------------------------------------------------
   real                 ::  delcm(nbly)         !  common/tbltmp/
   integer              ::  season,lseason,jtime(5),ljtime,jdnmc,ljdnmc
   integer              ::  ixxxx,lixxxx
   real                 ::  fcstda,daz(12),fjdnmc,tslag,rlag,timin
   real                 ::  tpi,hpi,year,day,dhr

   data delcm  /                                                               &
      0.300000e+02,  0.110000e+03,  0.600000e+02,  0.400000e+02,               &
      0.200000e+02,  0.500000e+02,  0.400000e+02,  0.500000e+02,               &
      0.110000e+03,  0.130000e+03,  0.100000e+03,  0.900000e+02,               &
      0.800000e+02,  0.130000e+03,  0.110000e+03/
!
! the following data are level-independent
! data rco2/3.3e-4/
! data ctauda/.5/
! data csolar/1.96/
! data ccosz/.5/
!
! for seasonal variation
!     season=1,2,3,4 for winter,spring,summer,fall only (not active)
!     season=5 - seasonal variation(i.e.interpolate to day of fcst)
!
   data season/5/
   data tslag/45.25/,  rlag/14.8125/
   data day/86400./,  year/365.25/
   data tpi/6.283185308/,  hpi/1.570796327/
   data jtime/0,1,0,0,0/
   data dhr/2./
   data daz/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
!
!     sea surface albedo data
!
   real,target          ::  albd(21,20)
   real                 ::  za(20),trn(21),dza(19)
   real,pointer         ::  alb1(:,:),alb2(:,:),alb3(:,:)
   data albd / .061,.062,.072,.087,.115,.163,.235,.318,.395,.472,.542,         &
    .604,.655,.693,.719,.732,.730,.681,.581,.453,.425,.061,.062,.070,          &
    .083,.108,.145,.198,.263,.336,.415,.487,.547,.595,.631,.656,.670,          &
    .652,.602,.494,.398,.370,.061,.061,.068,.079,.098,.130,.174,.228,          &
    .290,.357,.424,.498,.556,.588,.603,.592,.556,.488,.393,.342,.325,          &
    .061,.061,.065,.073,.086,.110,.150,.192,.248,.306,.360,.407,.444,          &
    .469,.480,.474,.444,.386,.333,.301,.290,.061,.061,.065,.070,.082,          &
    .101,.131,.168,.208,.252,.295,.331,.358,.375,.385,.377,.356,.320,          &
    .288,.266,.255,.061,.061,.063,.068,.077,.092,.114,.143,.176,.210,          &
    .242,.272,.288,.296,.300,.291,.273,.252,.237,.266,.220,.061,.061,          &
    .062,.066,.072,.084,.103,.127,.151,.176,.198,.219,.236,.245,.250,          &
    .246,.235,.222,.211,.205,.200,                                             &
              .061,.061,.061,.065,.071,.079,.094,.113,.134,.154,.173,          &
    .185,.190,.193,.193,.190,.188,.185,.182,.180,.178,.061,.061,.061,          &
    .064,.067,.072,.083,.099,.117,.135,.150,.160,.164,.165,.164,.162,          &
    .160,.159,.158,.157,.157,.061,.061,.061,.062,.065,.068,.074,.084,          &
    .097,.111,.121,.127,.130,.131,.131,.130,.129,.127,.126,.125,.122,          &
    .061,.061,.061,.061,.062,.064,.070,.076,.085,.094,.101,.105,.107,          &
    .106,.103,.100,.097,.096,.095,.095,.095,.061,.061,.061,.060,.061,          &
    .062,.065,.070,.075,.081,.086,.089,.090,.088,.084,.080,.077,.075,          &
    .074,.074,.074,.061,.061,.060,.060,.060,.061,.063,.065,.068,.072,          &
    .076,.077,.076,.074,.071,.067,.064,.062,.061,.061,.061,.061,.061,          &
    .060,.060,.060,.060,.061,.062,.065,.068,.069,.069,.068,.065,.061,          &
    .058,.055,.054,.053,.052,.052,                                             &
              .061,.061,.060,.060,.060,.060,.060,.060,.062,.065,.065,          &
    .063,.060,.057,.054,.050,.047,.046,.045,.044,.044,.061,.061,.060,          &
    .060,.060,.059,.059,.059,.059,.059,.058,.055,.051,.047,.043,.039,          &
    .035,.033,.032,.031,.031,.061,.061,.060,.060,.060,.059,.059,.058,          &
    .057,.056,.054,.051,.047,.043,.039,.036,.033,.030,.028,.027,.026,          &
    .061,.061,.060,.060,.060,.059,.059,.058,.057,.055,.052,.049,.045,          &
    .040,.036,.032,.029,.027,.026,.025,.025,.061,.061,.060,.060,.060,          &
    .059,.059,.058,.056,.053,.050,.046,.042,.038,.034,.031,.028,.026,          &
    .025,.025,.025,.061,.061,.060,.060,.059,.058,.058,.057,.055,.053,          &
    .050,.046,.042,.038,.034,.030,.028,.029,.025,.025,.025/
   data za/90.,88.,86.,84.,82.,80.,78.,76.,74.,70.,66.,62.,58.,54.,            &
           50.,40.,30.,20.,10.,0.0/
   data trn/.00,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,           &
            .70,.75,.80,.85,.90,.95,1.00/
   data dza/8*2.0,6*4.0,5*10.0/
!
   real                 ::  sc
   real                 ::  abcff(nb),pwts(nb),cfco2,cfo3,reflo3,rrayav
!  data sc/2.0/
!
!  changed to conform amip-ii specificatin of 1365 w/m**2
!
   data sc/1.9964498/
!
!---specification of data statements:
!         abcff=absorption coefficients for bands in k-distri-
!     bution. originally given by lacis and hansen, revised by
!     ramaswamy
!         pwts=corresponding weights assigned to bands in the
!     k-distribution
!         reflo3,rrayav= reflection coefficients given by
!     lacis and hansen to account for effects of rayleigh scattering
!     in the visible frequencies (band 1)
!         cfco2,cfo3=conversion factors from gm/cm**2 to cm-atm(stp)
!
!---the following are the coefficients for the 12-band shortwave
!   radiation code, specified by ramaswamy.
!
   data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190.,                     &
              989.,2706.,39011./
   data pwts/.5000,.121416,.0698,.1558,.0631,.0362,.0243,.0158,.0087,          &
             .001467,.002342,.001075/
!
!---the original 9-band lacis-hansen coefficients are given here; it
!   the user insists on using these values, she must also change
!   the parameter nb from 12 to 9. this parameter is defined in
!   rdparm.h . no other changes are required!
!     data abcff/2*4.0e-5,.002,.035,.377,1.95,9.40,44.6,190./
!     data pwts/.5000,.1470,.0698,.1443,.0584,.0335,.0225,.0158,.0087/
!
   data cfco2,cfo3/508.96,466.64/
   data reflo3/1.9/
   data rrayav/0.144/
!
   contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine init_rad_gfdl_blockdata
!-------------------------------------------------------------------------------
   alb1=>albd(1:21,1:7)
   alb2=>albd(1:21,8:14)
   alb3=>albd(1:21,15:20)
!
   return
   end subroutine init_rad_gfdl_blockdata
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine init_coef_rad_gfdl_blockdata
!-------------------------------------------------------------------------------
!
!   b0,b1,b2,b3 are coefficients used to correct for the use of 250k in
!   the planck function used in evaluating planck-weighted co2
!   transmission functions. (see ref. 4)
!
   b0=-.51926410e-4
   b1=-.18113332e-3
   b2=-.10680132e-5
   b3=-.67303519e-7
!
   return
   end subroutine init_coef_rad_gfdl_blockdata
!-------------------------------------------------------------------------------
   end module blockdata_rad_gfdl
