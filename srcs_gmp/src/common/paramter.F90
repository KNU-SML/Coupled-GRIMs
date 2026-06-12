#include "define.h"
   module paramter
!-------------------------------------------------------------------------------
#ifdef SMP_RA2SFC
   use paramodel, only : idim  => lonfscm_                                    ,&
                         jdim  => latgscm_                                    ,&
#else
   use paramodel, only : idim  => lonf_                                       ,&
                         jdim  => latg_                                       ,&
#endif
                         mwave => jcap_                                       ,&
                         kdim  => levs_                                       ,&
                         kdimq => levs_                                       ,&
                         kdims => levs_                                       ,&
                         nptdm => lpnt_                                       ,&
                         nstdm => ltstp_ 
!
   end module paramter
