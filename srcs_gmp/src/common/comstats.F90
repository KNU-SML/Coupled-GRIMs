#include <machine.h>
   module comstats
!-------------------------------------------------------------------------------
   real     ::  avgv(21,2), stdv(21,2),                                        &
                avgh(21,2), stdh(21,2),                                        &
                avgi(21,2), stdi(21,2),                                        &
                avghy(21),  stdhy(21),                                         &
                avgb,       stdb,                                              &
                skv(21,2), xkv(21,2),                                          &
                skh(21,2), xkh(21,2),                                          &
                ski(21,2), xki(21,2),                                          &
                skhy(21),  xkhy(21),                                           &
                skb,       xkb
   integer  ::  nv(21,2),                                                      &
                nh(21,2),                                                      &
                ni(21,2),                                                      &
                nhy(21),                                                       &
                nb        
!
   end module comstats
