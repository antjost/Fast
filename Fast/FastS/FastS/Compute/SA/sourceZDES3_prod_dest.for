!!  amulam :: μ kinematic viscosity based on Sutherland
#include  "FastS/Compute/mulam.for"
!! anutild :: ν turbulent SA nu tilde  || chi:: ν/v(kinematic viscosity)
#include  "FastS/Compute/SA/chi.for"
!choix distance
!!ad1 :: d_w wall distance
          ad1     = max( dlng(l), 1.e-27)                   !! d_w

          testfa  = 0.5+ sign(.5,SA_RCTEDES*adelta1 - ad1 ) !if testfa =1, everything makes sense

          !!testfa=0.0 || zgris ==1 !LES
          !!testfa=1.0 || zgris ==0 !RANS
          testfa=1-zgris(l)          
#include  "FastS/Compute/SA/fv1Gen.for"          
#include  "FastS/Compute/SA/xmut.for"
          fvv1    = (1.-testfa)+testfa*fvv1

          adtild1 = ad1+zgris(l)*(min( ad1 , SA_RCTEDES*adelta1) - ad1)  !! d_DES^III
          dist    = max( adtild1, 1.e-27)                                !! d_DES^III

          !CALCUL DU TERME DE PRODUCTION
          f2       = fv2(chi,fvv1)*testfa                             ! fv_2=1-\chi/(1+\chi*fv_1) [Eqn. (8)]
#include  "FastS/Compute/SA/zdes2_zdes3_common.for"
          fwg       = (1.-testfa) + testfa*fw(g)                      ! [Eqn. (9)]
#include  "FastS/Compute/SA/zdes2_zdes3_common_prt2.for"



!!Papers for equations: Recent improvements in the formulation of mode III of ZDES for WMLES use at Re_theta>10^4
!!note: fi     = fdelta
!!      fdelta = zgris          
!!      fdelta = 1-testfa           
