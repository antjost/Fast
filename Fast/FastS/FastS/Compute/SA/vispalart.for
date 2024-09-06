c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine vispalart(ndom,  param_int, param_real, ind_loop,
     &                     xmut,  rop, dlng, ro_cutOff)
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Calcul de la contribution des termes sources 
c_A    pour l'equation de Spalart Allmaras 
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    Navier-Stokes
c
c     INP
c_I   rop: var primitive
c
c     OUT
c     xmut
c***********************************************************************
!! xmut = mu + mu_t
!! Notes: This calculates mu_t=rho*nutilde*fv1 for SA only. Not for ZDES
!! ZDES is done in the source term approach.
      implicit  none
      
#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6), param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E  rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E dlng(param_int(NDIMDX)), ro_cutOff(param_int(NDIMDX))

      REAL_E param_real(0:*)
c Var loc 
      INTEGER_E incmax,i,j,k,l,ltij,lij,lt,lvo

      REAL_E amulam,anulam,fv1,fvv1,anutild,amut,xmuprov,
     &     chi,ad1, s,temp01,cmus1,coesut,t1_1
      REAL_E dist, cutOffDist, DistCoef
      INTEGER_E cutOffExist
      
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

      include "omp_lib.h"
            

c.....formulation originelle
      fv1(s)     = 1./(1.+SA_CV1/(s*s*s))

      cutOffExist = param_int(CUTOFF)
      
      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

      !if (OMP_get_thread_num()==0) then
      !   write(*,*)'cutOffExist=',cutOffExist
      !end if
      !
      !
      !Modele Spalart
      !
      !

      IF (param_int(SA_INT + SA_IDES-1 ).le.1) THEN !SA

#include     "FastS/Compute/loop_begin.for" 
       
#include     "FastS/Compute/mulam.for" 
#include     "FastS/Compute/SA/chi.for"
              if (cutOffExist==1) then
                 DistCoef = ro_cutOff(l)
              else
                 DistCoef = 1
              end if
              
              fvv1     = DistCoef*fv1(DistCoef*chi)
#include     "FastS/Compute/SA/xmut.for"

#include     "FastS/Compute/loop_end.for" 

      ENDIF


      end
