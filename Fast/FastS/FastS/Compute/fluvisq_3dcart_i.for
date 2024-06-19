        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        il   = l  -  inci
        ir   = l

        l220 = l  -  inci  -  incj
        l020 = l           -  incj
        l210 = l  -  inci  +  incj
        l010 = l           +  incj
        l202 = l  -  inci            -  inck
        l002 = l                     -  inck
        l201 = l  -  inci            +  inck
        l001 = l                     +  inck

CCCCCC
CCCCCC
CCCCCC   Facette I
CCCCCC
CCCCCC
      !Calcul des vecteurs surfaces I
      tix  = tcx
      tjy  = .25 * tcy
      tkz  = .25 * tcz

      ! --- Gradients en U_i
      f1 = rop(il +v2)
      f2 = rop(ir +v2)
      fv = f1 + f2
      f3 = fv + rop(l220 +v2)+rop(l020 +v2)
      f4 = fv + rop(l210 +v2)+rop(l010 +v2)
      f5 = fv + rop(l202 +v2)+rop(l002 +v2)
      f6 = fv + rop(l201 +v2)+rop(l001 +v2)

      gradU_nx=(f2-f1)*tix
      gradU_ny=(f4-f3)*tjy
      gradU_nz=(f6-f5)*tkz

      !--- Gradients en v_i
      f1 = rop(il +v3)
      f2 = rop(ir +v3)
      fv = f1 + f2
      f3 = fv + rop(l220 +v3)+rop(l020 +v3)
      f4 = fv + rop(l210 +v3)+rop(l010 +v3)

      gradV_nx=(f2-f1)*tix
      gradV_ny=(f4-f3)*tjy

      !--- Gradients en W_i
      f1 = rop(il +v4)
      f2 = rop(ir +v4)
      fv = f1 + f2
      f5 = fv + rop(l202 +v4)+rop(l002 +v4)
      f6 = fv + rop(l201 +v4)+rop(l001 +v4)

      gradW_nx=(f2-f1)*tix
      gradW_nz=(f6-f5)*tkz

      !--- Gradients en T_i
      f1 = rop(il +v5)
      f2 = rop(ir +v5)

      gradT_nx=(f2-f1)*tix

      !--- assemblage flux
      div = ( gradU_nx+gradV_ny+ gradW_nz)*cvisq
  
      f1 =  gradU_nx - div
      f2 =  gradU_ny + gradV_nx
      f3 =  gradU_nz + gradW_nx

#include  "FastS/Compute/mut_interface.for"
      fvWMLES1=0
      fvWMLES2=0
      fvWMLES3=0
      if (int(param_real(WL_IBM_SWTCH))>30) then
         gradtmodel_nxLocal= 0.5*(rop_tij_model(ir+v1)+
     &                            rop_tij_model(il+v1)) !t11
         fvWMLES1          = -(gradtmodel_nxLocal*tcx)
     
         gradtmodel_nxLocal= 0.5*(rop_tij_model(ir+v2)+
     &                            rop_tij_model(il+v2)) !t21 
         fvWMLES2          = -(gradtmodel_nxLocal*tcx)
     
         gradtmodel_nxLocal= 0.5*(rop_tij_model(ir+v4)+
     &                            rop_tij_model(il+v4)) !t31 
         fvWMLES3          = -(gradtmodel_nxLocal*tcx)
      end if
      fv     = -(2.*f1*tcx )*xmutvol + fvWMLES1
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     = -(   f2*tcx )*xmutvol + fvWMLES2
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     = -(   f3*tcx )*xmutvol + fvWMLES3
      fv5    = fv5  + fv * (rop(ir+v4)+rop(il +v4))
      flu4   = flu4 + fv

      fv     = fv5*0.5  - ( gradT_nx*tcx  )*xktvol
      flu5   = flu5 + fv

