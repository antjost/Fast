            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent

            tcx  = tijk(lmtr,ic)*ci_mtr
            tcy  = tijk(lmtr,jc)*cj_mtr
            tcz  = tijk(lmtr,kc)*ck_mtr

          !on calcule l'inverse de la surface de la facette
          s_1  = c_pertu/sqrt(tcx*tcx + tcy*tcy +tcz*tcz)

          call random_number(rnd)
          rnd = (rnd-0.5)*ampli_transi*0.05                       !pertu random
     &       +2.*ampli_transi*sin(omega*param_real(TEMPS))
     &                       *sin(lambda*z(ldx))                  !pertu sinuoidale temps et z + gaussien (x,y)
     &        *exp(-10.*xlong*( (x(ldx)-x0)*(x(ldx)-x0)
     &                         +(y(ldx)-y0)*(y(ldx)-y0)))

            u  =  rop(ldjr,2) - tcx*rnd*s_1
            v  =  rop(ldjr,3) - tcy*rnd*s_1
            w  =  rop(ldjr,4) -     rnd*c_pertu

            rop(l,1) = rop(ldjr,1)
            rop(l,2) = ventx*c_ale - u
            rop(l,3) = venty*c_ale - v
            rop(l,4) = ventz*c_ale - w

            rop(l,5) = rop(ldjr,5)
