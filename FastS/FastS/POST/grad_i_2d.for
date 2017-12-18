        l5 = l+inci
        l6 = l-inci

        tix = ti(lt,1)
        tiy = ti(lt,2)

        !!GradU
        u =(rop(l +V2)+rop(l6 +V2))*c1- (rop(l5+V2)+rop(l-inc2i+V2))*c2

        dvardc( l  + Vdudx ) = dvardc( l  + Vdudx ) - u*tix
        dvardc( l0 + Vdudx ) = dvardc( l0 + Vdudx ) + u*tix
        dvardc( l  + Vdudy ) = dvardc( l  + Vdudy ) - u*tiy
        dvardc( l0 + Vdudy ) = dvardc( l0 + Vdudy ) + u*tiy

        !!GradV
        u =(rop(l +V3)+rop(l6 +V3))*c1- (rop(l5+V3)+rop(l-inc2i+V3))*c2

        dvardc( l  + Vdvdx ) = dvardc( l  + Vdvdx ) - u*tix
        dvardc( l0 + Vdvdx ) = dvardc( l0 + Vdvdx ) + u*tix
        dvardc( l  + Vdvdy ) = dvardc( l  + Vdvdy ) - u*tiy
        dvardc( l0 + Vdvdy ) = dvardc( l0 + Vdvdy ) + u*tiy


