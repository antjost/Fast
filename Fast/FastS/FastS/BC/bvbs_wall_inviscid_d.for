C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.13 (r6666M) - 28 May 2018 09:28
C
C  Differentiation of bvbs_wall_inviscid in forward (tangent) mode:
C   variations   of useful results: rop
C   with respect to varying inputs: rop
C   RW status of diff variables: rop:in-out
C
C
C
C
C
C
C
C
C
C
C
C***********************************************************************
C     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
C     $Revision: 35 $
C     $Author: IvanMary $
C***********************************************************************
      SUBROUTINE BVBS_WALL_INVISCID_D(idir, lrhs, neq_mtr, mobile_coef, 
     +                                param_int, ind_loop, ventijk, tijk
     +                                , rop, ropd)
      IMPLICIT NONE
C
CSeule la valeur de k_vent et ck_vent a un sens dans cet appel
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
      INTEGER*4 idir, lrhs, neq_mtr, ind_loop(6), param_int(0:*)
C
      REAL*8 rop(param_int(41), param_int(36))
      REAL*8 ropd(param_int(41), param_int(36))
      REAL*8 ventijk(param_int(44), param_int(40))
      REAL*8 tijk(param_int(43), neq_mtr)
      REAL*8 mobile_coef
C Var local
      INTEGER*4 inc2, inc3, li1, li2, l, iref, jref, kref, njf, nkf, 
     +          ldjr, ic, jc, kc, i, j, k, li3, ldp, kc_vent, lmtr
      REAL*8 ci_mtr, cj_mtr, ck_mtr, ck_vent, c_ale, tcx, tcy, tcz, 
     +       ventx, venty, ventz, r, u, v, w, ut, vt, wt, ua, va, wa, 
     +       s_1, qn, surf
      REAL*8 ud, vd, wd, utd, vtd, wtd, uad, vad, wad, qnd
C
C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER_E inddm, i_1, j_1, k_1
C    adresse interface pour tableau metric
      INTEGER_E indmtr, i_3, j_3, k_3
C    adresse interface pour tableau vitesse entrainement
      INTEGER_E indven, i_4, j_4, k_4
      EXTERNAL SHAPE_TAB_MTR
      INTRINSIC SQRT
      INTRINSIC MAX
      REAL*8 arg1
C
C      write(*,*)'idir=', idir,nijk(4),nijk(5),ndimdx
C      write(*,*)'nijk=', nijk
C      write(*,*)'loop=', ind_loop
C
C
C......determine la forme des tableuz metrique en fonction de la nature du domaine
      CALL SHAPE_TAB_MTR(neq_mtr, param_int, idir, ic, jc, kc, kc_vent, 
     +                   ci_mtr, cj_mtr, ck_mtr, ck_vent, c_ale)
C
      c_ale = c_ale*2.*mobile_coef
      IF (lrhs .EQ. 1) c_ale = 0.
Cidir
C
      IF (idir .EQ. 1) THEN
C
        iref = 2*ind_loop(2) + 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (iref-i+param_int(0+3)-1) + (j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (ind_loop(2)+1+param_int(15+3)-1)*param_int(15
     +            ) + (j+param_int(15+3)-1)*param_int(15+1) + (k+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (ind_loop(2)+1+param_int(5+3)-1)*param_int(5)
     +            + (j+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5
     +            +4)-1)*param_int(5+2)
C
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (iref-i+param_int(0+3)-1) + (j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (ind_loop(2)+1+param_int(15+3)-1)*param_int(15
     +            ) + (j+param_int(15+3)-1)*param_int(15+1) + (k+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (ind_loop(2)+1+param_int(5+3)-1)*param_int(5)
     +            + (j+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5
     +            +4)-1)*param_int(5+2)
C
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      ELSE IF (idir .EQ. 2) THEN
C
C
        iref = 2*ind_loop(1) - 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (iref-i+param_int(0+3)-1) + (j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (ind_loop(1)+param_int(15+3)-1)*param_int(15) 
     +            + (j+param_int(15+3)-1)*param_int(15+1) + (k+param_int
     +            (15+4)-1)*param_int(15+2)
                lmtr = 1 + (ind_loop(1)+param_int(5+3)-1)*param_int(5) +
     +            (j+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5+4
     +            )-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (iref-i+param_int(0+3)-1) + (j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (ind_loop(1)+param_int(15+3)-1)*param_int(15) 
     +            + (j+param_int(15+3)-1)*param_int(15+1) + (k+param_int
     +            (15+4)-1)*param_int(15+2)
                lmtr = 1 + (ind_loop(1)+param_int(5+3)-1)*param_int(5) +
     +            (j+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5+4
     +            )-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      ELSE IF (idir .EQ. 3) THEN
C
C
C
        jref = 2*ind_loop(4) + 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (jref-j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (
     +            ind_loop(4)+1+param_int(15+3)-1)*param_int(15+1) + (k+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (ind_loop
     +            (4)+1+param_int(5+3)-1)*param_int(5+1) + (k+param_int(
     +            5+4)-1)*param_int(5+2)
C
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (jref-j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (
     +            ind_loop(4)+1+param_int(15+3)-1)*param_int(15+1) + (k+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (ind_loop
     +            (4)+1+param_int(5+3)-1)*param_int(5+1) + (k+param_int(
     +            5+4)-1)*param_int(5+2)
C
C
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      ELSE IF (idir .EQ. 4) THEN
C
C
        jref = 2*ind_loop(3) - 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
CDEC$ IVDEP
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (jref-j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(0+3)-1) + (ind_loop(3)+param_int(
     +            0+3)-1)*param_int(0) + (k+param_int(0+4)-1)*param_int(
     +            0)*param_int(0+1)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (ind_loop
     +            (3)+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5+
     +            4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (jref-j+param_int(0+3)
     +            -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (
     +            ind_loop(3)+param_int(15+3)-1)*param_int(15+1) + (k+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (ind_loop
     +            (3)+param_int(5+3)-1)*param_int(5+1) + (k+param_int(5+
     +            4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      ELSE IF (idir .EQ. 5) THEN
C
C
C
        kref = 2*ind_loop(6) + 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (kref-k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (j+
     +            param_int(15+3)-1)*param_int(15+1) + (ind_loop(6)+1+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (j+
     +            param_int(5+3)-1)*param_int(5+1) + (ind_loop(6)+1+
     +            param_int(5+4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (kref-k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (j+
     +            param_int(15+3)-1)*param_int(15+1) + (ind_loop(6)+1+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (j+
     +            param_int(5+3)-1)*param_int(5+1) + (ind_loop(6)+1+
     +            param_int(5+4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      ELSE
C
C
C
        kref = 2*ind_loop(5) - 1
Cparam_int(36)
C
        IF (param_int(36) .EQ. 5) THEN
C
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
CDEC$ IVDEP
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (kref-k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (j+
     +            param_int(15+3)-1)*param_int(15+1) + (ind_loop(5)+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (j+
     +            param_int(5+3)-1)*param_int(5+1) + (ind_loop(5)+
     +            param_int(5+4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
              ENDDO
            ENDDO
          ENDDO
        ELSE
C
          DO k=ind_loop(5),ind_loop(6)
            DO j=ind_loop(3),ind_loop(4)
              DO i=ind_loop(1),ind_loop(2)
C
                l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldjr = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +            param_int(0) + (kref-k+param_int(0+4)-1)*param_int(0)*
     +            param_int(0+1)
                ldp = 1 + (i+param_int(15+3)-1)*param_int(15) + (j+
     +            param_int(15+3)-1)*param_int(15+1) + (ind_loop(5)+
     +            param_int(15+4)-1)*param_int(15+2)
                lmtr = 1 + (i+param_int(5+3)-1)*param_int(5) + (j+
     +            param_int(5+3)-1)*param_int(5+1) + (ind_loop(5)+
     +            param_int(5+4)-1)*param_int(5+2)
                tcx = tijk(lmtr, ic)*ci_mtr
                tcy = tijk(lmtr, jc)*cj_mtr
                tcz = tijk(lmtr, kc)*ck_mtr
C
                arg1 = tcx*tcx + tcy*tcy + tcz*tcz
                surf = SQRT(arg1)
                IF (surf .LT. 1e-30) THEN
                  surf = 1e-30
                ELSE
                  surf = surf
                END IF
                s_1 = 1./surf
C
                tcx = tcx*s_1
                tcy = tcy*s_1
                tcz = tcz*s_1
C
                ventx = ventijk(ldp, 1)
                venty = ventijk(ldp, 2)
Con calcule l'inverse de la surface de la facette
                ventz = ventijk(ldp, kc_vent)*ck_vent
C
C
                ud = ropd(ldjr, 2)
                u = rop(ldjr, 2)
                vd = ropd(ldjr, 3)
                v = rop(ldjr, 3)
                wd = ropd(ldjr, 4)
                w = rop(ldjr, 4)
C
                qnd = tcx*ud + tcy*vd + tcz*wd
                qn = u*tcx + v*tcy + w*tcz
C
                utd = ud - tcx*qnd
                ut = u - qn*tcx
                vtd = vd - tcy*qnd
                vt = v - qn*tcy
                wtd = wd - tcz*qnd
                wt = w - qn*tcz
C
                uad = 2.*utd - ud
                ua = 2.*ut - u
                vad = 2.*vtd - vd
                va = 2.*vt - v
                wad = 2.*wtd - wd
                wa = 2.*wt - w
C
                ropd(l, 1) = ropd(ldjr, 1)
                ropd(l, 2) = uad
                ropd(l, 3) = vad
                ropd(l, 4) = wad
C
                ropd(l, 5) = ropd(ldjr, 5)
                ropd(l, 6) = ropd(ldjr, 6)
              ENDDO
            ENDDO
          ENDDO
        END IF
      END IF
      END

