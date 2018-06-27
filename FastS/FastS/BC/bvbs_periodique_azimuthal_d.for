C        Generated by TAPENADE     (INRIA, Ecuador team)
C  Tapenade 3.13 (r6666M) - 28 May 2018 09:28
C
C  Differentiation of bvbs_periodique_azimuthal in forward (tangent) mode:
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
C     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
C     $ $Revision: 56 $ 
C     $Author: IvanMary $
C*****a*****************************************************************
      SUBROUTINE BVBS_PERIODIQUE_AZIMUTHAL_D(idir, lrhs, param_int, 
     +                                       ind_loop, rop, ropd, 
     +                                       data_per)
      IMPLICIT NONE
C
C matrice de rotation
C
C
C
C      write(*,*)'idir=', idir,inc1,nijk(4),nijk(5),ndimdx
C      write(*,*)'nijk=', nijk
C      write(*,*)'loop=', ind_loop
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
      INTEGER*4 idir, lrhs, ind_loop(6), param_int(0:*)
C
      REAL*8 data_per(*)
      REAL*8 rop(param_int(41), param_int(36))
      REAL*8 ropd(param_int(41), param_int(36))
C
C Var local
      INTEGER*4 l, ijkm, im, jm, km, ldjr, i, j, k, inc1, ne, lij
      INTEGER*4 dir_axe
      REAL*8 c1
      REAL*8 pi, angle, cosa, sina
      REAL*8 ropr2, ropr3, ropr4
      REAL*8 ropr2d, ropr3d, ropr4d
      REAL*8, DIMENSION(3, 3) :: rot
C
C    adresse point courant pour tableau de la taille d'un domaine 
      INTEGER*4 inddm, i_1, j_1, k_1
      INTRINSIC ACOS
      INTRINSIC COS
      INTRINSIC SIN
Ccorrection angle
C
      IF (data_per(2) .EQ. 0. .AND. data_per(3) .EQ. 0.) THEN
        angle = data_per(1)
        dir_axe = 1
      ELSE IF (data_per(1) .EQ. 0. .AND. data_per(3) .EQ. 0.) THEN
        angle = data_per(2)
        dir_axe = 2
      ELSE IF (data_per(1) .EQ. 0. .AND. data_per(2) .EQ. 0.) THEN
        angle = data_per(3)
        dir_axe = 3
      END IF
      angle = angle*(-1.)
C
      pi = ACOS(-1.)
      cosa = COS(pi*angle/180.)
      sina = SIN(pi*angle/180.)
C mise a zero si tableau drodm du RHS implicit
      SELECT CASE (dir_axe) 
      CASE (1) 
        rot(1, 1) = 1.
        rot(1, 2) = 0.
        rot(1, 3) = 0.
        rot(2, 1) = 0.
        rot(2, 2) = cosa
        rot(2, 3) = -sina
        rot(3, 1) = 0.
        rot(3, 2) = sina
        rot(3, 3) = cosa
      CASE (2) 
        rot(1, 1) = cosa
        rot(1, 2) = -sina
        rot(1, 3) = 0.
        rot(2, 1) = 0.
        rot(2, 2) = 1.
        rot(2, 3) = 0.
        rot(3, 1) = sina
        rot(3, 2) = cosa
        rot(3, 3) = 0.
      CASE (3) 
        rot(1, 1) = cosa
        rot(1, 2) = -sina
        rot(1, 3) = 0.
        rot(2, 1) = sina
        rot(2, 2) = cosa
        rot(2, 3) = 0.
        rot(3, 1) = 0.
        rot(3, 2) = 0.
        rot(3, 3) = 1.
      END SELECT
      IF (lrhs .EQ. 1) THEN
        c1 = 0.
      ELSE
Cextrapolation ordre 0 si variable conservative ou primitive
        c1 = 1.
      END IF
C
      IF (idir .EQ. 1 .OR. idir .EQ. 2) THEN
C
        DO k=ind_loop(5),ind_loop(6)
          DO j=ind_loop(3),ind_loop(4)
            DO i=ind_loop(1),ind_loop(2)
              l = 1 + (i+param_int(0+3)-1) + (j+param_int(0+3)-1)*
     +          param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +          param_int(0+1)
              ropr2d = rot(1, 1)*ropd(l, 2) + rot(1, 2)*ropd(l, 3) + rot
     +          (1, 3)*ropd(l, 4)
              ropr2 = rop(l, 2)*rot(1, 1) + rop(l, 3)*rot(1, 2) + rop(l
     +          , 4)*rot(1, 3)
              ropr3d = rot(2, 1)*ropd(l, 2) + rot(2, 2)*ropd(l, 3) + rot
     +          (2, 3)*ropd(l, 4)
              ropr3 = rop(l, 2)*rot(2, 1) + rop(l, 3)*rot(2, 2) + rop(l
     +          , 4)*rot(2, 3)
              ropr4d = rot(3, 1)*ropd(l, 2) + rot(3, 2)*ropd(l, 3) + rot
     +          (3, 3)*ropd(l, 4)
              ropr4 = rop(l, 2)*rot(3, 1) + rop(l, 3)*rot(3, 2) + rop(l
     +          , 4)*rot(3, 3)
              ropd(l, 1) = c1*ropd(l, 1)
              ropd(l, 2) = c1*ropr2d
              ropd(l, 3) = c1*ropr3d
              ropd(l, 4) = c1*ropr4d
              ropd(l, 5) = c1*ropd(l, 5)
            ENDDO
          ENDDO
        ENDDO
      ELSE IF (idir .EQ. 3 .OR. idir .EQ. 4) THEN
C
C
        DO k=ind_loop(5),ind_loop(6)
          DO j=ind_loop(3),ind_loop(4)
            lij = 1 + (ind_loop(1)+param_int(0+3)-1) + (j+param_int(0+3)
     +        -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +        param_int(0+1)
C
CDIR$ IVDEP
            DO l=lij,lij+ind_loop(2)-ind_loop(1)
C
              ropr2d = rot(1, 1)*ropd(l, 2) + rot(1, 2)*ropd(l, 3) + rot
     +          (1, 3)*ropd(l, 4)
              ropr2 = rop(l, 2)*rot(1, 1) + rop(l, 3)*rot(1, 2) + rop(l
     +          , 4)*rot(1, 3)
              ropr3d = rot(2, 1)*ropd(l, 2) + rot(2, 2)*ropd(l, 3) + rot
     +          (2, 3)*ropd(l, 4)
              ropr3 = rop(l, 2)*rot(2, 1) + rop(l, 3)*rot(2, 2) + rop(l
     +          , 4)*rot(2, 3)
              ropr4d = rot(3, 1)*ropd(l, 2) + rot(3, 2)*ropd(l, 3) + rot
     +          (3, 3)*ropd(l, 4)
              ropr4 = rop(l, 2)*rot(3, 1) + rop(l, 3)*rot(3, 2) + rop(l
     +          , 4)*rot(3, 3)
              ropd(l, 1) = c1*ropd(l, 1)
              ropd(l, 2) = c1*ropr2d
              ropd(l, 3) = c1*ropr3d
              ropd(l, 4) = c1*ropr4d
              ropd(l, 5) = c1*ropd(l, 5)
            ENDDO
          ENDDO
        ENDDO
      ELSE
C
C
C
        DO k=ind_loop(5),ind_loop(6)
          DO j=ind_loop(3),ind_loop(4)
C
            lij = 1 + (ind_loop(1)+param_int(0+3)-1) + (j+param_int(0+3)
     +        -1)*param_int(0) + (k+param_int(0+4)-1)*param_int(0)*
     +        param_int(0+1)
C
CDIR$ IVDEP
            DO l=lij,lij+ind_loop(2)-ind_loop(1)
C
              ropr2d = rot(1, 1)*ropd(l, 2) + rot(1, 2)*ropd(l, 3) + rot
     +          (1, 3)*ropd(l, 4)
              ropr2 = rop(l, 2)*rot(1, 1) + rop(l, 3)*rot(1, 2) + rop(l
     +          , 4)*rot(1, 3)
              ropr3d = rot(2, 1)*ropd(l, 2) + rot(2, 2)*ropd(l, 3) + rot
     +          (2, 3)*ropd(l, 4)
              ropr3 = rop(l, 2)*rot(2, 1) + rop(l, 3)*rot(2, 2) + rop(l
     +          , 4)*rot(2, 3)
              ropr4d = rot(3, 1)*ropd(l, 2) + rot(3, 2)*ropd(l, 3) + rot
     +          (3, 3)*ropd(l, 4)
              ropr4 = rop(l, 2)*rot(3, 1) + rop(l, 3)*rot(3, 2) + rop(l
     +          , 4)*rot(3, 3)
              ropd(l, 1) = c1*ropd(l, 1)
              ropd(l, 2) = c1*ropr2d
              ropd(l, 3) = c1*ropr3d
              ropd(l, 4) = c1*ropr4d
              ropd(l, 5) = c1*ropd(l, 5)
            ENDDO
          ENDDO
        ENDDO
      END IF
      END

