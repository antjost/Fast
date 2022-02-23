c.....Metrique
 
#include  "FastS/Compute/Normale/normale_2d_i.for"

        nm  = l -  inci
        nm2 = l -2*inci
        np  = l +  inci
!        np2  = l +2*inci

!senseur a toute les iter pour le moment a changer 
#include  "FastS/Compute/SENSORHYPER/wiggle_i.for"
!ducros calcule dans src_term

! pente (qm) a l'interface droite et  (qp) a l'interface gauche
        vslp = v1
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm1 = qm
        qp1 = qp

        vslp = v2
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm2 = qm
        qp2 = qp

        vslp = v3
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm3 = qm
        qp3 = qp


        vslp = v5
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm5 = qm
        qp5 = qp

!FLU SENSEUR
!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_GD_2d.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_ale_2d_i.for"

       ! modification de vitesse normale par ajout
       ! de stabilisation de type Rhie-Chow
       us  = 0.25*(qn1+qn2)- (1.-wig(l+wigd))*(
     &  c2*si*(p2-p1)*(opt0*wig( l+ wig_i)+1.-opt0))
       tdu = max(abs(us),c1*si)*wig( l+ wig_i)
       !Calcul du flux total
       p1p2= (p1+p2)*0.5

#include "FastS/Compute/Vit_ent/fluhyper_ale_2d_i.for"
