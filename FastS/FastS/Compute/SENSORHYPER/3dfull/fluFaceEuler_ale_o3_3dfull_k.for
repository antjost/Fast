c.....Metrique
 
#include  "FastS/Compute/Normale/normale_3dfull_k.for"

        nm  = l -  inck
        nm2 = l -2*inck
        np  = l +  inck
!        np2  = l +2*inck

!senseur a toute les iter pour le moment a changer 
#include  "FastS/Compute/SENSORHYPER/wiggle_k.for"
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

        vslp = v4                    
#include  "FastS/Compute/Slope/o3_slope_var.for"   
        qm4 = qm                     
        qp4 = qp                     

        vslp = v5
#include  "FastS/Compute/Slope/o3_slope_var.for"
        qm5 = qm
        qp5 = qp

!FLU SENSEUR
!determination etat gauche (rou1) et droit (rou2): ro, roui, roe+p
#include  "FastS/Compute/etat_GD.for"

!determination vitesse normale interface
#include "FastS/Compute/Vit_ent/qn_ale_3dfull_k.for"

       ! modification de vitesse normale par ajout
       ! de stabilisation de type Rhie-Chow
       us  = 0.25*(qn1+qn2)- (1.-wig(l+wigd))*(
     &  c2*sk*(p2-p1)*(opt0*wig( l+ wig_k)+1.-opt0))
       tdu = max(abs(us),c1*sk)*wig( l+ wig_k)
       !Calcul du flux total
       p1p2= (p1+p2)*0.5

#include "FastS/Compute/Vit_ent/fluhyper_ale_3dfull_k.for"
