c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 aout 2013) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine fluausm_select_d(ndom, nitcfg, ithread, nptpsi, 
     &                        param_int, param_real,
     &                        ind_dm, ind_loop, ijkv_thread,ijkv_sdm,
     &                        synchro_send_sock, synchro_send_th, 
     &                        synchro_receive_sock,synchro_receive_th,
     &                        ibloc , jbloc , kbloc ,
     &                        icache, jcache, kcache,
     &                        psi, wig, stat_wig,
     &                        rop, ropd, drodm, drodmd,
     &                        ti, ti_df,tj,tj_df,tk,tk_df,vol,vol_df,
     &                        venti, ventj, ventk, xmut, xmutd)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Appel du calcul des flux explicites
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    tijk     : vecteur normale aux facettes des mailles
c_I    ventijk     : vitesses d entrainement aux facettes preced.
c_I    qm,qp    : etats droit et gauche aux interfaces d une maille
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    drodm    : increment de la solution
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, nitcfg,ithread, nptpsi,
     & icache, jcache, kcache,ibloc, jbloc, kbloc,
     & ijkv_thread(3), ijkv_sdm(3), ind_loop(6),ind_dm(6),
     & synchro_send_sock(3),synchro_send_th(3),
     & synchro_receive_sock(3), synchro_receive_th(3), param_int(0:*)


      REAL_E rop(*),xmut(*),drodm(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),stat_wig(*),
     & ti_df(*),tj_df(*),tk_df(*),vol_df(*)

      REAL_E drodmd(*), ropd(*), xmutd(*)

      REAL_E param_real(0:*)
      REAL_E psi(nptpsi)

C Var loc
      INTEGER_E option

      if(ind_loop(1).gt.ind_loop(2)) return 
      if(ind_loop(3).gt.ind_loop(4)) return 
      if(ind_loop(5).gt.ind_loop(6)) return

      option =  1000*param_int(LALE)
     &        +  100*param_int(SLOPE)
     &        +   10*param_int(IFLOW)
     &        +      param_int(ITYPZONE)


       IF  (option.eq.223) THEN
                                               
           call fluausm_lamin_o3_2d_d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
       ELSEIF (option.eq.233) THEN
                                               
           call fluausm_SA_o3_2d_d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut, xmutd)
                                               
       ELSEIF (option.eq.213) THEN
                                               
           call fluausm_euler_o3_2d_d(ndom, ithread,
     &                 param_int, param_real,
     &                 ind_dm, ind_loop, ijkv_thread, ijkv_sdm,
     &                 synchro_send_sock, synchro_send_th,
     &                 synchro_receive_sock, synchro_receive_th,
     &                 ibloc , jbloc , kbloc ,
     &                 icache, jcache, kcache,
     &                 rop, ropd, drodm, drodmd, wig,
     &                 venti, ventj, ventk,
     &                 ti, tj, tk, vol, xmut)
                                               
      ELSE
         write(*,*) ' option = ' , option 
            write(*,*)'Unknown flux options'
           call error('fluselect_d$',70,1)

      ENDIF
 
      end
