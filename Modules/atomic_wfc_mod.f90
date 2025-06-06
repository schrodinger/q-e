!
! Copyright (C) 2023-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE atomic_wfc_acc( xk, npw, igk_k, nat, nsp, ityp, tau, &
     noncolin, domag, angle1, angle2, starting_spin_angle, &
     npwx, npol, natomwfc, wfcatom )
  !-----------------------------------------------------------------------
  !! This routine computes the superposition of atomic wavefunctions
  !! See below for input variables, output on wfcatom (ACC variable)
  !! Computation is performed on GPU if available
  !! Can be called by CP as well (does not use PW-specific modules)
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE cell_base,        ONLY : omega, tpiba
  USE gvect,            ONLY : mill, eigts1, eigts2, eigts3, g
  USE uspp_param,       ONLY : upf, nwfcm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: xk(3)
  !! k-point
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: nsp
  !! number of types of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! indices of the type of atom  for each atom
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (in units of alat)
  INTEGER, INTENT(IN) :: npw
  !! number of plane waves
  INTEGER, INTENT(IN) :: igk_k(npw)
  !! index of G in the k+G list
  LOGICAL, INTENT(IN) ::  noncolin
  !! true if calculation noncolinear
  LOGICAL, INTENT(IN) :: domag
  !! true if nonzero noncolinear magnetization
  LOGICAL, INTENT(IN) :: starting_spin_angle
  !! true if initial spin direction is set
  REAL(DP), INTENT(IN) :: angle1(nsp)
  !! angle theta of initial spin direction
  REAL(DP), INTENT(IN) ::  angle2(nsp)
  !! angle phi of initial spin direction
  INTEGER, INTENT(IN) :: npol
  !! npol = 2 for noncolinear calculations
  INTEGER, INTENT(IN) :: npwx
  !! max number of plane waves
  INTEGER, INTENT(IN) :: natomwfc
  !! number of atomic wavefunctions
  COMPLEX(DP), INTENT(OUT) :: wfcatom(npwx,npol,natomwfc)
  !! Superposition of atomic wavefunctions
  !
  ! ... local variables
  !
  INTEGER :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3
  COMPLEX(DP) :: kphase
  REAL(DP)    :: arg, px, ux, vx, wx
  !
  REAL(DP) :: xk1, xk2, xk3, qgr
  REAL(DP), ALLOCATABLE :: chiq(:,:,:), qg(:)
  REAL(DP), ALLOCATABLE :: ylm(:,:), gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  !
  !
  ! calculate max angular momentum required in wavefunctions
  lmax_wfc = 0
  DO nt = 1, nsp
     lmax_wfc = MAX( lmax_wfc, MAXVAL( upf(nt)%lchi(1:upf(nt)%nwfc) ) )
  END DO
  !
  ALLOCATE( ylm(npw,(lmax_wfc+1)**2), chiq(npw,nwfcm,nsp)) 
  ALLOCATE( qg(npw), gk(3,npw), sk(npw) )
  !$acc data create (ylm, chiq, gk, qg, sk) &
  !$acc      present(g, igk_k, eigts1, eigts2, eigts3, mill, wfcatom)
  !
  xk1 = xk(1)
  xk2 = xk(2)
  xk3 = xk(3)
  !
  !$acc parallel loop
  DO ig = 1, npw
     iig = igk_k(ig)
     gk(1,ig) = xk1 + g(1,iig)
     gk(2,ig) = xk2 + g(2,iig)
     gk(3,ig) = xk3 + g(3,iig)
     qg(ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
  END DO
  !
  !  ylm = spherical harmonics
  !
  CALL ylmr2( (lmax_wfc+1)**2, npw, gk, qg, ylm )
  !
  ! set now q=|k+G| in atomic units
  !
  !$acc parallel loop
  DO ig = 1, npw
     qg(ig) = SQRT( qg(ig) )*tpiba
  END DO
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  CALL interp_atwfc ( npw, qg, nwfcm, chiq )
  !
  !$acc kernels
  wfcatom(:,:,:) = (0.0_dp, 0.0_dp)
  !$acc end kernels
  !
  n_starting_wfc = 0
  !
  DO na = 1, nat
     arg = (xk1*tau(1,na) + xk2*tau(2,na) + xk3*tau(3,na)) * tpi
     kphase = CMPLX( COS(arg), - SIN(arg) ,KIND=DP)
     !
     !     sk is the structure factor
     !
     !$acc parallel loop
     DO ig = 1, npw
        iig = igk_k(ig)
        sk(ig) = kphase * eigts1(mill(1,iig),na) * &
                          eigts2(mill(2,iig),na) * &
                          eigts3(mill(3,iig),na)
     END DO
     !
     nt = ityp(na)
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           !
           !  the factor i^l MUST BE PRESENT in order to produce
           !  wavefunctions for k=0 that are real in real space
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 IF (starting_spin_angle.OR..NOT.domag) THEN
                    CALL atomic_wfc_so ( npw, npwx, npol, natomwfc, nsp, nt, &
                         nb, lmax_wfc, ylm, chiq, sk, n_starting_wfc, wfcatom )
                 ELSE
                    CALL atomic_wfc_so_mag ( npw, npwx, npol, natomwfc, nsp, &
                         nt, nb, angle1, angle2, lmax_wfc, ylm, chiq, sk, &
                         n_starting_wfc, wfcatom )
                 END IF
                 !
              ELSE
                 !
                 CALL atomic_wfc_nc ( npw, npwx, npol, natomwfc, nsp, &
                         nt, nb, angle1, angle2, lmax_wfc, ylm, chiq, sk, &
                         n_starting_wfc, wfcatom )
                 !
              END IF
              !
           ELSE
              !
              CALL atomic_wfc_lsda  ( npw, npwx, npol, natomwfc, &
                         nsp, nt, nb, lmax_wfc, ylm, chiq, sk, &
                         n_starting_wfc, wfcatom )
              !
           END IF
           !
        END IF
        !
     END DO
     !
  END DO

  IF ( n_starting_wfc /= natomwfc) call errore ('atomic_wfc', &
       'internal error: some wfcs were lost ', 1 )

  !$acc end data
  DEALLOCATE( sk, gk, qg, chiq, ylm ) 
  
  RETURN

END SUBROUTINE atomic_wfc_acc
!
!----------------------------------------------------------------
SUBROUTINE atomic_wfc_so( npw, npwx, npol, natomwfc, nsp, nt, &
     nb, lmax_wfc, ylm, chiq, sk, n_starting_wfc, wfcatom )
   !------------------------------------------------------------
   !! Spin-orbit case, no magnetization
   !
   USE kinds,            ONLY : DP
   USE upf_spinorb,      ONLY : rot_ylm, lmaxx
   USE uspp_param,       ONLY : upf, nwfcm
   !
   IMPLICIT NONE
   INTEGER,  INTENT(IN)  :: nsp, nt, nb, natomwfc, npw, npwx, npol, lmax_wfc
   REAL(DP), INTENT(IN) :: chiq(npw,nwfcm,nsp)
   REAL(DP), INTENT(IN) :: ylm(npw,(lmax_wfc+1)**2)
   COMPLEX(DP), INTENT(IN) :: sk(npw)
   INTEGER, INTENT(INOUT) :: n_starting_wfc
   COMPLEX(DP), INTENT(INOUT) :: wfcatom(npwx,npol,natomwfc)
   !
   REAL(DP) :: fact(2), fact_is, j
   COMPLEX(DP) :: rot_ylm_in1, lphase
   REAL(DP), EXTERNAL :: spinor
   INTEGER,  EXTERNAL :: sph_ind
   INTEGER :: l, ind, ind1, n1, is, m, ig
   !
   j = upf(nt)%jchi(nb)
   l = upf(nt)%lchi(nb)
   lphase = (0.d0,1.d0)**l
   !
   DO m = -l-1, l
      fact(1) = spinor(l,j,m,1)
      fact(2) = spinor(l,j,m,2)
      IF ( ABS(fact(1)) > 1.d-8 .OR. ABS(fact(2)) > 1.d-8 ) THEN
         n_starting_wfc = n_starting_wfc + 1
         IF (n_starting_wfc > natomwfc) CALL errore &
              ('atomic_wfc_so', 'internal error: too many wfcs', 1)
         !     
         DO is = 1, 2
            fact_is = fact(is)
            IF (ABS(fact(is)) > 1.d-8) THEN
               ind = lmaxx + 1 + sph_ind(l,j,m,is)
               DO n1 = 1, 2*l+1
                  ind1 = l**2+n1
                  rot_ylm_in1 = rot_ylm(ind,n1)
                  IF (ABS(rot_ylm_in1) > 1.d-8) THEN
                     !$acc parallel loop
                     DO ig = 1, npw
                        wfcatom(ig,is,n_starting_wfc) = &
                             wfcatom(ig,is,n_starting_wfc) + &
                             lphase * rot_ylm_in1 * sk(ig) * &
                                CMPLX(ylm(ig,ind1)*fact_is* &
                                chiq(ig,nb,nt), KIND=DP)
                     END DO
                  ENDIF
               ENDDO
            END IF
         END DO
         !
      END IF
   END DO
   !
  END SUBROUTINE atomic_wfc_so
  ! 
  SUBROUTINE atomic_wfc_so_mag( npw, npwx, npol, natomwfc, nsp, nt, &
       nb, angle1, angle2, lmax_wfc, ylm, chiq, sk, &
       n_starting_wfc, wfcatom )
   !------------------------------------------------------------
   !
   !! Spin-orbit case, magnetization along "angle1" and "angle2"
   !! In the magnetic case we always assume that magnetism is much larger
   !! than spin-orbit and average the wavefunctions at l+1/2 and l-1/2
   !! filling then the up and down spinors with the average wavefunctions,
   !! according to the direction of the magnetization, following what is
   !! done in the noncollinear case.
   !
   USE kinds,            ONLY : DP
   USE constants,        ONLY : pi
   USE uspp_param,       ONLY : upf, nwfcm
   !
   IMPLICIT NONE
   INTEGER,  INTENT(IN)  :: nsp, nt, nb, natomwfc, npw, npwx, npol, lmax_wfc
   REAL(DP), INTENT(IN) :: chiq(npw,nwfcm,nsp)
   REAL(DP), INTENT(IN) :: ylm(npw,(lmax_wfc+1)**2)
   COMPLEX(DP), INTENT(IN) :: sk(npw)
   REAL(DP), INTENT(IN) :: angle1(*)
   !! angle theta of initial spin direction
   REAL(DP), INTENT(IN) ::  angle2(*)
   !! angle phi of initial spin direction
   INTEGER, INTENT(INOUT) :: n_starting_wfc
   COMPLEX(DP), INTENT(INOUT) :: wfcatom(npwx,npol,natomwfc)
   !
   REAL(DP) :: alpha, gamman, j
   COMPLEX(DP) :: fup, fdown, aux, lphase
   INTEGER :: nc, ib, ig, l, m, lm
   !
   j = upf(nt)%jchi(nb)
   l = upf(nt)%lchi(nb)
   lphase = (0.d0,1.d0)**l
   !
   !
   !  This routine creates two functions only in the case j=l+1/2 or exit in the
   !  other case 
   !    
   IF (ABS(j-l+0.5_DP)<1.d-4) RETURN
   !
   !  Find the functions j=l-1/2
   !
   nc = nb
   IF (l > 0)  THEN
      DO ib = 1, upf(nt)%nwfc
         IF ((upf(nt)%lchi(ib) == l).AND. &
              (ABS(upf(nt)%jchi(ib)-l+0.5_DP)<1.d-4)) THEN
            nc = ib
            EXIT
         ENDIF
      ENDDO
   END IF
   !
   !  and construct the starting wavefunctions as in the noncollinear case.
   !
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2*l+1
      lm = l**2+m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc ) CALL errore &
            ('atomic_wfc_so_mag', 'internal error: too many wfcs', 1)
      !
      !$acc parallel loop
      DO ig = 1, npw
         !
         !  Average the two functions
         !
         aux = lphase * sk(ig)* CMPLX( ylm(ig,lm) * (chiq(ig,nb,nt)*DBLE(l+1)+&
              chiq(ig,nc,nt)*l)/DBLE(2*l+1), KIND=DP )
         !
         ! now, rotate wfc as needed
         ! first : rotation with angle alpha around (OX)
         !
         fup = CMPLX(COS(0.5d0*alpha), KIND=DP)*aux
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha), KIND=DP)*aux
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
              -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)), KIND=DP)*aux
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)), KIND=DP)*aux
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0 *gamman), KIND=DP))*fup
         wfcatom(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
      END DO
   END DO
   !
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
  END SUBROUTINE atomic_wfc_so_mag
  !
  !
  SUBROUTINE atomic_wfc_nc( npw, npwx, npol, natomwfc, nsp,  nt, &
       nb, angle1, angle2, lmax_wfc, ylm, chiq, sk, &
       n_starting_wfc, wfcatom )
   !
   !! noncolinear case, magnetization along "angle1" and "angle2"
   !
   USE kinds,            ONLY : DP
   USE constants,        ONLY : pi 
   USE uspp_param,       ONLY : upf, nwfcm
   !
   IMPLICIT NONE
   INTEGER,  INTENT(IN)  :: nsp, nt, nb, natomwfc, npw, npwx, npol, lmax_wfc
   REAL(DP), INTENT(IN) :: chiq(npw,nwfcm,nsp)
   REAL(DP), INTENT(IN) :: ylm(npw,(lmax_wfc+1)**2)
   COMPLEX(DP), INTENT(IN) :: sk(npw)
   REAL(DP), INTENT(IN) :: angle1(*)
   !! angle theta of initial spin direction
   REAL(DP), INTENT(IN) ::  angle2(*)
   !! angle phi of initial spin direction
   INTEGER, INTENT(INOUT) :: n_starting_wfc
   COMPLEX(DP), INTENT(INOUT) :: wfcatom(npwx,npol,natomwfc)
   !
   REAL(DP) :: alpha, gamman
   COMPLEX(DP) :: fup, fdown, aux, lphase
   INTEGER :: m, lm, ig, l  
   !
   l = upf(nt)%lchi(nb)
   lphase = (0.d0,1.d0)**l
   alpha = angle1(nt)
   gamman = - angle2(nt) + 0.5d0*pi
   !
   DO m = 1, 2*l+1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc + 2*l+1 > natomwfc) CALL errore &
            ('atomic_wfc_nc', 'internal error: too many wfcs', 1)
      !$acc parallel loop
      DO ig = 1, npw
         aux = lphase*sk(ig)*CMPLX(ylm(ig,lm)*chiq(ig,nb,nt), KIND=DP)
         !
         ! now, rotate wfc as needed
         ! first : rotation with angle alpha around (OX)
         !
         fup = CMPLX(COS(0.5d0*alpha), KIND=DP)*aux
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*alpha), KIND=DP)*aux
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         wfcatom(ig,1,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom(ig,2,n_starting_wfc) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                        -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
         !
         ! second: rotation with angle gamma around (OZ)
         !
         ! Now, build the orthogonal wfc
         ! first rotation with angle (alpha+pi) around (OX)
         !
         fup = CMPLX(COS(0.5d0*(alpha+pi)), KIND=DP)*aux
         fdown = (0.d0,1.d0)*CMPLX(SIN(0.5d0*(alpha+pi)), KIND=DP)*aux
         !
         ! second, rotation with angle gamma around (OZ)
         !
         wfcatom(ig,1,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  +(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fup
         wfcatom(ig,2,n_starting_wfc+2*l+1) = (CMPLX(COS(0.5d0*gamman), KIND=DP) &
                  -(0.d0,1.d0)*CMPLX(SIN(0.5d0*gamman), KIND=DP))*fdown
      END DO
   END DO
   n_starting_wfc = n_starting_wfc + 2*l+1
   !
  END SUBROUTINE atomic_wfc_nc
  !
  !
  SUBROUTINE atomic_wfc_lsda( npw, npwx, npol, natomwfc, nsp, nt, &
       nb, lmax_wfc, ylm, chiq, sk, n_starting_wfc, wfcatom )
   !
   !! LSDA or nonmagnetic case
   !
   USE kinds,            ONLY : DP
   USE uspp_param,       ONLY : upf, nwfcm
   !
   IMPLICIT NONE
   INTEGER,  INTENT(IN)  :: nsp, nt, nb, natomwfc, npw, npwx, npol, lmax_wfc
   REAL(DP), INTENT(IN) :: chiq(npw,nwfcm,nsp)
   REAL(DP), INTENT(IN) :: ylm(npw,(lmax_wfc+1)**2)
   COMPLEX(DP), INTENT(IN) :: sk(npw)
   INTEGER, INTENT(INOUT) :: n_starting_wfc
   COMPLEX(DP), INTENT(INOUT) :: wfcatom(npwx,npol,natomwfc)
   !
   COMPLEX(DP) :: lphase
   INTEGER :: m, lm, ig, l
   !
   l = upf(nt)%lchi(nb)
   lphase = (0.d0,1.d0)**l
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
      IF ( n_starting_wfc > natomwfc) CALL errore &
         ('atomic_wfc_lsda', 'internal error: too many wfcs', 1)
      !
      !$acc parallel loop
      DO ig = 1, npw
         wfcatom(ig,1,n_starting_wfc) = lphase * &
            sk(ig) * CMPLX(ylm(ig,lm) * chiq(ig,nb,nt), KIND=DP)
      ENDDO
      !
   END DO
   !
  END SUBROUTINE atomic_wfc_lsda
