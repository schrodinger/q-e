!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
program dynmat
  !--------------------------------------------------------------------
  !! This program:
  !
  !! * reads a dynamical matrix file produced by the phonon code;
  !! * adds the nonanalytical part (if Z* and epsilon are read from file),
  !!   applies the chosen Acoustic Sum Rule (if q=0);
  !! * diagonalise the dynamical matrix; 
  !! * calculates IR and Raman cross sections (if Z* and Raman tensors
  !!   are read from file, respectively);
  !! * writes the results to files, both for inspection and for plotting.
  !
  !! Input data (namelist "input"):
  !
  !! * \(\text{fildyn} [character]: input file containing the dynamical matrix
  !!   (default: fildyn='matdyn')
  !! * \(q(3)) - [real]: calculate LO modes (add nonanalytic terms) along
  !!   the direction q (cartesian axis, default: q=(0,0,0) )
  !! * \(\text{amass}(\text{nt})\) - [real]: mass for atom type nt, amu
  !!   (default: amass is read from file fildyn)
  !! * \(\text{asr}\) - [character]: indicates the type of Acoustic Sum Rule imposed:
  !!    * 'no': no Acoustic Sum Rules imposed (default)
  !!    * 'simple':  previous implementation of the asr used
  !!      (3 translational asr imposed by correction of
  !!      the diagonal elements of the dynamical matrix)
  !!    * 'crystal': 3 translational asr imposed by optimized
  !!      correction of the dyn. matrix (projection).
  !!    * 'one-dim': 3 translational asr + 1 rotational asr
  !!      imposed by optimized correction of the dyn. mat. (the
  !!      rotation axis is the direction of periodicity; it
  !!      will work only if this axis considered is one of
  !!      the cartesian axis).
  !!    * 'zero-dim': 3 translational asr + 3 rotational asr
  !!      imposed by optimized correction of the dyn. mat.
  !!      Note that in certain cases, not all the rotational asr
  !!      can be applied (e.g. if there are only 2 atoms in a
  !!      molecule or if all the atoms are aligned, etc.).
  !!      In these cases the supplementary asr are cancelled
  !!      during the orthonormalization procedure (see below).
  !!      Finally, in all cases except 'no' a simple correction
  !!      on the effective charges is performed (same as in the
  !!      previous implementation).
  !! * \(\text{axis}\) - [integer]: indicates the rotation axis for a 1D system
  !!   (1=Ox, 2=Oy, 3=Oz ; default =3)
  !! * \(\text{lperm}\) - [logical]: TRUE to calculate Gamma-point mode contributions to
  !!   dielectric permittivity tensor (default: lperm=.false.)
  !! * \(\text{lplasma}\) - [logical]: TRUE to calculate Gamma-point mode effective plasma 
  !!   frequencies, automatically triggers lperm = TRUE
  !!   (default: lplasma=.false.)
  !! * \(\text{filout} - [character]: output file containing phonon frequencies and normalized
  !!   phonon displacements (i.e. eigenvectors divided by the
  !!   square root of the mass and then normalized; they are
  !!   not orthogonal). Default: filout='dynmat.out'
  !! * \(\text{fileig}\) - [character]: output file containing phonon frequencies and eigenvectors
  !!   of the dynamical matrix (they are orthogonal). Default: fileig=' '
  !! * \(\text{filmol}\) - [character]: as above, in a format suitable for 'molden'
  !!   (default: filmol='dynmat.mold')
  !! * \(\text{filxsf}\) - [character]: as above, in axsf format suitable for xcrysden
  !!   (default: filxsf='dynmat.axsf')
  !! * \(\text{loto_2d}\) - [logical]: set to TRUE to activate two-dimensional treatment of
  !!   LO-TO splitting.
  !
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_bcast
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE mp_world,    ONLY : world_comm
  USE io_global,   ONLY : ionode, ionode_id, stdout
  USE environment, ONLY : environment_start, environment_end
  USE io_dyn_mat,  ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail
  USE constants,   ONLY : amu_ry
  USE dynamical,  ONLY : dyn, m_loc, ityp, tau, zstar, dchi_dtau  
  USE rigid,       ONLY : dyndiag, nonanal
  !
  implicit none
  !
  integer, parameter :: ntypx = 10
  character(len=256):: fildyn, filout, filmol, filxsf, fileig
  character(len=256) :: prefix, filspm, filvib
  character(len=3) :: atm(ntypx)
  character(len=10) :: asr
  logical :: lread, gamma, loto_2d
  complex(DP), allocatable :: z(:,:)
  real(DP) :: amass(ntypx), amass_(ntypx), eps0(3,3), a0, omega, &
       at(3,3), bg(3,3), q(3), q_(3)
  real(DP), allocatable :: w2(:)
  integer :: nat, na, nt, ntyp, iout, axis, nspin_mag, ios
  real(DP) :: celldm(6)
  logical :: xmldyn, lrigid, lraman, lperm, lplasma, remove_interaction_blocks
  logical, external :: has_xml
  integer :: ibrav, nqs
  integer, allocatable :: itau(:)
  namelist /input/ amass, asr, axis, fildyn, filout, filmol, filxsf, &
                   fileig, lperm, lplasma, q, loto_2d, &
                   remove_interaction_blocks, prefix, filspm, filvib
  !
  ! code is parallel-compatible but not parallel
  !
  CALL mp_startup()
  CALL environment_start('DYNMAT')
  !
  IF (ionode) CALL input_from_file ( )
  !
  asr  = 'no'
  axis = 3
  fildyn='matdyn'
  filout='dynmat.out'
  filmol='dynmat.mold'
  filxsf='dynmat.axsf'
  prefix=' '
  filspm=' '
  filvib=' '
  fileig=' '
  amass(:)=0.0d0
  q(:)=0.0d0
  lperm=.false.
  lplasma=.false.
  loto_2d=.false.
  remove_interaction_blocks = .false. 
  !
  IF (ionode) read (5,input, iostat=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('dynmat', 'reading input namelist', ABS(ios))
  !
  CALL mp_bcast(asr,ionode_id, world_comm)
  CALL mp_bcast(axis,ionode_id, world_comm)
  CALL mp_bcast(amass,ionode_id, world_comm)
  CALL mp_bcast(fildyn,ionode_id, world_comm)
  CALL mp_bcast(filout,ionode_id, world_comm)
  CALL mp_bcast(filmol,ionode_id, world_comm)
  CALL mp_bcast(fileig,ionode_id, world_comm)
  CALL mp_bcast(filxsf,ionode_id, world_comm)
  CALL mp_bcast(q,ionode_id, world_comm)
  CALL mp_bcast(prefix,ionode_id, world_comm)
  CALL mp_bcast(filspm,ionode_id, world_comm)
  CALL mp_bcast(filvib,ionode_id, world_comm)
  !
  IF ( trim( prefix ) /= ' ' ) THEN
     fildyn = trim(prefix) // '.save/' // trim(fildyn)
  END IF
  !
  IF (ionode) inquire(file=fildyn,exist=lread)
  CALL mp_bcast(lread, ionode_id, world_comm)
  IF (lread) THEN
     IF (ionode) WRITE(6,'(/5x,a,a)') 'Reading Dynamical Matrix from file '&
                                     , TRIM(fildyn)
  ELSE
     CALL errore('dynmat', 'File '//TRIM(fildyn)//' not found', 1)
  END IF
  !
  ntyp = ntypx ! avoids spurious out-of-bound errors
  xmldyn=has_xml(fildyn)
  IF (xmldyn) THEN
     CALL read_dyn_mat_param(fildyn,ntyp,nat)
     ALLOCATE (m_loc(3,nat))
     ALLOCATE (tau(3,nat))
     ALLOCATE (ityp(nat))
     ALLOCATE (zstar(3,3,nat))
     ALLOCATE (dchi_dtau(3,3,3,nat) )
     CALL read_dyn_mat_header(ntyp, nat, ibrav, nspin_mag, &
             celldm, at, bg, omega, atm, amass_, tau, ityp, &
             m_loc, nqs, lrigid, eps0, zstar, lraman, dchi_dtau)
     IF (nqs /= 1) CALL errore('dynmat','only q=0 matrix allowed',1)
     a0=celldm(1) ! define alat
     ALLOCATE (dyn(3,3,nat,nat) )
     CALL read_dyn_mat(nat,1,q_,dyn(:,:,:,:))
     CALL read_dyn_mat_tail(nat)
     IF(asr.ne.'no') THEN
         CALL set_asr ( asr, axis, nat, tau, dyn, zstar )
     END IF
     IF (ionode) THEN
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)
        END DO
     END IF
  ELSE
     IF (ionode) THEN
        CALL readmat2 ( fildyn, asr, axis, nat, ntyp, atm, a0, &
                        at, omega, amass_, eps0, q_ )
        DO nt=1, ntyp
           IF (amass(nt) <= 0.0d0) amass(nt)=amass_(nt)/amu_ry
        END DO
     END IF
  ENDIF
  IF (remove_interaction_blocks)  CALL remove_interaction(dyn, nat) 
  !
  IF (ionode) THEN
     !
     ! from now on, execute on a single processor
     !
     gamma = ( abs( q_(1)**2+q_(2)**2+q_(3)**2 ) < 1.0d-8 )
     !
     IF (gamma .and. .not.loto_2d) THEN
        ALLOCATE (itau(nat))
        DO na=1,nat
           itau(na)=na
        END DO
        CALL nonanal ( nat, nat, itau, eps0, q, zstar, omega, dyn )
        DEALLOCATE (itau)
     END IF
     !
     ALLOCATE ( z(3*nat,3*nat), w2(3*nat) )
     CALL dyndiag(nat,ntyp,amass,ityp,dyn,w2,z)
     !
     IF (filout.eq.' ') then
        iout=6
     ELSE
        iout=4
        OPEN (unit=iout,file=filout,status='unknown',form='formatted')
     END IF
     CALL writemodes(nat,q_,w2,z,iout)
     IF(iout .ne. 6) close(unit=iout)
     IF (fileig .ne. ' ') THEN
       OPEN (unit=15,file=TRIM(fileig),status='unknown',form='formatted')
       CALL write_eigenvectors (nat,ntyp,amass,ityp,q_,w2,z,15)
       CLOSE (unit=15)
     ENDIF
     CALL writemolden (filmol, gamma, nat, atm, a0, tau, ityp, w2, z)
     CALL writexsf (filxsf, gamma, nat, atm, a0, at, tau, ityp, z)
     IF (gamma) THEN 
        CALL RamanIR (nat, omega, w2, z, zstar, eps0, dchi_dtau)
        IF (lperm .OR. lplasma) THEN
            CALL polar_mode_permittivity(nat,eps0,z,zstar,w2,omega, &
                                         lplasma)
            IF ( ABS( q(1)**2+q(2)**2+q(3)**2 ) > 1.0d-8 ) &
               WRITE(6,'(5x,a)') 'BEWARE: phonon contribution to &
               & permittivity computed with TO-LO splitting'
        ENDIF
        !
     ENDIF
     !
     CALL dump_ir_raman(prefix, nat, omega, w2, z, zstar, eps0, dchi_dtau)
  ENDIF
  !
  IF (xmldyn) THEN
     DEALLOCATE (m_loc)
     DEALLOCATE (tau)
     DEALLOCATE (ityp)
     DEALLOCATE (zstar)
     DEALLOCATE (dchi_dtau)
     DEALLOCATE (dyn)
  ENDIF
  !
  CALL environment_end('DYNMAT')
  !
  CALL mp_global_end()
  !
  CONTAINS 
    subroutine remove_interaction(d,na_) 
      !! this routine removes from the dynamical matrix the columsn and the rows 
      !! for the atoms with vanishing diagonal blocks (i,j,ia,ia) 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)      :: na_
      COMPLEX(DP), INTENT(INOUT)  :: d(3,3,na_,na_)
      REAL(DP)                :: norm 
      ! 
      INTEGER ia, ipol, jpol  
      COMPLEX(dp) :: z(3,3) 
      DO ia = 1, na_ 
        norm = 0._dp 
        z = d(:,:,ia,ia) 
        DO ipol =1, 3
          norm = norm + REAL(z(ipol,ipol))**2  + AIMAG(z(ipol,ipol))**2   
          DO jpol =ipol+1, 3
             norm = norm + 2._DP * REAL(z(ipol,jpol))**2  + AIMAG(z(ipol,jpol))**2
          END DO 
        END DO 
        IF (norm .lt. 1.e-8_DP ) THEN 
          d(:,:,ia,:) = 0._DP 
          d(:,:,:,ia)  =  0._DP 
        END IF 
      END DO  
   END SUBROUTINE remove_interaction  
   !
   SUBROUTINE dump_ir_raman(prefix, nat, omega, w2, z, zstar, eps0, dchi_dtau)
      !
      USE kinds, ONLY: DP
      USE constants, ONLY : RY_TO_CMM1, amu_ry, eps8
      implicit none
      ! input
      CHARACTER(len=*), INTENT(in) :: prefix
      integer, intent(in) :: nat
      real(DP), intent(in) :: omega, w2(3*nat), zstar(3,3,nat), eps0(3,3), &
           dchi_dtau(3,3,3,nat)
      complex(DP), intent(in) :: z(3*nat,3*nat)
      !
      INTEGER, EXTERNAL :: find_free_unit
      ! local
      integer na, nu, ipol, jpol, lpol, nat3, iunit
      logical noraman
      real(DP) :: freq(3*nat), infrared(3*nat), raman(3, 3, 3*nat), &
         raman_act(3*nat), polar(3), irfac, alpha, beta2
      character(len=256) :: filename
      !
      nat3 = 3 * nat
      !
      !   conversion factor for IR cross sections from
      !   (Ry atomic units * e^2)  to  (Debye/A)^2/amu
      !   1 Ry mass unit = 2 * mass of one electron = 2 amu
      !   1 e = 4.80324x10^(-10) esu = 4.80324 Debye/A
      !     (1 Debye = 10^(-18) esu*cm = 0.2081928 e*A)
      !
      irfac = 4.80324d0**2/2.d0*amu_ry
      !
      noraman = .true.
      do nu = 1, nat3
         !
         freq(nu) = sqrt(abs(w2(nu)))*RY_TO_CMM1
         if (w2(nu).lt.0.0) freq(nu) = -freq(nu)
         !
         polar(:) = 0._DP
         !
         do na=1,nat
            do ipol=1,3
               do jpol=1,3
                  polar(ipol) = polar(ipol) +  &
                       zstar(ipol,jpol,na)*z((na-1)*3+jpol,nu)
               end do
            end do
         end do
         !
         infrared(nu) = 2.d0*(polar(1)**2+polar(2)**2+polar(3)**2)*irfac
         !
         ! Check for nan
         if(infrared(nu) /= infrared(nu)) infrared(nu) = 0.0
         !
         do ipol=1,3
            do jpol=1,3
               raman(ipol,jpol,nu)=0.0d0
               do na=1,nat
                  do lpol=1,3
                     raman(ipol,jpol,nu) = raman(ipol,jpol,nu) + &
                          dchi_dtau(ipol,jpol,lpol,na) * z((na-1)*3+lpol,nu)
                  end do
               end do
               noraman=noraman .and. abs(raman(ipol,jpol,nu)).lt.1.d-12
            end do
         end do
         !   Raman cross sections are in units of bohr^4/(Ry mass unit)
      end do
      !
      filename = TRIM(prefix) // ".vib.spm"
      iunit = find_free_unit()
      OPEN(unit=iunit, file=TRIM(filename), status='unknown', form='formatted')
      CALL writespm(nat, freq, infrared, .FALSE., iunit)
      CLOSE(unit=iunit)
      !
      IF(noraman) THEN
         raman_act(:) = 0._DP
      ELSE
         !
         do nu = 1, nat3
            !
            alpha = (raman(1,1,nu) + raman(2,2,nu) + raman(3,3,nu))/3.d0
            beta2 = ( (raman(1,1,nu) - raman(2,2,nu))**2 + &
                      (raman(1,1,nu) - raman(3,3,nu))**2 + &
                      (raman(2,2,nu) - raman(3,3,nu))**2 + 6.d0 * &
                      (raman(1,2,nu)**2 + raman(1,3,nu)**2 + raman(2,3,nu)**2) )/2.d0
            raman_act(nu) = (45.d0*alpha**2 + 7.0d0*beta2)*amu_ry
         end do
         !
         filename = TRIM(prefix) // ".raman_vib.spm"
         iunit = find_free_unit()
         OPEN(unit=iunit, file=TRIM(filename), status='unknown', form='formatted')
         CALL writespm(nat, freq, raman_act, .TRUE., iunit)
         CLOSE(unit=iunit)
      END IF
      !
      filename = TRIM(prefix) // ".vib"
      iunit = find_free_unit()
      OPEN(unit=iunit, file=TRIM(filename), status='unknown', form='formatted')
      CALL writevib(nat, freq, infrared, raman_act, ntyp, amass, ityp, z, iunit)
      CLOSE(unit=iunit)
      !
      IF ( ANY(zstar .GT. eps8) ) THEN
         ! Check if Born charges are present (occupations = fixed)
         filename = 'born.charges'
         iunit = find_free_unit()
         OPEN(unit=iunit, file=TRIM(filename), status='unknown', form='formatted')
         WRITE(iunit, '(9f24.12)') ((eps0(ipol,jpol), jpol=1,3), ipol=1,3)
         WRITE(iunit, '(*(f24.12))') &
            (((zstar(ipol,jpol,na), jpol=1,3), ipol=1,3), na=1,nat)
         CLOSE(unit=iunit)
      END IF
      !
   END SUBROUTINE dump_ir_raman
   !
end program dynmat
