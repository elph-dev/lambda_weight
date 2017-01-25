!
! Copyright (C) 2002-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
PROGRAM lambda_weight
  !
  ! read files 'filelph' produced by phonon (one for each q-point)
  ! sum over q-points to produce the electron-phonon coefficients:
  ! lambda using weighted-average method
  ! T_c using Allen-Dynes formula
  ! see T.Koretsune and R.Arita, arXiv:1610.09441
  !
  ! INPUT from standard input:
  ! &input
  !   ...
  ! /
  ! nq
  ! q(1,1)   q(2,1)   q(3,1)   wq(1)
  !   ...      ...      ...     ...
  ! q(1,nq)  q(2,nq)  q(3,nq)  wq(nq)
  !
  !
  ! Input cards: namelist &lambda
  !
  !  DOStetra        tetrahedron DOS including spin degrees of freedom (states/Ry) 
  !
  !  filelph         prefix of el-ph files produced by ph.x
  !
  !  flfrq           input file fro phonon frequency produced by matdyn.x
  !
  !  read_flfrq      if .true. phonon frequecies calculated by matdyn.x are used.
  !                  otherwise those calculated by ph.x are used.
  !
  !  flambda         output file for lambda
  !
  !  mustar          mu^* for Allen-Dynes formula
  !
  !  emax_a2F        a2F is plotted from 0 to emax_a2F (cm-1)
  !  degauss_a2F     smearing width for sum over q (cm-1) default 1
  !  ngauss_a2F      default 0
  !
  USE constants, ONLY : pi, ry_to_cmm1, ry_to_ghz, rytoev, ry_to_kelvin
  USE kinds, ONLY : DP
  USE environment, ONLY : environment_start, environment_end
  IMPLICIT NONE
  INTEGER:: ios, i, iq, nq, mu, nu, nmodes, nmodes1, ng, nsig, nsig1, nbnd, nks
  INTEGER:: ngauss, ngauss1, ngauss_a2F, nemax_a2F
  INTEGER:: iuelph, iufrq
  LOGICAL:: read_flfrq, use_tetra
  REAL(DP), PARAMETER:: eps = 2.d0  ! in cm^-1
  REAL(DP):: DOStetra, mustar, q_in(3), degauss1, cmm1_to_kelvin
  REAL(DP):: lambda_a2F, omegalog_a2F, emax_a2F, degauss_a2F, de_a2F
  REAL(DP):: e, phase_space, gammaq, dosef1, ef1, x, Tc
  REAL(DP), ALLOCATABLE:: q(:,:), wq(:), q_matdyn(:,:)
  REAL(DP), ALLOCATABLE:: w2(:), omega(:,:), omega_matdyn(:,:), omegalog(:)
  REAL(DP), ALLOCATABLE:: degauss(:), phase_space_sum(:)
  REAL(DP), ALLOCATABLE:: lambda(:), lambdaq(:,:,:), dosef(:), ef(:), alpha2F(:,:)
  REAL(DP), ALLOCATABLE:: sum_alpha2F(:)
  CHARACTER(LEN=256) :: filelph, filelph_n, flfrq, fmt_str, flambda
  !
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(LEN=6) :: int_to_char
  REAL(DP), EXTERNAL :: w0gauss
  !
  namelist /input/ DOStetra, filelph, flfrq, read_flfrq, mustar, flambda
  namelist /plot/ nbnd, nks
  !
  cmm1_to_kelvin = ry_to_kelvin/ry_to_cmm1
  !
  CALL environment_start('LAMBDA_WGT')
  !
  CALL input_from_file ( )
  !
  ! set namelist default
  !
  DOStetra = -1
  filelph = ''
  flfrq = ''
  flambda = 'lambda_weight.dat'
  read_flfrq = .false.
  mustar = 0
  emax_a2F = -1
  degauss_a2F = 5.d0
  ngauss_a2F = 0
  !
  READ ( 5, input, IOSTAT =ios )
  CALL errore('lambda_weight','error reading input namelist', abs(ios))
  !
  use_tetra = .true.
  IF (DOStetra < 0) THEN
     use_tetra = .false.
     !CALL errore('lambda_weight', 'DOStetra is not specified', 1)
  END IF
  DOStetra = DOStetra*rytoev/2
  !
  READ(5,*) nq
  ALLOCATE ( q(3,nq), wq(nq) )
  DO iq = 1, nq
     READ(5,*) (q(i,iq), i=1,3), wq(iq)
  END DO
  !
  !  normalization of weight
  !
  wq(1:nq) = wq(1:nq)/sum(wq(1:nq))
  !
  !  read flfrq
  !
  IF (read_flfrq) THEN
     iufrq = find_free_unit()
     OPEN( UNIT = iufrq, FILE = TRIM(flfrq), IOSTAT = ios)
     CALL errore ('lambda_weight', 'opening file '//flfrq, ABS (ios) )
     !
     WRITE(6,*) ' reading ', TRIM(flfrq)
     READ(iufrq, plot, IOSTAT = ios)
     IF(nq /= nks) &
        CALL errore('lambda_weight', 'inconsistent nq and nks', nks)
  ALLOCATE( q_matdyn(3,nks), omega_matdyn(nbnd, nks) )
     DO iq=1, nks
        READ(iufrq, *) q_matdyn(1:3, iq)
        READ(iufrq, *) (omega_matdyn(mu,iq), mu=1, nbnd)
     END DO
     CLOSE(iufrq)
  END IF
  !
  DO iq = 1, nq
     !
     ! read filelph
     !
     filelph_n=TRIM(filelph)//TRIM(int_to_char(iq))//'.elph.'//TRIM(int_to_char(iq))
     iuelph = find_free_unit()
     OPEN (UNIT = iuelph, FILE = TRIM(filelph_n), IOSTAT = ios)
     CALL errore ('lambda_weight', 'opening file '//filelph_n, ABS (ios) )
     WRITE(6,*) ' reading ', TRIM(filelph_n)
     !
     READ(iuelph,*) q_in(1:3), nsig, nmodes
     IF( sum(abs(q_in(1:3) - q(1:3,iq))) > 1d-3 ) &
        CALL errore('lambda_weight', 'inconsistent q', iq)
     IF(iq == 1) THEN
        ALLOCATE ( w2(nmodes), lambdaq(nmodes,nsig,nq), omega(nmodes, nq) )
        ALLOCATE ( degauss(nsig), dosef(nsig), ef(nsig), phase_space_sum(nsig) )
        w2 = 0
        lambdaq = 0
        degauss = 0
        dosef = 0
        ef = 0
        phase_space_sum = 0
        omega = 0
        nmodes1 = nmodes
        nsig1 = nsig
     ELSE
        IF(nsig1 /= nsig .or. nmodes1 /= nmodes) &
           CALL errore('lambda_weight', 'inconsistent nsig or nmodes', iq)
     END IF
     !
     READ(iuelph,*) (w2(nu), nu=1, nmodes)
     !
     ! set omega in unit of cm-1
     !
     IF(read_flfrq) THEN
        IF(nmodes /= nbnd) &
           CALL errore('lambda_weight', 'inconsistent nmodes and nbnd', nbnd)
        omega(1:nmodes,iq) = omega_matdyn(1:nmodes,iq)
        write(6,*) '|delta omega| =', sum(abs(sqrt(w2(1:nmodes)) * ry_to_cmm1 - omega_matdyn(1:nmodes,iq)))
     ELSE
        omega(1:nmodes,iq) = sqrt(w2(1:nmodes)) * ry_to_cmm1
     END IF
     !
     DO ng=1, nsig
        READ(iuelph,'(26x,f7.3,12x,i4)') degauss1, ngauss1
        READ(iuelph,'(10x,f10.6,32x,f10.6)') dosef1, ef1
        READ(iuelph,'(25x,f20.10)') phase_space
        phase_space_sum(ng) = phase_space_sum(ng) + phase_space * wq(iq)
        DO mu=1, nmodes
           READ(iuelph,'(12x,i5,2x,f15.8,9x,f12.4)', iostat=ios) nu, lambdaq(mu,ng,iq), gammaq
           if(ios > 0) then
             if(iq == 1) then
               write(6,*) 'read lambda at q=0 mode : failed'
             else
               write(6,*) 'read lambda : failed'
               stop
             end if
           end if
           ! rescaling using matdyn phonon frequency (if needed)
           !IF(omega(nu,iq) > eps) &
              !lambdaq(mu,ng,iq) = lambdaq(mu,ng,iq) * w2(nu) / (omega(nu,iq)/ry_to_cmm1)**2
        END DO
        IF (iq == 1) THEN
           degauss(ng) = degauss1
           ngauss = ngauss1
           dosef(ng) = dosef1
           ef(ng) = ef1
        ELSE
           IF(degauss(ng) /= degauss1 .or. ngauss /= ngauss1) &
              CALL errore('lambda_weight', 'inconsistent gaussian broadening', iq)
           IF(abs(dosef(ng) - dosef1) > 1d-1 .or. abs(ef(ng) - ef1) > 1d-1) &
              CALL errore('lambda_weight', 'inconsistent DOS(EF) or EF', iq)
        END IF
     END DO ! ng
     CLOSE(UNIT=iuelph)
  END DO  ! iq
  !
  ! rescale using tetrahedron DOS
  !
  IF(use_tetra) THEN
     WRITE(6,*) "# rescaling a2F"
     WRITE(6,*) "# degauss,  DOS(smearing) => DOS(tetrahedron) (states/Ry/spin)"
     DO ng=1, nsig
        !write(6,'(f9.5, 2f15.8)') degauss(ng), dosef(ng), DOStetra
        write(6,'(f9.5, 3f15.8)') degauss(ng), 2*dosef(ng)**2, phase_space_sum(ng)
        ! phase_space_sum ~ 2 dosef(ng)^2
        lambdaq(:,ng,:) = lambdaq(:,ng,:)/phase_space_sum(ng) * dosef(ng) * DOStetra * 2
     END DO
  END IF
  !
  ! calculate lambda, alpha2F
  !
  de_a2F = degauss_a2F/10
  IF(emax_a2F < 0) emax_a2F = maxval(omega(1:nmodes,1:nq))
  nemax_a2F = int(emax_a2F/de_a2F+20)
  !
  ALLOCATE(lambda(nsig), omegalog(nsig), alpha2F(nemax_a2F, nsig), sum_alpha2F(nsig))
  lambda = 0
  omegalog = 0
  alpha2F = 0
  DO ng=1, nsig
     DO iq=1, nq
        !write(6,*) iq, ng
        DO mu=1, nmodes
           ! lambdaq already contains a factor of 2 (spin)
           IF(omega(mu,iq) > eps) THEN
              lambda(ng) = lambda(ng) + wq(iq) * lambdaq(mu,ng,iq)
              omegalog(ng) = omegalog(ng) + wq(iq) * lambdaq(mu,ng,iq) * log(omega(mu,iq))
              DO i=1,nemax_a2F
                 e = (i-1)*de_a2F
                 alpha2F(i,ng) = alpha2F(i,ng) + &
                      wq(iq) * lambdaq(mu,ng,iq) * omega(mu,iq) * 0.5d0 * &
                      w0gauss((e-omega(mu,iq))/degauss_a2F,ngauss_a2F)/degauss_a2F
              END DO
           END IF
        END DO
     END DO ! iq
     omegalog(ng) = exp(omegalog(ng)/lambda(ng))*cmm1_to_kelvin
  END DO  ! ng
  !
  ! output lambda_weight.dat
  !
  OPEN(UNIT=iuelph,FILE=TRIM(flambda),STATUS='unknown')
  write(iuelph, '("# degauss(Ry),    lambda, lambda_a2F,     w_log(K), w_log_a2F(K),      DOS(EF)")')
  DO ng=1, nsig
     lambda_a2F = 0.d0
     omegalog_a2F = 0.d0
     do i=2, nemax_a2F
        e = (i-1)*de_a2F
        lambda_a2F = lambda_a2F + alpha2F(i,ng)/e
        omegalog_a2F = omegalog_a2F + alpha2F(i,ng)*log(e)/e
     END DO
     lambda_a2F = lambda_a2F * 2 * de_a2F
     omegalog_a2F = exp(omegalog_a2F * 2 * de_a2F / lambda_a2F) * cmm1_to_kelvin
     WRITE(iuelph,'(3f12.6, 3f14.6)') &
        degauss(ng), lambda(ng), lambda_a2F, omegalog(ng), omegalog_a2F, dosef(ng)
  END DO
  CLOSE(UNIT=iuelph)
  !
  ! output lambda, omega_log, Tc
  !
  WRITE(6,*)
  WRITE(6,'("lambda", 8x, "omega_log(K)", 10x, "T_c(K) at mu*=", f8.4)') mustar
  DO ng=1, nsig
     WRITE(6,'(f10.5,5x,f9.3,10x,f9.3)')  &
        lambda(ng), omegalog(ng), Ad_Tc(lambda(ng), omegalog(ng), mustar)
  END DO
  !
  ! output alpha2F.dat
  !
  OPEN(UNIT=iuelph,FILE='alpha2F.dat',STATUS='unknown')
  WRITE(fmt_str, '(i3,"f10.5")') nsig
  WRITE(iuelph, '("# degauss(Ry):",' // TRIM(fmt_str) // ')') (degauss(ng),ng=1,nsig)
  WRITE(iuelph, '("# energy (cm-1),  a2F(ng), ng = 1,", i5)') nsig
  DO i=1,nemax_a2F
     e=(i-1)*de_a2F
     WRITE(iuelph,'(f16.6,' // TRIM(fmt_str) // ')') e,(alpha2F(i,ng),ng=1,nsig)
  END DO
  CLOSE(UNIT=iuelph)
  !
  ! output sum_alpha2F.dat
  !
  OPEN(UNIT=iuelph,FILE='sum_alpha2F.dat',STATUS='unknown')
  WRITE(fmt_str, '(i3,"f10.5")') nsig
  WRITE(iuelph, '("# degauss(Ry):",' // TRIM(fmt_str) // ')') (degauss(ng),ng=1,nsig)
  WRITE(iuelph, '("# energy (cm-1),  sum a2F(ng), ng = 1,", i5)') nsig
  sum_alpha2F(1:nsig) = 0
  DO i=2,nemax_a2F
     e=(i-1)*de_a2F
     sum_alpha2F(1:nsig) = sum_alpha2F(1:nsig) + alpha2F(i,1:nsig)/e * 2 * de_a2F
     WRITE(iuelph,'(f16.6,' // TRIM(fmt_str) // ')') e,(sum_alpha2F(ng),ng=1,nsig)
  END DO
  CLOSE(UNIT=iuelph)
  !
  CALL environment_end('LAMBDA_WGT')
  !
  STOP
  !
  !
CONTAINS
  !
  ! Allen-Dynes formula
  !
  function AD_Tc(lambda, omegalog, mustar)
    real(DP):: AD_Tc, lambda, omegalog, mustar
    AD_Tc = omegalog/1.2d0*exp(-1.04d0*(1+lambda)/(lambda-mustar*(1+0.62d0*lambda)))
  end function
  !
END PROGRAM
