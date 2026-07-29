!***********************************************************************
! project: MHDG
! file: venus_physics.f90
! description: Venus-specific helpers for the modular neutral coupling.
!              Atomic rates, Dnn and the source assembly live in
!              physics.f90 (compute_neutral_sources & co.) and are NOT
!              duplicated here: this module only carries plasma-state
!              admissibility checks, wall BC coefficients and the
!              finite-difference self-test of the shared kernels.
!***********************************************************************
MODULE venus_physics
  USE globals
  USE physics
  USE types, ONLY: bc_Bohm, bc_BohmPump, bc_BohmPuff
  IMPLICIT NONE

  REAL*8, PARAMETER :: venus_state_tol = 1.d-20
  REAL*8, PARAMETER :: venus_rho_floor = 1.d-12
  REAL*8, PARAMETER :: venus_energy_floor = 1.d-20

CONTAINS

  SUBROUTINE venus_neutral_wall_bc_coefficients(bc_type, recycling_coeff, cryopump_coeff, puff_coeff)
    INTEGER, INTENT(IN)  :: bc_type
    REAL*8, INTENT(OUT) :: recycling_coeff, cryopump_coeff, puff_coeff

    SELECT CASE (bc_type)
    CASE (bc_Bohm)
       recycling_coeff = phys%Re
       cryopump_coeff = 0.d0
       puff_coeff = 0.d0
    CASE (bc_BohmPump)
       recycling_coeff = phys%Re_pump
       puff_coeff = 0.d0
       IF (Mesh%pump_area > 1.d-10) THEN
          cryopump_coeff = phys%cryopump_power / (Mesh%pump_area * phys%lscale**2) &
               / simpar%refval_diffusion * phys%lscale
       ELSE
          cryopump_coeff = 0.d0
       END IF
    CASE (bc_BohmPuff)
       recycling_coeff = phys%Re
       cryopump_coeff = 0.d0
       IF (Mesh%puff_area > 1.d-10) THEN
          puff_coeff = phys%puff / simpar%refval_density / (Mesh%puff_area * phys%lscale**2) &
               / simpar%refval_diffusion * phys%lscale
       ELSE
          puff_coeff = 0.d0
       END IF
    CASE DEFAULT
       recycling_coeff = phys%Re
       cryopump_coeff = 0.d0
       puff_coeff = 0.d0
    END SELECT
  END SUBROUTINE venus_neutral_wall_bc_coefficients

  FUNCTION venus_plasma_state_is_admissible(U) RESULT(is_ok)
    REAL*8, INTENT(IN) :: U(:)
    LOGICAL :: is_ok
    REAL*8 :: ne, nEe, e_int_ion

    is_ok = .TRUE.
    IF (ANY(ISNAN(U(1:4)))) is_ok = .FALSE.
    IF (.NOT. is_ok) RETURN

    ne = U(1)
    nEe = U(4)
    IF (ne <= venus_rho_floor) is_ok = .FALSE.
    IF (nEe <= venus_energy_floor) is_ok = .FALSE.
    IF (.NOT. is_ok) RETURN

    e_int_ion = U(3) - 0.5d0*U(2)**2/ne
    IF (e_int_ion <= venus_energy_floor) is_ok = .FALSE.
  END FUNCTION venus_plasma_state_is_admissible

  SUBROUTINE venus_check_plasma_solution(u, is_ok, bad_index)
    REAL*8, INTENT(IN)  :: u(:)
    LOGICAL, INTENT(OUT) :: is_ok
    INTEGER, INTENT(OUT) :: bad_index
    INTEGER :: neq, Npel, iel, i, j
    REAL*8 :: Uloc(4)

    is_ok = .TRUE.
    bad_index = 0
    neq = phys%Neq
    Npel = Mesh%Nnodesperelem
    DO iel = 1, Mesh%Nelems
       DO i = 1, Npel
          DO j = 1, neq
             Uloc(j) = u(((iel - 1)*Npel + i - 1)*neq + j)
          END DO
          IF (.NOT. venus_plasma_state_is_admissible(Uloc)) THEN
             is_ok = .FALSE.
             bad_index = ((iel - 1)*Npel + i - 1)*neq + 1
             RETURN
          END IF
       END DO
    END DO
  END SUBROUTINE venus_check_plasma_solution

  SUBROUTINE venus_report_inadmissible_dof(u, dof_index)
    REAL*8, INTENT(IN) :: u(:)
    INTEGER, INTENT(IN) :: dof_index
    INTEGER :: neq, Npel, node_idx, iel, inode, j
    REAL*8 :: Uloc(4), e_int_ion

    neq = phys%Neq
    Npel = Mesh%Nnodesperelem
    node_idx = (dof_index - 1) / neq
    iel = node_idx / Npel + 1
    inode = MOD(node_idx, Npel) + 1
    DO j = 1, neq
       Uloc(j) = u(node_idx*neq + j)
    END DO

    WRITE (6, '("   [ADMISSIBILITY] element=", I6, " node=", I4, " dof_index=", I8)') &
         iel, inode, dof_index
    WRITE (6, '("   [ADMISSIBILITY] U = [n, Gamma, nEi, nEe] = ", 4(E14.6, 1X))') Uloc(1:4)

    IF (Uloc(1) <= venus_rho_floor) THEN
       WRITE (6, '("   [ADMISSIBILITY] fail: n <= ", E12.5)') venus_rho_floor
    END IF
    IF (Uloc(4) <= venus_energy_floor) THEN
       WRITE (6, '("   [ADMISSIBILITY] fail: nEe <= ", E12.5)') venus_energy_floor
    END IF
    IF (Uloc(1) > venus_rho_floor) THEN
       e_int_ion = Uloc(3) - 0.5d0*Uloc(2)**2/Uloc(1)
       WRITE (6, '("   [ADMISSIBILITY] ion internal energy nEi-0.5*Gamma^2/n = ", E14.6)') e_int_ion
       IF (e_int_ion <= venus_energy_floor) THEN
          WRITE (6, '("   [ADMISSIBILITY] fail: ion internal energy <= ", E12.5)') venus_energy_floor
       END IF
    END IF
  END SUBROUTINE venus_report_inadmissible_dof

#if defined(NEUTRAL) || defined(VENUS)
  SUBROUTINE venus_verify_source_jacobian(Utest, max_rel_err)
    ! Finite-difference check of compute_neutral_sources (physics.f90).
    ! The quasilinear convention gives back the physical source exactly:
    !   S = -(Sn0 + Sn.U)   since   Sn = -dS/dU, Sn0 = -S + (dS/dU).U
    ! so dS/dU_j by central differences must match -Sn(:,j).
    REAL*8, INTENT(IN)  :: Utest(5)
    REAL*8, INTENT(OUT) :: max_rel_err
    REAL*8 :: Up(5), Um(5)
    REAL*8 :: Sn(5, 5), Sn0(5), Snp(5, 5), Sn0p(5), Snm(5, 5), Sn0m(5)
    REAL*8 :: S_plus(5), S_minus(5), fd, eps, rel_err, sn_scale
    INTEGER :: i, j

    max_rel_err = 0.d0
    CALL compute_neutral_sources(Utest, Sn, Sn0)
    sn_scale = MAXVAL(ABS(Sn))

    DO j = 1, 5
       eps = 1.d-6 * MAX(ABS(Utest(j)), 1.d0)
       Up = Utest
       Up(j) = Up(j) + eps
       CALL compute_neutral_sources(Up, Snp, Sn0p)
       S_plus = -(Sn0p + MATMUL(Snp, Up))
       Um = Utest
       Um(j) = Um(j) - eps
       CALL compute_neutral_sources(Um, Snm, Sn0m)
       S_minus = -(Sn0m + MATMUL(Snm, Um))
       DO i = 1, 5
          fd = (S_plus(i) - S_minus(i)) / (2.d0*eps)
          rel_err = ABS(fd + Sn(i, j)) &
               / MAX(ABS(Sn(i, j)), 1.d-8*sn_scale, 1.d-30)
          max_rel_err = MAX(max_rel_err, rel_err)
       END DO
    END DO
  END SUBROUTINE venus_verify_source_jacobian
#endif

  SUBROUTINE venus_run_jacobian_selftest()
#if defined(NEUTRAL) || defined(VENUS)
    REAL*8 :: Utest(5), max_rel_err
    INTEGER :: idx_rhon_save, mpi_rank

    ! compute_neutral_sources reads U(phys%idx_rhon_eq); in the 4-equation
    ! build it is 0, so point it at the 5th slot of the extended test vector
    ! for the duration of the check only.
    ! Test n_n: any ordinary non-dimensional value is sufficient here, since
    ! the point is to differentiate the sources, not to represent a physical
    ! state.
    Utest(1:4) = (/1.d0, 0.1d0, 0.5d0, 0.4d0/)
    Utest(5) = 1.d-2
    idx_rhon_save = phys%idx_rhon_eq
    phys%idx_rhon_eq = 5
    CALL venus_verify_source_jacobian(Utest, max_rel_err)
    phys%idx_rhon_eq = idx_rhon_save
#ifdef PARALL
    mpi_rank = MPIvar%glob_id
#else
    mpi_rank = 0
#endif
    IF (mpi_rank == 0) THEN
       WRITE (6, '("   [VENUS FD TEST] max relative Jacobian error: ", E12.5)') max_rel_err
       IF (max_rel_err > 1.d-2) THEN
          PRINT *, 'ERROR: compute_neutral_sources finite-difference test failed.'
          STOP
       END IF
    END IF
#else
    PRINT *, 'ERROR: neutral coupling self-test requires a build with neutral physics'
    PRINT *, '       (compile with MODE=VENUS or a NEUTRAL model).'
    STOP
#endif
  END SUBROUTINE venus_run_jacobian_selftest

END MODULE venus_physics
