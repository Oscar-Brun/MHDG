!***********************************************************************
! project: MHDG
! file: neutral_coupling.f90
! description: Modular coupling interface between HDG plasma and neutrals.
!***********************************************************************
MODULE neutral_coupling
  USE globals
  USE MPI_OMP
  USE types, ONLY: bc_Bohm, bc_BohmPump, bc_BohmPuff
  USE venus_types
  USE venus_mesh
  USE venus_physics
#if defined(NEUTRAL) || defined(VENUS)
  ! Atomic rates, Dnn and sources: THE routines of the monolithic model
  ! (physics.f90), no local reimplementation.
  ! Diffusive: compute_sigmavcx/vnn are not called directly - CX enters
  ! inside compute_Dnn, through the transport denominator.
  ! Euler: they ARE called directly, because the friction term of the
  ! momentum equation must rebuild EXACTLY that denominator
  ! (nu_fric + nu_iz = denominator of compute_Dnn), which is what makes the
  ! collisional limit exact by construction.
  USE physics, ONLY: compute_Ti, compute_Te, compute_sigmaviz, &
       &compute_sigmavrec, compute_sigmavcx, compute_sigmavnn, &
       &compute_Dnn, compute_neutral_sources, setLocalDiff, &
       &compute_neutral_free_streaming_speed
#else
  USE physics, ONLY: compute_Ti, compute_Te
#endif
  IMPLICIT NONE
  SAVE

  REAL*8, PARAMETER :: nn_floor = 1.d-20
  ! Gauss-Seidel tolerance of the neutral solve. The number of sweeps scales
  ! like ln(tol) (geometric GS convergence), and 1e-4 is enough here: the
  ! neutral field is then accurate to 1e-4 relative, far below the threshold
  ! at which the fields are compared (~1%), and the coupling still converges
  ! under picard_tol = 1e-3. Tighten to 1e-8 for a machine-accuracy steady
  ! state.
  REAL*8, PARAMETER :: fv_gs_tol = 1.d-4
  ! Cap on the number of Gauss-Seidel sweeps. The FV system is quasi-Neumann
  ! (boundary faces carry a flux only, no diffusive coefficient), so the level
  ! of n_n is set solely by ionization and by the BDF term, and GS contracts
  ! very slowly there: several thousand sweeps are needed at a tight tolerance.
  ! 50000 leaves a factor 10 of margin for ~1 s of compute, which is nothing
  ! next to the PASTIX factorization that every Newton-Raphson iteration pays
  ! for more than ten times over.
  ! A true Krylov solver (BiCGSTAB) remains desirable for larger meshes, but is
  ! NOT a blocking issue at this size.
  INTEGER, PARAMETER :: fv_gs_max_iter = 50000

  ! ----- Euler model: local pseudo-time marching -----
  ! CFL of the explicit first-order HLL step: the positivity proof requires
  ! dt <= vol / sum(smax_f * L_f * R_f), i.e. CFL <= 1; 0.8 leaves some margin
  ! without needlessly slowing convergence towards the steady state.
  REAL*8, PARAMETER :: euler_cfl = 0.8d0
  ! Same tolerance as the diffusive Gauss-Seidel: the criterion is the relative
  ! change of the state between two pseudo-steps (homogeneous norm, with Gamma
  ! weighted by 1/c_n so that it is comparable to n_n).
  REAL*8, PARAMETER :: euler_tol = 1.d-8
  INTEGER, PARAMETER :: euler_max_iter = 100000
  ! Near-vacuum safeguard: |u_n| <= euler_umax_cn * c_n. It only bites in core
  ! cells where n_n sits at the floor (1e-20) and where Gamma/n_n would lose
  ! all meaning; everywhere else the neutral velocities stay of order c_n.
  REAL*8, PARAMETER :: euler_umax_cn = 20.d0

  ! Pseudo-transient continuation: ceiling on the time step for the SER ramp
  ! (compute_dt), enabled under active coupling in the time loop (MHDG.f90).
  ! The steady state of the segregated coupling is reached transiently with a
  ! GROWING but BOUNDED dt, not with steady = .true. (pure Newton, which breaks
  ! the Picard loop by removing the dt anchor). Bounded, because too large a dt
  ! makes the BDF diagonal of the FV neutral solver vanish -> Gauss-Seidel then
  ! stops converging. 1e3 = 10x the time-init dt0: conservative, safe; raise it
  ! for speed once convergence is established.
  REAL*8, PARAMETER :: venus_dt_max = 1.d3

  LOGICAL :: coupling_active = .FALSE.
  ! True as long as the FV solver has never solved for n_n: the very first
  ! solve must not be under-relaxed against the seed, which means nothing.
  LOGICAL :: nn_never_solved = .TRUE.

  REAL*8, ALLOCATABLE, TARGET :: neutral_density_gauss(:, :)
  REAL*8, ALLOCATABLE, TARGET :: neutral_vx_gauss(:, :)
  REAL*8, ALLOCATABLE, TARGET :: neutral_vy_gauss(:, :)
  REAL*8, ALLOCATABLE, TARGET :: neutral_source_gauss(:, :, :)
  REAL*8, ALLOCATABLE, TARGET :: neutral_jac_source_gauss(:, :, :, :)

  REAL*8, ALLOCATABLE :: u_plasma_old(:)
  ! Time history of n_n (P0 per cell) for the BDF term.
  ! Same convention as sol%u0: slot 1 = previous time step, slot 2 = the one
  ! before that, etc. Advanced ONCE per time step, in
  ! neutral_coupling_begin_timestep - never in solve_step, which may be called
  ! again at every Picard iteration of the same step.
  REAL*8, ALLOCATABLE :: nn_time_hist(:, :)
  ! Wall sources per cell, integrated at the 1D Gauss points of the exterior
  ! faces by neutral_coupling_wall_moments.
  !   wall_influx(K) = INT_dK [Re*(outgoing plasma flux) + puff] dline(R)
  !   wall_pump(K)   = INT_dK pump dline(R)   (multiplies n_n on the diagonal)
  REAL*8, ALLOCATABLE :: wall_influx(:), wall_pump(:)
  TYPE(venus_mesh_t) :: v_mesh

CONTAINS

  LOGICAL FUNCTION neutral_coupling_is_active()
    neutral_coupling_is_active = coupling_active
  END FUNCTION neutral_coupling_is_active

  ! -------------------------------------------------------------------
  ! Venus banner (segregated coupling): purely cosmetic, rank 0 only.
  ! Called EARLY (right after read_input, MHDG.f90) so that it shows up BEFORE
  ! the heavy initialization (mesh, symbolic factorization) - the first thing
  ! one sees when starting a Venus run. Venus-only: it is printed only when a
  ! neutral model is coupled (Diffusion/Venus), so the legacy monolithic path
  ! prints nothing. "Venus glow" gradient (gold -> pink), ANSI 256 colors.
  ! -------------------------------------------------------------------
  SUBROUTINE neutral_coupling_banner()
    CHARACTER(LEN=1), PARAMETER :: esc = ACHAR(27)   ! ESC for the ANSI color codes

    IF (TRIM(ADJUSTL(simpar%neutral_model_type)) /= 'Diffusion' .AND. &
         & TRIM(ADJUSTL(simpar%neutral_model_type)) /= 'Venus') RETURN
    IF (MPIvar%glob_id /= 0) RETURN

    WRITE(6,'(A)') ''
    ! --- SOLEDGE-HDG (cyan gradient) ---
    WRITE(6,'(A)') esc//'[38;5;39m'//' ____   ___  _     _____ ____   ____ _____      _   _ ____   ____ '//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;45m'//'/ ___| / _ \| |   | ____|  _ \ / ___| ____|    | | | |  _ \ / ___|'//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;51m'//'\___ \| | | | |   |  _| | | | | |  _|  _| _____| |_| | | | | |  _ '//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;45m'//' ___) | |_| | |___| |___| |_| | |_| | |___|_____|  _  | |_| | |_| |'//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;39m'//'|____/ \___/|_____|_____|____/ \____|_____|    |_| |_|____/ \____|'//esc//'[0m'
    WRITE(6,'(A)') ''
    ! --- VENUS (gold gradient) + Venus planet (amber gradient) side by side ---
    WRITE(6,'(A)') esc//'[38;5;227m'//'__     _______ _   _ _   _ ____  '//esc//'[0m'//'   '//esc//'[38;5;230m'//'  .-~~~~-.  '//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;221m'//'\ \   / / ____| \ | | | | / ___| '//esc//'[0m'//'   '//esc//'[38;5;222m'//' / ~ ~ ~  \ '//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;215m'//' \ \ / /|  _| |  \| | | | \___ \ '//esc//'[0m'//'   '//esc//'[38;5;214m'//'| ~ ~ ~ ~  |'//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;209m'//'  \ V / | |___| |\  | |_| |___) | '//esc//'[0m'//'   '//esc//'[38;5;208m'//' \ ~ ~ ~  / '//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;213m'//'   \_/  |_____|_| \_|\___/|____/ '//esc//'[0m'//'   '//esc//'[38;5;172m'//'  ''-~~~~-''  '//esc//'[0m'
    WRITE(6,'(A)') esc//'[1;38;5;213m'//'   segregated plasma / neutral coupling'//esc//'[0m'
    WRITE(6,'(A)') esc//'[38;5;245m'//'   neutral model : '//TRIM(ADJUSTL(simpar%neutral_model_type))//esc//'[0m'
    WRITE(6,'(A)') ''
  END SUBROUTINE neutral_coupling_banner

  SUBROUTINE neutral_coupling_init()
    INTEGER :: neq_plasma

    coupling_active = (TRIM(ADJUSTL(simpar%neutral_model_type)) == 'Diffusion' &
         & .OR. TRIM(ADJUSTL(simpar%neutral_model_type)) == 'Venus')
    neq_plasma = phys%Neq

    PRINT *, '--- Neutral Coupling: Initializing API ---'
    PRINT *, '   -> Selected model: ', TRIM(ADJUSTL(simpar%neutral_model_type))
    IF (.NOT. coupling_active) THEN
       PRINT *, '--- Neutral Coupling: inactive (plasma-only path) ---'
       RETURN
    END IF

    ! The segregated coupling assumes a 4-equation plasma: in a build where
    ! the neutral is already an HDG equation (idx_rhon_eq > 0), the two models
    ! would fight over the same field.
    IF (phys%idx_rhon_eq /= 0) THEN
       PRINT *, 'ERROR: the segregated neutral coupling requires the 4-equation plasma build'
       PRINT *, '       (this build already carries a native HDG neutral equation).'
       STOP
    END IF

    IF (ALLOCATED(neutral_density_gauss)) DEALLOCATE(neutral_density_gauss)
    IF (ALLOCATED(neutral_vx_gauss)) DEALLOCATE(neutral_vx_gauss)
    IF (ALLOCATED(neutral_vy_gauss)) DEALLOCATE(neutral_vy_gauss)
    IF (ALLOCATED(neutral_source_gauss)) DEALLOCATE(neutral_source_gauss)
    IF (ALLOCATED(neutral_jac_source_gauss)) DEALLOCATE(neutral_jac_source_gauss)
    IF (ALLOCATED(u_plasma_old)) DEALLOCATE(u_plasma_old)
    IF (ALLOCATED(nn_time_hist)) DEALLOCATE(nn_time_hist)
    IF (ALLOCATED(wall_influx)) DEALLOCATE(wall_influx)
    IF (ALLOCATED(wall_pump)) DEALLOCATE(wall_pump)

    ALLOCATE(neutral_density_gauss(Mesh%Nelems, refElPol%NGauss2d))
    ALLOCATE(neutral_vx_gauss(Mesh%Nelems, refElPol%NGauss2d))
    ALLOCATE(neutral_vy_gauss(Mesh%Nelems, refElPol%NGauss2d))
    ALLOCATE(neutral_source_gauss(neq_plasma, Mesh%Nelems, refElPol%NGauss2d))
    ALLOCATE(neutral_jac_source_gauss(neq_plasma, neq_plasma, Mesh%Nelems, refElPol%NGauss2d))
    ! Slots zeroed as for sol%u0 (initialize_solu0_uiter_uconv): the slots
    ! beyond the history actually elapsed are never read, the BDF ramp-up
    ! (order min(time%ik, time%tis)) takes care of that. Slot 1 is filled by
    ! the first neutral_coupling_begin_timestep, with the seed value.
    ALLOCATE(nn_time_hist(Mesh%Nelems, time%tis))
    nn_time_hist = 0.d0
    ALLOCATE(wall_influx(Mesh%Nelems))
    ALLOCATE(wall_pump(Mesh%Nelems))
    wall_influx = 0.d0
    wall_pump = 0.d0

    ! The wall source replicates assembly_bohm_bc (hdg_BC.f90:2256-2268)
    ! WITHOUT three of its optional branches. Refusing them explicitly is
    ! better than ignoring them silently:
    !   - transport_1d: carries the pinch (APinch) and modifies the diffusion;
    !   - import_diffusion_1D: modifies the diffusion at the boundary;
    !   - apply_trim: reflection coefficient RN(E, theta) instead of 1.
    IF (switch%transport_1d .OR. switch%import_diffusion_1D .OR. phys%apply_trim) THEN
       PRINT *, 'ERROR: the segregated neutral coupling does not support transport_1d,'
       PRINT *, '       import_diffusion_1D or apply_trim (wall recycling would'
       PRINT *, '       silently diverge from the monolithic formula).'
       STOP
    END IF

    neutral_density_gauss = 0.d0
    neutral_vx_gauss = 0.d0
    neutral_vy_gauss = 0.d0
    neutral_source_gauss = 0.d0
    neutral_jac_source_gauss = 0.d0

    CALL venus_mesh_build(v_mesh)
    CALL venus_run_jacobian_selftest()
    ! Manufactured-solution check of the FV scheme. A no-op unless VENUS_MMS=1
    ! is set in the environment, in which case it measures, prints and STOPs.
    CALL venus_run_mms_selftest()
    nn_never_solved = .TRUE.
    CALL neutral_coupling_seed_initial_density()

    ALLOCATE(u_plasma_old(SIZE(sol%u)))
    u_plasma_old = sol%u

    PRINT *, '--- Neutral Coupling: API successfully initialized ---'
  END SUBROUTINE neutral_coupling_init

  SUBROUTINE neutral_coupling_begin_timestep(u_current)
    REAL*8, INTENT(IN) :: u_current(:)
    INTEGER :: iord, iel

    IF (.NOT. coupling_active) RETURN
    IF (.NOT. ALLOCATED(u_plasma_old)) THEN
       ALLOCATE(u_plasma_old(SIZE(u_current)))
    ELSEIF (SIZE(u_plasma_old) /= SIZE(u_current)) THEN
       DEALLOCATE(u_plasma_old)
       ALLOCATE(u_plasma_old(SIZE(u_current)))
    END IF
    u_plasma_old = u_current

    ! Shift of the neutral BDF history, at the same point of the cycle as the
    ! plasma (update_solution shifts sol%u0 at the end of an accepted step;
    ! here we are at the beginning of the next step, where v_mesh%cells(:)%U(1)
    ! holds exactly the converged value of the previous step - the seed at the
    ! very first step). The two instants are equivalent, because under coupling
    ! there is no rollback path that keeps the loop running: a non-converged NR
    ! or Picard iteration ends in a STOP.
    ! Full-window shift where update_solution only shifts the filled slots:
    ! identical in content, since the empty slots are 0 on both sides.
    IF (ALLOCATED(nn_time_hist)) THEN
       DO iord = SIZE(nn_time_hist, 2), 2, -1
          nn_time_hist(:, iord) = nn_time_hist(:, iord - 1)
       END DO
       DO iel = 1, v_mesh%ncells
          nn_time_hist(iel, 1) = v_mesh%cells(iel)%U(1)
       END DO
    END IF
  END SUBROUTINE neutral_coupling_begin_timestep


  SUBROUTINE neutral_coupling_seed_initial_density()
    INTEGER :: iel
    REAL*8 :: nn_seed

    IF (.NOT. coupling_active) RETURN
    ! Seed value for n_n. The monolithic model starts from zero
    ! (analytical.f90:100: up(ind,inn) = 0.) and lets recycling fill the
    ! domain: we do the same, except that the FV solver needs a value strictly
    ! above nn_floor (Dnn reads n_n, and the relaxation treats a very small
    ! nn_old as "not initialized yet").
    nn_seed = 100.d0*nn_floor
    DO iel = 1, v_mesh%ncells
       v_mesh%cells(iel)%U(1) = nn_seed
       v_mesh%cells(iel)%U(2:3) = 0.d0
    END DO
  END SUBROUTINE neutral_coupling_seed_initial_density

  SUBROUTINE neutral_coupling_clear_fields()
    IF (.NOT. ALLOCATED(neutral_density_gauss)) RETURN
    neutral_density_gauss = 0.d0
    neutral_vx_gauss = 0.d0
    neutral_vy_gauss = 0.d0
    neutral_source_gauss = 0.d0
    neutral_jac_source_gauss = 0.d0
  END SUBROUTINE neutral_coupling_clear_fields

  SUBROUTINE neutral_coupling_solve_step(ures)
    REAL*8, INTENT(IN) :: ures(:)

    SELECT CASE (TRIM(ADJUSTL(simpar%neutral_model_type)))
    CASE ('None')
       RETURN
    CASE ('Venus')
       CALL neutral_euler_solve_step(ures)
    CASE ('Diffusion')
       CALL neutral_diffusion_solve_step(ures)
    CASE DEFAULT
       PRINT *, 'ERROR: Unknown neutral_model_type in neutral_coupling_solve_step'
       STOP
    END SELECT
  END SUBROUTINE neutral_coupling_solve_step

  SUBROUTINE neutral_diffusion_solve_step(ures)
    REAL*8, INTENT(IN) :: ures(:)

    INTEGER :: iel, i, g, j, Npel, idx_start
    INTEGER :: iter, ifa, iface_loc, cell_L, cell_R, neighbor
    INTEGER :: neq_plasma
    REAL*8 :: nn_old, n_n_relaxed, vt_n, rel_res, field_norm, new_norm
    REAL*8 :: dx, dy, d_centers, D_face, sum_offdiag
    REAL*8 :: U_cell(4), U_g(4), n_n
    REAL*8 :: U5_cell(5), U5_g(5), Sn5(5, 5), Sn05(5)
    INTEGER :: idx_rhon_save
    REAL*8, ALLOCATABLE :: ue_el(:, :)
    REAL*8, ALLOCATABLE :: diag_vec(:), rhs_vec(:), Dnn(:), n_n_vec(:), n_n_prev(:)
    ! Free-streaming speed per cell (for the limiter bound) and transport
    ! coefficient per face - precomputed ONCE and consumed identically by the
    ! assembly and by the GS sweep, which removes the duplicated formula
    ! between the two (fragile: any discrepancy means two different systems).
    REAL*8, ALLOCATABLE :: cn_cell(:), coeff_face(:)
    REAL*8 :: nn_L, nn_R, n_n_face, cn_face, grad_n
    REAL*8 :: gamma_unlim, gamma_max, gamma_abs, ratio_lim, phi_lim
    LOGICAL :: limiter_on
    ! Moments integrated with the element quadrature, NOT cell-averaged rates.
    !   mom_iz(K)  = INT_K ne*sigmaviz    dOmega   [neutral sink, per unit n_n]
    !   mom_rec(K) = INT_K ne^2*sigmavrec dOmega   [neutral source]
    REAL*8, ALLOCATABLE :: mom_iz(:), mom_rec(:)
    ! Cell volume INT_K R dOmega, accumulated with the SAME quadrature as the
    ! moments - the time term must live in the same geometry as the reactions
    ! and the transport.
    REAL*8, ALLOCATABLE :: cell_vol(:)
    REAL*8 :: ktis(time%tis + 1)
    INTEGER :: iord
    REAL*8 :: Xel(Mesh%Nnodesperelem, 2)
    REAL*8 :: xyg(refElPol%Ngauss2d, 2)
    REAL*8 :: J11(refElPol%Ngauss2d), J12(refElPol%Ngauss2d)
    REAL*8 :: J21(refElPol%Ngauss2d), J22(refElPol%Ngauss2d)
    REAL*8 :: detJ(refElPol%Ngauss2d), dvolu, r_face
    REAL*8 :: sviz_g, svrec_g
    REAL*8, ALLOCATABLE :: ueg(:, :)
    REAL*8 :: bal_plasma, bal_neutral
    ! Second-order MUSCL (switch%neutral_muscl): limited cell gradient (lagged,
    ! Green-Gauss + Barth-Jespersen) and linear reconstruction of n_n. grad = 0
    ! (switch off) => the first-order scheme, bit for bit. Conservation is kept:
    ! n_n is reconstructed on BOTH sides (neutral sink off_iz AND plasma source
    ! nn_g).
    REAL*8, ALLOCATABLE :: grad_nn(:, :)   ! (2, ncells) limited gradient of n_n
    REAL*8, ALLOCATABLE :: off_iz(:)       ! ionization offset of the linear term
    REAL*8 :: gx, gy, nnf, cx, cy, nn_lo, nn_hi, dvert, phi_bj, nn_g, xnx, xny
    REAL*8 :: Xv(3, 2)
    INTEGER :: ivert

#if defined(NEUTRAL) || defined(VENUS)
    IF (MPIvar%glob_id == 0) THEN
       PRINT *, '   [NEUTRAL SOLVER] Solving 2D Finite Volumes Neutral Diffusion...'
    END IF

    neq_plasma = phys%Neq

    ! Wall sources. BEFORE the idx_rhon_eq bracket below - setLocalDiff writes
    ! d_iso(inn,inn,:) under the guard inn > 0, and inn = 5 would overflow the
    ! arrays dimensioned with Neq = 4.
    CALL neutral_coupling_wall_moments(ures)

    ! The compute_* routines of physics.f90 read U(phys%idx_rhon_eq). In a
    ! 4-equation build it is 0 (set_model_layout): we point it at the 5th
    ! component of the extended vector [U1..U4, n_n] FOR THE DURATION OF THIS
    ! ROUTINE only - setting it globally would enable unguarded runtime readers
    ! that index arrays dimensioned with Neq.
    idx_rhon_save = phys%idx_rhon_eq
    phys%idx_rhon_eq = neq_plasma + 1
    Npel = Mesh%Nnodesperelem
    ALLOCATE(ue_el(Npel, neq_plasma))
    ALLOCATE(diag_vec(v_mesh%ncells))
    ALLOCATE(rhs_vec(v_mesh%ncells))
    ALLOCATE(Dnn(v_mesh%ncells))
    ALLOCATE(n_n_vec(v_mesh%ncells))
    ALLOCATE(n_n_prev(v_mesh%ncells))
    ALLOCATE(mom_iz(v_mesh%ncells))
    ALLOCATE(mom_rec(v_mesh%ncells))
    ALLOCATE(cell_vol(v_mesh%ncells))
    ALLOCATE(cn_cell(v_mesh%ncells))
    ALLOCATE(coeff_face(v_mesh%nfaces))
    limiter_on = (TRIM(ADJUSTL(phys%neutral_flux_limiter_mode)) == 'lagged_flux_limiter')
    ALLOCATE(ueg(refElPol%Ngauss2d, neq_plasma))
    ALLOCATE(grad_nn(2, v_mesh%ncells))
    ALLOCATE(off_iz(v_mesh%ncells))

    ! ---- Second-order MUSCL: limited cell gradient, built on the INCOMING n_n
    !      (lagged) -> stays linear inside the implicit solve. Green-Gauss, then
    !      a Barth-Jespersen limiter on the 3 vertices (bounds the linear part
    !      => positivity, no new extremum). grad_nn = 0 if the switch is off. --
    grad_nn = 0.d0
    off_iz  = 0.d0
    IF (switch%neutral_muscl) THEN
       DO iel = 1, v_mesh%ncells
          n_n = MAX(v_mesh%cells(iel)%U(1), nn_floor)
          cx = v_mesh%cells(iel)%center(1)
          cy = v_mesh%cells(iel)%center(2)
          gx = 0.d0; gy = 0.d0
          nn_lo = n_n; nn_hi = n_n
          DO iface_loc = 1, 3
             ifa = v_mesh%cell_faces(iface_loc, iel)
             cell_L = v_mesh%faces(ifa)%cell_L
             cell_R = v_mesh%faces(ifa)%cell_R
             IF (cell_L == iel) THEN
                neighbor = cell_R
                xnx =  v_mesh%faces(ifa)%normal(1)   ! outgoing normal (L->R)
                xny =  v_mesh%faces(ifa)%normal(2)
             ELSE
                neighbor = cell_L
                xnx = -v_mesh%faces(ifa)%normal(1)   ! iel = R: flip the sign
                xny = -v_mesh%faces(ifa)%normal(2)
             END IF
             IF (neighbor > 0) THEN
                nnf = MAX(v_mesh%cells(neighbor)%U(1), nn_floor)
                nn_lo = MIN(nn_lo, nnf)
                nn_hi = MAX(nn_hi, nnf)
                nnf = 0.5d0*(n_n + nnf)
             ELSE
                nnf = n_n                            ! wall: zero gradient
             END IF
             gx = gx + nnf * xnx * v_mesh%faces(ifa)%length
             gy = gy + nnf * xny * v_mesh%faces(ifa)%length
          END DO
          gx = gx / MAX(v_mesh%cells(iel)%area, 1.d-30)
          gy = gy / MAX(v_mesh%cells(iel)%area, 1.d-30)
          ! Barth-Jespersen limiter on the 3 vertices of the triangle
          Xv(1, :) = Mesh%X(Mesh%Tlin(iel, 1), :)
          Xv(2, :) = Mesh%X(Mesh%Tlin(iel, 2), :)
          Xv(3, :) = Mesh%X(Mesh%Tlin(iel, 3), :)
          phi_bj = 1.d0
          DO ivert = 1, 3
             dvert = gx*(Xv(ivert,1) - cx) + gy*(Xv(ivert,2) - cy)
             IF (dvert > 1.d-30) THEN
                phi_bj = MIN(phi_bj, MIN(1.d0, (nn_hi - n_n)/dvert))
             ELSE IF (dvert < -1.d-30) THEN
                phi_bj = MIN(phi_bj, MIN(1.d0, (nn_lo - n_n)/dvert))
             END IF
          END DO
          grad_nn(1, iel) = phi_bj * gx
          grad_nn(2, iel) = phi_bj * gy
       END DO
    END IF

    DO iel = 1, Mesh%Nelems
       idx_start = (iel - 1) * Npel + 1
       DO i = 1, Npel
          ue_el(i, 1) = ures((idx_start + i - 2) * neq_plasma + 1)
          ue_el(i, 2) = ures((idx_start + i - 2) * neq_plasma + 2)
          ue_el(i, 3) = ures((idx_start + i - 2) * neq_plasma + 3)
          ue_el(i, 4) = ures((idx_start + i - 2) * neq_plasma + 4)
       END DO

       ! The cell average is now only used for Dnn (transport, not sources).
       ! Sources go through the quadrature; the wall goes through
       ! neutral_coupling_wall_moments.
       DO j = 1, neq_plasma
          U_cell(j) = SUM(ue_el(:, j)) / Npel
       END DO

       nn_old = MAX(v_mesh%cells(iel)%U(1), nn_floor)
       U5_cell(1:4) = U_cell(1:4)
       U5_cell(5) = nn_old
       CALL compute_Dnn(U5_cell, Dnn(iel))
       ! cn = SQRT(Mref*Ti_limited) (physics.f90:895), reads the plasma only -
       ! the free-streaming bound of the limiter needs it at the faces.
       CALL compute_neutral_free_streaming_speed(U5_cell, cn_cell(iel))

       ! ----- Moments on the SAME quadrature as the HDG solver -----
       ! n_n is P0 over the cell, hence constant: it factors out of the
       ! integral and the source stays linear in n_n_cell. The neutral sink is
       ! then exactly the opposite of what the plasma gains, Gauss point by
       ! Gauss point. Pattern taken from
       ! account_neutral_reaction_source_totals (hdg_ComputeJacobian.f90:1394):
       ! same N2D, same weights, same detJ, same axisymmetric factor - that is
       ! the condition for conservation.
       Xel = Mesh%X(Mesh%T(iel, :), :)
       xyg = MATMUL(refElPol%N2D, Xel)
       ueg = MATMUL(refElPol%N2D, ue_el)
       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12

       mom_iz(iel) = 0.d0
       mom_rec(iel) = 0.d0
       cell_vol(iel) = 0.d0
       cx = v_mesh%cells(iel)%center(1)
       cy = v_mesh%cells(iel)%center(2)
       DO g = 1, refElPol%Ngauss2d
          dvolu = refElPol%gauss_weights2D(g) * detJ(g)
          IF (switch%axisym) dvolu = dvolu * xyg(g, 1)
          U5_g(1:4) = ueg(g, 1:4)
          U5_g(5) = nn_old
          CALL compute_sigmaviz(U5_g, sviz_g)
          CALL compute_sigmavrec(U5_g, svrec_g)
          mom_iz(iel)  = mom_iz(iel)  + ueg(g, 1) * sviz_g * dvolu
          mom_rec(iel) = mom_rec(iel) + ueg(g, 1)**2 * svrec_g * dvolu
          cell_vol(iel) = cell_vol(iel) + dvolu
          ! MUSCL: offset of the reconstructed linear term (lagged gradient,
          ! = 0 if off). Neutral ionization sink = mom_iz*n_n_cell + off_iz;
          ! the constant term off_iz goes to the right-hand side. The plasma
          ! sees the SAME nn_g (source loop below) -> conservation preserved.
          off_iz(iel) = off_iz(iel) + ueg(g, 1) * sviz_g * &
               & (grad_nn(1, iel)*(xyg(g,1) - cx) + grad_nn(2, iel)*(xyg(g,2) - cy)) * dvolu
       END DO
    END DO

    diag_vec = 0.d0
    rhs_vec = 0.d0
    DO iel = 1, v_mesh%ncells
       diag_vec(iel) = diag_vec(iel) + mom_iz(iel)
       rhs_vec(iel) = rhs_vec(iel) + mom_rec(iel) - off_iz(iel)
    END DO

    ! ----- BDF time term, same discretization as the plasma -----
    ! The monolithic solver assembles  Auu += ktis(1)*NNi/dt  and
    ! rhs += sum ktis(iord+1)*u0(:,iord)/dt  (hdg_ComputeJacobian.f90:2923
    ! and :3277), and skips all of it if switch%steady. Since n_n is P0, the FV
    ! projection of the same discretization is
    !   vol_K*ktis(1)/dt  on the diagonal,
    !   vol_K/dt * sum ktis(iord+1)*n_n^(step-iord)  on the right-hand side,
    ! with vol_K = INT_K R dOmega: same Gauss weights, same detJ, same
    ! axisymmetric factor as the reaction moments. Same start-up ramp too
    ! (ktis depends on time%ik, read at the same step on both sides): the FV
    ! solver and the plasma raise their BDF order together.
    IF (.NOT. switch%steady) THEN
       CALL venus_time_integration_coefficients(ktis)
       DO iel = 1, v_mesh%ncells
          diag_vec(iel) = diag_vec(iel) + ktis(1) * cell_vol(iel) / time%dt
          DO iord = 1, time%tis
             rhs_vec(iel) = rhs_vec(iel) &
                  & + ktis(iord + 1) * nn_time_hist(iel, iord) * cell_vol(iel) / time%dt
          END DO
       END DO
    END IF

    ! The moments above carry the axisymmetric factor R (dOmega = R.dA, the
    ! 2.pi being a global factor that cancels out). Transport and BCs must
    ! carry it too, otherwise the FV equation mixes two geometries.
    !
    ! The coefficient of every interior face is precomputed HERE and consumed
    ! as is by the assembly AND by the GS sweep: a single formula.
    !
    ! Neutral flux limiter, an FV transcription of
    ! compute_neutral_flux_limiter (physics.f90:904) - the shared kernel is NOT
    ! callable from here: it does RESHAPE(Q,(Ndim,simpar%Neq)) and then reads
    ! Qpr(:,inn) with inn = 5, out of bounds in a 4-equation build. A line by
    ! line transcription is used instead. "Lagged": phi is evaluated on the
    ! INCOMING n_n (previous time step / previous Picard iterate), hence frozen
    ! during the sweeps - the matrix stays linear in n_n, as in the monolithic
    ! model where phi multiplies Dnn at the Gauss point with the state of
    ! iteration k-1. Accepted difference: the two-point FV gradient only has
    ! the component NORMAL to the face, whereas the monolithic model takes the
    ! norm of the full vector. Wherever a tangential gradient exists,
    ! phi_FV >= phi_mono (less limiting). The two-point FV scheme does not see
    ! those gradients anyway.
    coeff_face = 0.d0
    DO ifa = 1, v_mesh%nfaces
       cell_L = v_mesh%faces(ifa)%cell_L
       cell_R = v_mesh%faces(ifa)%cell_R
       IF (cell_R <= 0) CYCLE   ! wall: handled per cell further down
       r_face = 1.d0
       IF (switch%axisym) r_face = v_mesh%faces(ifa)%midpoint(1)
       dx = v_mesh%cells(cell_L)%center(1) - v_mesh%cells(cell_R)%center(1)
       dy = v_mesh%cells(cell_L)%center(2) - v_mesh%cells(cell_R)%center(2)
       d_centers = MAX(SQRT(dx**2 + dy**2), 1.d-30)
       D_face = 0.5d0 * (Dnn(cell_L) + Dnn(cell_R))

       IF (limiter_on) THEN
          nn_L = MAX(v_mesh%cells(cell_L)%U(1), nn_floor)
          nn_R = MAX(v_mesh%cells(cell_R)%U(1), nn_floor)
          n_n_face = 0.5d0*(nn_L + nn_R)
          cn_face = 0.5d0*(cn_cell(cell_L) + cn_cell(cell_R))
          grad_n = (nn_R - nn_L)/d_centers
          ! physics.f90:923: gamma_unlim = -Dnn*Qpr(:,inn)
          gamma_unlim = D_face*ABS(grad_n)
          ! :934: gamma_max = MAX(fs_fraction*MAX(n_n,0)*cn, fs_flux_min)
          gamma_max = MAX(phys%neutral_flux_limiter_fs_fraction*n_n_face*cn_face, &
               &phys%neutral_flux_limiter_fs_flux_min)
          ! :936: norm smoothed by eps
          gamma_abs = SQRT(gamma_unlim**2 + phys%neutral_flux_limiter_eps**2)
          ! :938-947: ratio and phi = 1/(1+ratio)
          IF (gamma_max > 0.d0) THEN
             ratio_lim = gamma_abs/gamma_max
             phi_lim = 1.d0/(1.d0 + ratio_lim)
          ELSE
             phi_lim = 0.d0
          END IF
          ! :951: diffusion floor diff_nn_min
          IF (D_face > 0.d0) phi_lim = MAX(phi_lim, MIN(1.d0, phys%diff_nn_min/D_face))
          D_face = phi_lim*D_face
       END IF

       coeff_face(ifa) = D_face * v_mesh%faces(ifa)%length * r_face / d_centers
       diag_vec(cell_L) = diag_vec(cell_L) + coeff_face(ifa)
       diag_vec(cell_R) = diag_vec(cell_R) + coeff_face(ifa)
    END DO

    DO iel = 1, v_mesh%ncells
       diag_vec(iel) = diag_vec(iel) + wall_pump(iel)
       rhs_vec(iel) = rhs_vec(iel) + wall_influx(iel)
    END DO

    DO iel = 1, v_mesh%ncells
       n_n_vec(iel) = MAX(v_mesh%cells(iel)%U(1), nn_floor)
       IF (diag_vec(iel) <= 1.d-30) THEN
          PRINT *, 'ERROR: neutral diffusion diagonal is zero at cell ', iel
          STOP
       END IF
    END DO

    ! Stopping criterion: norm of the CHANGE divided by the norm of the field.
    ! Comparing the change of the NORM instead would be wrong in both
    ! directions: a sweep that redistributes n_n at constant norm would pass as
    ! converged, and the effective tolerance would be fv_gs_tol^2.
    !
    ! The normalization uses the field AFTER the sweep, not before. Normalizing
    ! by the incoming field would let the very first Picard iteration of the
    ! very first step declare convergence after ONE sweep: n_n is then the seed
    ! (~1e-18), so the norm falls below any absolute guard and the residual is
    ! reported as zero on a field that has solved nothing. Normalizing by the
    ! current field gives rel_res ~ 1 on the first sweep from zero - the truth -
    ! and is only ~0 when the solution itself is zero, where "converged" is the
    ! right answer.
    rel_res = 1.d0
    DO iter = 1, fv_gs_max_iter
       n_n_prev = n_n_vec
       DO iel = 1, v_mesh%ncells
          sum_offdiag = 0.d0
          DO iface_loc = 1, 3
             ifa = v_mesh%cell_faces(iface_loc, iel)
             cell_L = v_mesh%faces(ifa)%cell_L
             cell_R = v_mesh%faces(ifa)%cell_R
             IF (cell_L == iel) THEN
                neighbor = cell_R
             ELSE
                neighbor = cell_L
             END IF
             IF (neighbor > 0) THEN
                ! THE SAME coefficient as the assembly - precomputed in
                ! coeff_face, flux limiter included: the sweep and the matrix
                ! solve the same system by construction, not by duplicating a
                ! formula.
                sum_offdiag = sum_offdiag + coeff_face(ifa) * n_n_vec(neighbor)
             END IF
          END DO
          n_n_vec(iel) = MAX((rhs_vec(iel) + sum_offdiag) / diag_vec(iel), nn_floor)
       END DO
       field_norm = SUM(n_n_vec**2)
       new_norm = SUM((n_n_vec - n_n_prev)**2)
       IF (field_norm > 1.d-30) THEN
          rel_res = SQRT(new_norm / field_norm)
       ELSE
          rel_res = 0.d0
       END IF
       IF (rel_res < fv_gs_tol) EXIT
    END DO

    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [NEUTRAL SOLVER] GS iterations: ", I6, "  residual: ", E12.5)') iter, rel_res
       IF (rel_res >= fv_gs_tol) THEN
          WRITE (6, '("   [NEUTRAL SOLVER] WARNING: Gauss-Seidel did not converge in ", I6, " sweeps")') fv_gs_max_iter
       END IF
    END IF

    DO iel = 1, v_mesh%ncells
       nn_old = v_mesh%cells(iel)%U(1)
       ! Under-relaxing the first solve against the seed would make no sense:
       ! the seed is not a previous state, it is just a filler value. An
       ! explicit flag is used rather than a magnitude test on nn_old, so that
       ! the behaviour does not depend on the gap between two constants.
       IF (nn_never_solved) nn_old = n_n_vec(iel)
       n_n_relaxed = numer%picard_relax * n_n_vec(iel) + (1.d0 - numer%picard_relax) * nn_old
       n_n_relaxed = MAX(n_n_relaxed, nn_floor)
       v_mesh%cells(iel)%U(1) = n_n_relaxed
       ! Neutral thermal speed at T_n = 3 eV (Franck-Condon), non-dimensionalized
       ! with the convention of the code: v_adim = SQRT(Mref * T_adim), the same
       ! one as compute_neutral_free_streaming_speed (physics.f90). Inert in the
       ! diffusive model - nobody reads U(2:3) there - but it matters as soon as
       ! a convective flux is connected.
       vt_n = SQRT(phys%Mref * 3.d0 / simpar%refval_temperature)
       v_mesh%cells(iel)%U(2) = n_n_relaxed * vt_n * v_mesh%cells(iel)%normal_wall(1)
       v_mesh%cells(iel)%U(3) = n_n_relaxed * vt_n * v_mesh%cells(iel)%normal_wall(2)
    END DO
    nn_never_solved = .FALSE.

    ! Neutral inventory INT n_n R dOmega - the signature of the transient.
    ! It should start from ~0 (seed), grow at the pace of recycling, and
    ! saturate towards the steady equilibrium (where the BDF term vanishes:
    ! ktis(1) = sum of ktis(2:), consistency of the scheme).
    ! The mean density makes the number interpretable: the integral alone does
    ! not tell whether n_n is 1e17 or 1e21 m-3. The wall flux closes the global
    ! balance: accumulation = wall - pump - (ionization - recombination).
    IF (MPIvar%glob_id == 0) THEN
       field_norm = SUM(cell_vol * (/ (v_mesh%cells(iel)%U(1), iel = 1, v_mesh%ncells) /))
       WRITE (6, '("   [NEUTRALS] inventory INT n_n R dOmega : ", E22.15)') field_norm
       WRITE (6, '("   [NEUTRALS] mean n_n : ", E12.5, " (adim)  = ", E12.5, " m-3")') &
            &field_norm/SUM(cell_vol), field_norm/SUM(cell_vol)*simpar%refval_density
       WRITE (6, '("   [NEUTRALS] wall flux (recycling+puff) : ", E15.8, &
            &"   pump (coef) : ", E15.8)') SUM(wall_influx), SUM(wall_pump)
    END IF

    ! ----- n_n seen by the plasma = n_n_cell, P0, WITHOUT smoothing -----
    ! A P0->P1 smoothing here (nodal average of n_n, then reinterpolation)
    ! would DESTROY conservation: the plasma would see an n_n(x) varying inside
    ! the element while the FV solver only resolves a constant value, and the
    ! neutral sink could no longer match the plasma source. The FV scheme
    ! represents a cell-wise constant field: that is the field the plasma must
    ! see. Smoothing is used for the HDF5 output only
    ! (save_coupled_neutral_density), where it is cosmetic.
    bal_plasma = 0.d0
    bal_neutral = 0.d0
    DO iel = 1, Mesh%Nelems
       idx_start = (iel - 1) * Npel + 1
       DO i = 1, Npel
          ue_el(i, 1) = ures((idx_start + i - 2) * neq_plasma + 1)
          ue_el(i, 2) = ures((idx_start + i - 2) * neq_plasma + 2)
          ue_el(i, 3) = ures((idx_start + i - 2) * neq_plasma + 3)
          ue_el(i, 4) = ures((idx_start + i - 2) * neq_plasma + 4)
       END DO

       n_n = MAX(v_mesh%cells(iel)%U(1), nn_floor)
       cx = v_mesh%cells(iel)%center(1)
       cy = v_mesh%cells(iel)%center(2)

       ! same quadrature as the moments above: this is what closes the balance
       Xel = Mesh%X(Mesh%T(iel, :), :)
       xyg = MATMUL(refElPol%N2D, Xel)
       ueg = MATMUL(refElPol%N2D, ue_el)
       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12

       ! what the neutrals lose in this cell: mom_iz*n_n + off_iz - mom_rec,
       ! exactly the coefficients fed into diag_vec and rhs_vec (MUSCL included)
       bal_neutral = bal_neutral + mom_iz(iel)*n_n + off_iz(iel) - mom_rec(iel)

       DO g = 1, refElPol%NGauss2d
          U_g(1:4) = ueg(g, 1:4)

          ! Atomic sources from compute_neutral_sources (physics.f90), with the
          ! quasilinear convention of the monolithic model: the HDG solver
          ! consumes these arrays as  Auu += (.)*NNi  and  rhs -= (.)*Ni ,
          ! exactly like  Auu += Sn*NNi  and  rhs -= Sn0*Ni  on the monolithic
          ! side. The n_n column is frozen (Picard):
          !   Sn0(i) <- Sn0(i) + Sn(i,inn)*n_n,  column excluded from the
          ! Jacobian.
          ! MUSCL: n_n reconstructed linearly at the Gauss point (the SAME one
          ! as in the off_iz sink above -> the plasma gains exactly what the
          ! neutrals lose). grad_nn = 0 if the switch is off => n_n constant.
          nn_g = n_n
          IF (switch%neutral_muscl) THEN
             nn_g = MAX(n_n + grad_nn(1,iel)*(xyg(g,1) - cx) &
                  &        + grad_nn(2,iel)*(xyg(g,2) - cy), nn_floor)
          END IF
          U5_g(1:4) = U_g(1:4)
          U5_g(5) = nn_g
          CALL compute_neutral_sources(U5_g, Sn5, Sn05)
          neutral_source_gauss(:, iel, g) = Sn05(1:4) + Sn5(1:4, 5)*nn_g
          neutral_jac_source_gauss(:, :, iel, g) = Sn5(1:4, 1:4)
          neutral_density_gauss(iel, g) = nn_g
          neutral_vx_gauss(iel, g) = v_mesh%cells(iel)%U(2) / MAX(n_n, nn_floor)
          neutral_vy_gauss(iel, g) = v_mesh%cells(iel)%U(3) / MAX(n_n, nn_floor)

          ! what the plasma gains: S_1 = -(Sn0(1) + Sn(1,:).U), an exact
          ! identity of the quasilinear convention (Sn = -dS/dU,
          ! Sn0 = -S + (dS/dU).U). Integrated with the same quadrature.
          dvolu = refElPol%gauss_weights2D(g) * detJ(g)
          IF (switch%axisym) dvolu = dvolu * xyg(g, 1)
          bal_plasma = bal_plasma - (Sn05(1) + DOT_PRODUCT(Sn5(1, :), U5_g)) * dvolu
       END DO
    END DO

    ! ----- Conservation check: INT S_plasma + INT S_neutral = 0 -----
    ! Exact BY CONSTRUCTION as soon as (a) both sides integrate the same
    ! quantity with the same quadrature and (b) n_n is P0. It is therefore not
    ! a physics test but a plumbing test: if it does not close, one of the two
    ! assumptions is violated somewhere.
    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [PARTICLE BALANCE] plasma gains : ", E22.15)') bal_plasma
       WRITE (6, '("   [PARTICLE BALANCE] neutrals lose : ", E22.15)') bal_neutral
       IF (ABS(bal_neutral) > 1.d-30) THEN
          WRITE (6, '("   [PARTICLE BALANCE] relative error : ", E12.5)') &
               &ABS(bal_plasma - bal_neutral)/ABS(bal_neutral)
       END IF
    END IF

    DEALLOCATE(ue_el)
    DEALLOCATE(diag_vec, rhs_vec, Dnn, n_n_vec, n_n_prev)
    DEALLOCATE(mom_iz, mom_rec, cell_vol, cn_cell, coeff_face, ueg)
    DEALLOCATE(grad_nn, off_iz)

    phys%idx_rhon_eq = idx_rhon_save
#else
    PRINT *, 'ERROR: Diffusion neutral coupling requires a build with neutral physics'
    PRINT *, '       (compile with MODE=VENUS or a NEUTRAL model).'
    STOP
#endif
  END SUBROUTINE neutral_diffusion_solve_step

  ! ------------------------------------------------------------------
  ! Isothermal 2-moment Euler solver for the neutrals.
  !
  ! Per-cell state U = [n_n, Gamma_x, Gamma_y] (venus_cell_t%U), pressure
  ! p_n = c_n^2 * n_n with c_n^2 = Mref*Ti_limited - which is EXACTLY the
  ! 0.5*a*Ti prefactor of compute_Dnn (Mref = 0.5*phys%a,
  ! adimensionalization.f90:121-124). In conservative form, the strongly
  ! collisional limit of the momentum equation gives
  !   Gamma -> n_n*u_i + [rec*u_i - grad(c_n^2 n_n)] / (nu_fric + nu_iz)
  ! with nu_fric + nu_iz = THE denominator of compute_Dnn, evaluated on the
  ! same cell-average state: the Fick coefficient of that limit is the
  ! diffusive D BY CONSTRUCTION.
  !
  ! Scheme: explicit first-order HLL flux + implicit local sources (IMEX),
  ! marched to steady state with a LOCAL time step (each cell has its own CFL).
  ! The neutrals are "slaved" to the plasma: solved to THEIR steady state at
  ! every Picard call. No Jacobian, no global solver - the source update is
  ! analytic per cell, and the collisional lock Gamma -> n_n*u_i comes out of
  ! the implicit scheme with no tuning (denominator 1 + dt*nu).
  !
  ! Accepted differences with the diffusive model:
  !   - n_n*u_i drift: CX friction targets the poloidal ion velocity
  !     u_i = u_par*b_(R,Z); the diffusive model has no advection term;
  !   - flux -grad(c^2 n)/nu instead of -(c^2/nu)*grad(n): the n*grad(Ti)
  !     term is carried by the conservative form;
  !   - no softplus clipping to [diff_nn_min, diff_nn] on the effective
  !     coefficient: where the diffusive model clips D, the Euler model
  !     transports freely (that is the POINT of a two-moment model).
  !
  ! Wall: mirror state -> the HLL mass flux is EXACTLY zero through the wall
  ! (symmetric states: S_L = -S_R, zero density jump), and
  ! recycling/puff/pump enter as a VOLUMETRIC source in the boundary cells
  ! through wall_influx/wall_pump - the same geometry as in the diffusive
  ! model, hence an identical mass balance between the two models.
  ! The recycled neutral is injected WITHOUT momentum: friction and pressure
  ! gradient set its velocity within one relaxation length.
  ! ------------------------------------------------------------------
  SUBROUTINE neutral_euler_solve_step(ures)
    REAL*8, INTENT(IN) :: ures(:)
#if defined(NEUTRAL) || defined(VENUS)
    INTEGER :: iel, i, g, j, Npel, idx_start, neq_plasma, idx_rhon_save
    INTEGER :: iter, ifa, cell_L, cell_R
    REAL*8 :: U_cell(4), U5_cell(5), U5_g(5), U_g(4)
    REAL*8 :: nn_old, n_n, upar, bnorm
    REAL*8 :: sviz_c, svcx_c, svnn_c, sviz_g, svrec_g, dvolu, r_face
    REAL*8 :: Xel(Mesh%Nnodesperelem, 2)
    REAL*8 :: xyg(refElPol%Ngauss2d, 2)
    REAL*8 :: J11(refElPol%Ngauss2d), J12(refElPol%Ngauss2d)
    REAL*8 :: J21(refElPol%Ngauss2d), J22(refElPol%Ngauss2d)
    REAL*8 :: detJ(refElPol%Ngauss2d)
    REAL*8 :: Bnod(Mesh%Nnodesperelem, 3), bmean(3)
    REAL*8, ALLOCATABLE :: ue_el(:, :), ueg(:, :)
    REAL*8, ALLOCATABLE :: mom_iz(:), mom_rec(:), cell_vol(:)
    REAL*8, ALLOCATABLE :: cn2(:), nu_fric(:), nu_izm(:), ui_cell(:, :)
    REAL*8, ALLOCATABLE :: Ueu(:, :), Unew(:, :), resid(:, :), wsum(:)
    REAL*8, ALLOCATABLE :: face_w(:)
    REAL*8 :: UL3(3), UR3(3), FL(3), FR(3), Fh(3), nhat(2)
    REAL*8 :: nL, nR, cL, cR, unL, unR, pL, pR, SL, SR, smax, w
    REAL*8 :: GnL, uxL, uyL, uxR, uyR
    REAL*8 :: dt_loc, dtvol, denom_g, gs_x, gs_y, gmag, gmax
    REAL*8 :: rel_res, res_num, res_den
    REAL*8 :: n_new, gx_new, gy_new
    REAL*8 :: mach_max, mach
    REAL*8 :: U_relaxed(3), U_prev(3)
    REAL*8 :: Sn5(5, 5), Sn05(5)
    REAL*8 :: bal_plasma, bal_neutral, field_norm

    IF (MPIvar%glob_id == 0) THEN
       PRINT *, '   [NEUTRAL SOLVER] Solving 2D Finite Volumes Neutral Euler (2 moments)...'
    END IF

    neq_plasma = phys%Neq

    ! Wall sources. BEFORE the idx_rhon_eq bracket (setLocalDiff, see the
    ! equivalent comment in the diffusive solver).
    CALL neutral_coupling_wall_moments(ures)

    ! Same idx_rhon_eq = 5 bracket as the diffusive solver: the compute_*
    ! routines of physics.f90 read U(phys%idx_rhon_eq).
    idx_rhon_save = phys%idx_rhon_eq
    phys%idx_rhon_eq = neq_plasma + 1
    Npel = Mesh%Nnodesperelem
    ALLOCATE(ue_el(Npel, neq_plasma))
    ALLOCATE(ueg(refElPol%Ngauss2d, neq_plasma))
    ALLOCATE(mom_iz(v_mesh%ncells), mom_rec(v_mesh%ncells), cell_vol(v_mesh%ncells))
    ALLOCATE(cn2(v_mesh%ncells), nu_fric(v_mesh%ncells), nu_izm(v_mesh%ncells))
    ALLOCATE(ui_cell(2, v_mesh%ncells))
    ALLOCATE(Ueu(3, v_mesh%ncells), Unew(3, v_mesh%ncells))
    ALLOCATE(resid(3, v_mesh%ncells), wsum(v_mesh%ncells))
    ALLOCATE(face_w(v_mesh%nfaces))

    ! ----- Coefficients FROZEN during the marching (Picard lag, exactly as
    ! Dnn is frozen at the incoming n_n in the diffusive solver) -----
    DO iel = 1, Mesh%Nelems
       idx_start = (iel - 1)*Npel + 1
       DO i = 1, Npel
          DO j = 1, neq_plasma
             ue_el(i, j) = ures((idx_start + i - 2)*neq_plasma + j)
          END DO
       END DO
       DO j = 1, neq_plasma
          U_cell(j) = SUM(ue_el(:, j))/Npel
       END DO
       nn_old = MAX(v_mesh%cells(iel)%U(1), nn_floor)
       U5_cell(1:4) = U_cell(1:4)
       U5_cell(5) = nn_old

       ! c_n^2 = Mref*Ti_limited: the prefactor of compute_Dnn (see header).
       CALL compute_neutral_free_streaming_speed(U5_cell, cL)
       cn2(iel) = cL*cL

       ! nu_fric + nu_izm = denominator of compute_Dnn
       ! (physics.f90:813-824: n_i*(sviz + svcx) + n_n*svnn), evaluated on THE
       ! SAME state U5_cell on which the diffusive solver evaluates Dnn. Split
       ! into friction (CX + neutral-neutral, target u_i) and ionization sink
       ! (target 0): their SUM is the coefficient of the diffusive limit.
       CALL compute_sigmaviz(U5_cell, sviz_c)
       CALL compute_sigmavcx(U5_cell, svcx_c)
       CALL compute_sigmavnn(U5_cell, svnn_c)
       nu_fric(iel) = U_cell(1)*svcx_c + nn_old*svnn_c
       nu_izm(iel) = U_cell(1)*sviz_c

       ! Friction target: poloidal u_i = u_par * b_(R,Z).
       Bnod = phys%B(Mesh%T(iel, :), :)
       bmean(1) = SUM(Bnod(:, 1))/Npel
       bmean(2) = SUM(Bnod(:, 2))/Npel
       bmean(3) = SUM(Bnod(:, 3))/Npel
       bnorm = MAX(SQRT(bmean(1)**2 + bmean(2)**2 + bmean(3)**2), 1.d-30)
       upar = U_cell(2)/MAX(U_cell(1), 1.d-30)
       ui_cell(1, iel) = upar*bmean(1)/bnorm
       ui_cell(2, iel) = upar*bmean(2)/bnorm

       ! ----- Moments on the SAME quadrature as the HDG solver - an exact copy
       ! of the diffusive solver, this is the condition for conservation -----
       Xel = Mesh%X(Mesh%T(iel, :), :)
       xyg = MATMUL(refElPol%N2D, Xel)
       ueg = MATMUL(refElPol%N2D, ue_el)
       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12

       mom_iz(iel) = 0.d0
       mom_rec(iel) = 0.d0
       cell_vol(iel) = 0.d0
       DO g = 1, refElPol%Ngauss2d
          dvolu = refElPol%gauss_weights2D(g)*detJ(g)
          IF (switch%axisym) dvolu = dvolu*xyg(g, 1)
          U5_g(1:4) = ueg(g, 1:4)
          U5_g(5) = nn_old
          CALL compute_sigmaviz(U5_g, sviz_g)
          CALL compute_sigmavrec(U5_g, svrec_g)
          mom_iz(iel) = mom_iz(iel) + ueg(g, 1)*sviz_g*dvolu
          mom_rec(iel) = mom_rec(iel) + ueg(g, 1)**2*svrec_g*dvolu
          cell_vol(iel) = cell_vol(iel) + dvolu
       END DO
    END DO

    ! Geometric face weight: L_f * R_f, the SAME axisymmetric measure as the
    ! face coefficient of the diffusive solver (coeff_face).
    DO ifa = 1, v_mesh%nfaces
       r_face = 1.d0
       IF (switch%axisym) r_face = v_mesh%faces(ifa)%midpoint(1)
       face_w(ifa) = v_mesh%faces(ifa)%length*r_face
    END DO

    ! Working state: the state of the previous Picard iterate (warm start),
    ! with the density above the floor by construction.
    DO iel = 1, v_mesh%ncells
       Ueu(1, iel) = MAX(v_mesh%cells(iel)%U(1), nn_floor)
       Ueu(2, iel) = v_mesh%cells(iel)%U(2)
       Ueu(3, iel) = v_mesh%cells(iel)%U(3)
    END DO

    ! ----- Local pseudo-time marching towards the steady state -----
    rel_res = 1.d0
    DO iter = 1, euler_max_iter
       resid = 0.d0
       wsum = 0.d0
       DO ifa = 1, v_mesh%nfaces
          cell_L = v_mesh%faces(ifa)%cell_L
          cell_R = v_mesh%faces(ifa)%cell_R
          nhat = v_mesh%faces(ifa)%normal
          UL3 = Ueu(:, cell_L)
          cL = SQRT(cn2(cell_L))
          IF (cell_R > 0) THEN
             UR3 = Ueu(:, cell_R)
             cR = SQRT(cn2(cell_R))
          ELSE
             ! Wall: mirror state (HLL mass flux exactly zero).
             GnL = UL3(2)*nhat(1) + UL3(3)*nhat(2)
             UR3(1) = UL3(1)
             UR3(2) = UL3(2) - 2.d0*GnL*nhat(1)
             UR3(3) = UL3(3) - 2.d0*GnL*nhat(2)
             cR = cL
          END IF

          nL = MAX(UL3(1), nn_floor)
          nR = MAX(UR3(1), nn_floor)
          uxL = UL3(2)/nL
          uyL = UL3(3)/nL
          uxR = UR3(2)/nR
          uyR = UR3(3)/nR
          unL = uxL*nhat(1) + uyL*nhat(2)
          unR = uxR*nhat(1) + uyR*nhat(2)
          pL = cn2(cell_L)*nL
          IF (cell_R > 0) THEN
             pR = cn2(cell_R)*nR
          ELSE
             pR = pL
          END IF

          ! Davis wave-speed bounds, then the standard HLL flux.
          SL = MIN(unL - cL, unR - cR)
          SR = MAX(unL + cL, unR + cR)

          FL(1) = nL*unL
          FL(2) = UL3(2)*unL + pL*nhat(1)
          FL(3) = UL3(3)*unL + pL*nhat(2)
          FR(1) = nR*unR
          FR(2) = UR3(2)*unR + pR*nhat(1)
          FR(3) = UR3(3)*unR + pR*nhat(2)

          IF (SL >= 0.d0) THEN
             Fh = FL
          ELSE IF (SR <= 0.d0) THEN
             Fh = FR
          ELSE
             Fh = (SR*FL - SL*FR + SL*SR*(UR3 - UL3))/(SR - SL)
          END IF
          smax = MAX(ABS(SL), ABS(SR))

          w = face_w(ifa)
          resid(:, cell_L) = resid(:, cell_L) - Fh*w
          wsum(cell_L) = wsum(cell_L) + smax*w
          IF (cell_R > 0) THEN
             resid(:, cell_R) = resid(:, cell_R) + Fh*w
             wsum(cell_R) = wsum(cell_R) + smax*w
          END IF
       END DO

       res_num = 0.d0
       res_den = 0.d0
       DO iel = 1, v_mesh%ncells
          ! Axisymmetric geometric term: (div T)_R carries -T_theta,theta/R
          ! = -p/R, i.e. +INT_K p dA (PLANE measure, without R) on the
          ! right-hand side of Gamma_R. The Z component has no curvature term.
          IF (switch%axisym) THEN
             resid(2, iel) = resid(2, iel) + cn2(iel)*Ueu(1, iel)*v_mesh%cells(iel)%area
          END IF

          dt_loc = euler_cfl*cell_vol(iel)/MAX(wsum(iel), 1.d-30)
          dtvol = dt_loc/cell_vol(iel)

          ! Continuity: the SAME integrated moments as the diffusive solver
          ! (mom_iz/mom_rec/wall_*), terms linear in n_n taken implicitly -
          ! unconditional positivity (denominator > 1, numerator > 0).
          n_new = (Ueu(1, iel) + dtvol*(resid(1, iel) + mom_rec(iel) + wall_influx(iel))) &
               & /(1.d0 + dtvol*(mom_iz(iel) + wall_pump(iel)))
          n_new = MAX(n_new, nn_floor)

          ! Momentum: friction + ionization + pump taken implicitly, analytic
          ! per cell (no Jacobian). The friction target uses n_new (local
          ! Gauss-Seidel): when nu -> inf, Gamma -> n_n*u_i within the step
          ! itself, which is the collisional lock.
          gs_x = Ueu(2, iel) + dtvol*resid(2, iel)
          gs_y = Ueu(3, iel) + dtvol*resid(3, iel)
          denom_g = 1.d0 + dt_loc*(nu_fric(iel) + nu_izm(iel)) + dtvol*wall_pump(iel)
          gx_new = (gs_x + dt_loc*nu_fric(iel)*n_new*ui_cell(1, iel) &
               & + dtvol*mom_rec(iel)*ui_cell(1, iel))/denom_g
          gy_new = (gs_y + dt_loc*nu_fric(iel)*n_new*ui_cell(2, iel) &
               & + dtvol*mom_rec(iel)*ui_cell(2, iel))/denom_g

          ! Near-vacuum safeguard (see euler_umax_cn).
          gmax = euler_umax_cn*SQRT(cn2(iel))*n_new
          gmag = SQRT(gx_new**2 + gy_new**2)
          IF (gmag > gmax) THEN
             gx_new = gx_new*gmax/gmag
             gy_new = gy_new*gmax/gmag
          END IF

          ! Homogeneous norm: Gamma weighted by 1/c_n to be comparable to n.
          res_num = res_num + (n_new - Ueu(1, iel))**2 &
               & + ((gx_new - Ueu(2, iel))**2 + (gy_new - Ueu(3, iel))**2)/cn2(iel)
          res_den = res_den + n_new**2 + (gx_new**2 + gy_new**2)/cn2(iel)

          Unew(1, iel) = n_new
          Unew(2, iel) = gx_new
          Unew(3, iel) = gy_new
       END DO
       Ueu = Unew

       IF (res_den > 1.d-300) THEN
          rel_res = SQRT(res_num/res_den)
       ELSE
          rel_res = 0.d0
       END IF
       IF (rel_res < euler_tol) EXIT
    END DO

    mach_max = 0.d0
    DO iel = 1, v_mesh%ncells
       mach = SQRT(Ueu(2, iel)**2 + Ueu(3, iel)**2) &
            & /MAX(Ueu(1, iel)*SQRT(cn2(iel)), 1.d-300)
       IF (mach > mach_max) mach_max = mach
    END DO

    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [NEUTRAL SOLVER] Euler pseudo-steps: ", I7, "  residual: ", E12.5, &
            &"  max Mach: ", F9.3)') MIN(iter, euler_max_iter), rel_res, mach_max
       IF (rel_res >= euler_tol) THEN
          WRITE (6, '("   [NEUTRAL SOLVER] WARNING: Euler marching did not converge in ", I7, &
               &" pseudo-steps")') euler_max_iter
       END IF
    END IF

    ! Picard under-relaxation on the FULL state [n, Gamma], same logic and same
    ! nn_never_solved flag as in the diffusive solver.
    DO iel = 1, v_mesh%ncells
       U_prev = v_mesh%cells(iel)%U
       IF (nn_never_solved) U_prev = Ueu(:, iel)
       U_relaxed = numer%picard_relax*Ueu(:, iel) + (1.d0 - numer%picard_relax)*U_prev
       U_relaxed(1) = MAX(U_relaxed(1), nn_floor)
       v_mesh%cells(iel)%U = U_relaxed
    END DO
    nn_never_solved = .FALSE.

    ! Neutral inventory (same diagnostic as the diffusive solver).
    IF (MPIvar%glob_id == 0) THEN
       field_norm = SUM(cell_vol*(/ (v_mesh%cells(iel)%U(1), iel = 1, v_mesh%ncells) /))
       WRITE (6, '("   [NEUTRALS] inventory INT n_n R dOmega : ", E22.15)') field_norm
       WRITE (6, '("   [NEUTRALS] mean n_n : ", E12.5, " (adim)  = ", E12.5, " m-3")') &
            &field_norm/SUM(cell_vol), field_norm/SUM(cell_vol)*simpar%refval_density
       WRITE (6, '("   [NEUTRALS] wall flux (recycling+puff) : ", E15.8, &
            &"   pump (coef) : ", E15.8)') SUM(wall_influx), SUM(wall_pump)
    END IF

    ! ----- Publication to the plasma - an exact copy of the diffusive solver
    ! (n_n P0, no smoothing, same quasilinear conventions), with ONE
    ! difference: the published velocity is the resolved Euler velocity, not
    ! the thermal placeholder of the diffusive solver. -----
    bal_plasma = 0.d0
    bal_neutral = 0.d0
    DO iel = 1, Mesh%Nelems
       idx_start = (iel - 1)*Npel + 1
       DO i = 1, Npel
          DO j = 1, neq_plasma
             ue_el(i, j) = ures((idx_start + i - 2)*neq_plasma + j)
          END DO
       END DO

       n_n = MAX(v_mesh%cells(iel)%U(1), nn_floor)

       Xel = Mesh%X(Mesh%T(iel, :), :)
       xyg = MATMUL(refElPol%N2D, Xel)
       ueg = MATMUL(refElPol%N2D, ue_el)
       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12

       bal_neutral = bal_neutral + mom_iz(iel)*n_n - mom_rec(iel)

       DO g = 1, refElPol%NGauss2d
          U_g(1:4) = ueg(g, 1:4)
          U5_g(1:4) = U_g(1:4)
          U5_g(5) = n_n
          CALL compute_neutral_sources(U5_g, Sn5, Sn05)
          neutral_source_gauss(:, iel, g) = Sn05(1:4) + Sn5(1:4, 5)*n_n
          neutral_jac_source_gauss(:, :, iel, g) = Sn5(1:4, 1:4)
          neutral_density_gauss(iel, g) = n_n
          neutral_vx_gauss(iel, g) = v_mesh%cells(iel)%U(2)/n_n
          neutral_vy_gauss(iel, g) = v_mesh%cells(iel)%U(3)/n_n

          dvolu = refElPol%gauss_weights2D(g)*detJ(g)
          IF (switch%axisym) dvolu = dvolu*xyg(g, 1)
          bal_plasma = bal_plasma - (Sn05(1) + DOT_PRODUCT(Sn5(1, :), U5_g))*dvolu
       END DO
    END DO

    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [PARTICLE BALANCE] plasma gains : ", E22.15)') bal_plasma
       WRITE (6, '("   [PARTICLE BALANCE] neutrals lose : ", E22.15)') bal_neutral
       IF (ABS(bal_neutral) > 1.d-30) THEN
          WRITE (6, '("   [PARTICLE BALANCE] relative error : ", E12.5)') &
               &ABS(bal_plasma - bal_neutral)/ABS(bal_neutral)
       END IF
    END IF

    DEALLOCATE(ue_el, ueg)
    DEALLOCATE(mom_iz, mom_rec, cell_vol)
    DEALLOCATE(cn2, nu_fric, nu_izm, ui_cell)
    DEALLOCATE(Ueu, Unew, resid, wsum, face_w)

    phys%idx_rhon_eq = idx_rhon_save
#else
    PRINT *, 'ERROR: the Euler neutral coupling requires a build with neutral physics'
    PRINT *, '       (compile with MODE=VENUS or a NEUTRAL model).'
    STOP
#endif
  END SUBROUTINE neutral_euler_solve_step

#if defined(NEUTRAL) || defined(VENUS)
  SUBROUTINE neutral_coupling_wall_moments(ures)
    ! Neutral wall source = recycling of the ACTUAL PLASMA FLUX, at the 1D
    ! Gauss points of every exterior face - a line by line replica of
    ! assembly_bohm_bc (hdg_BC.f90:2256-2268):
    !   parallel  : Re * Gamma*bn                     (bn = b.n, the incidence
    !               angle of the field on the wall)
    !   diffusive : -Re * (D_iso*(Qpr.n)_1 - D_ani*bn*(Qpr.b)_1)
    !   pinch     : 0 here - transport_1d is refused at initialization.
    ! Puff and pump are integrated with the SAME quadrature (dline, axisym R).
    !
    ! MUST be called OUTSIDE the idx_rhon_eq = 5 bracket of solve_step:
    ! setLocalDiff writes d_iso(inn,inn,:) under the guard inn > 0, and inn = 5
    ! would overflow the arrays dimensioned with Neq = 4.
    REAL*8, INTENT(IN) :: ures(:)
    INTEGER :: ie, iel, ifl, g, k, c, fl, bc_type, Fi, Npfl, neq_l, Npel_l
    INTEGER :: nod(refElPol%Nfacenodes)
    REAL*8 :: Xf(refElPol%Nfacenodes, 2), Bfl(refElPol%Nfacenodes, 3)
    REAL*8 :: Bmod_nod(refElPol%Nfacenodes), b_nod(refElPol%Nfacenodes, 3)
    REAL*8 :: uef(refElPol%Nfacenodes, phys%Neq)
    REAL*8 :: qf(refElPol%Nfacenodes, 2*phys%Neq)
    REAL*8 :: uf_tr(refElPol%Nfacenodes, phys%Neq)
    REAL*8 :: xyg_f(refElPol%Ngauss1d, 2), xyDer(refElPol%Ngauss1d, 2)
    REAL*8 :: uefg(refElPol%Ngauss1d, phys%Neq), qfg(refElPol%Ngauss1d, 2*phys%Neq)
    REAL*8 :: ufg_tr(refElPol%Ngauss1d, phys%Neq)
    REAL*8 :: bg(refElPol%Ngauss1d, 3)
    REAL*8 :: diff_iso_fac(phys%Neq, phys%Neq, refElPol%Ngauss1d)
    REAL*8 :: diff_ani_fac(phys%Neq, phys%Neq, refElPol%Ngauss1d)
    REAL*8 :: Qpr(2, phys%Neq)
    REAL*8 :: t_g(2), n_g(2), bn, dline, xyDerNorm
    REAL*8 :: recycling_coeff, cryopump_coeff, puff_coeff
    REAL*8 :: par_flux, diff_flux
#ifdef KEQUATION
    REAL*8 :: q_cylfl(refElPol%Nfacenodes), q_cylg(refElPol%Ngauss1d)
#endif

    wall_influx = 0.d0
    wall_pump = 0.d0
    neq_l = phys%Neq
    Npfl = refElPol%Nfacenodes
    Npel_l = Mesh%Nnodesperelem

    DO ie = 1, Mesh%Nextfaces
       fl = Mesh%boundaryFlag(ie)
#ifdef PARALL
       IF (fl .EQ. 0) CYCLE
       IF (Mesh%ghostFaces(Mesh%Nintfaces + ie) .NE. 0) CYCLE
#endif
       bc_type = phys%bcflags(fl)
       ! same coefficients as the monolithic model (RN = 1: apply_trim is
       ! refused at initialization); the DEFAULT branch recycles with Re.
       CALL venus_neutral_wall_bc_coefficients(bc_type, recycling_coeff, cryopump_coeff, puff_coeff)
       IF (switch%neutral_wall_sources_in_elements) THEN
          cryopump_coeff = 0.d0
          puff_coeff = 0.d0
       END IF

       iel = Mesh%extfaces(ie, 1)
       ifl = Mesh%extfaces(ie, 2)
       nod = refElPol%face_nodes(ifl, :)
       Xf = Mesh%X(Mesh%T(iel, nod), :)

       ! element solution at the face nodes (same source as solve_step)
       DO k = 1, Npfl
          DO c = 1, neq_l
             uef(k, c) = ures(((iel - 1)*Npel_l + nod(k) - 1)*neq_l + c)
          END DO
       END DO
       ! gradient at the face nodes - sol%q layout checked against hdg_BC.f90:
       ! ind_qf = (iel-1)*Ndim*Neq*Npel + Neq*Ndim blocks of the face nodes,
       ! component c = dim + (eq-1)*Ndim (dim varying fastest).
       DO k = 1, Npfl
          DO c = 1, 2*neq_l
             qf(k, c) = sol%q((iel - 1)*2*neq_l*Npel_l + 2*neq_l*(nod(k) - 1) + c)
          END DO
       END DO
       ! trace at the face nodes: this is what the monolithic model passes to
       ! setLocalDiff (hdg_BC.f90:670)
       Fi = Mesh%Nintfaces + ie
       DO k = 1, Npfl
          DO c = 1, neq_l
             uf_tr(k, c) = sol%u_tilde((Fi - 1)*neq_l*Npfl + (k - 1)*neq_l + c)
          END DO
       END DO

       Bfl = phys%B(Mesh%T(iel, nod), :)
       Bmod_nod = SQRT(Bfl(:, 1)**2 + Bfl(:, 2)**2 + Bfl(:, 3)**2)
       b_nod(:, 1) = Bfl(:, 1)/Bmod_nod
       b_nod(:, 2) = Bfl(:, 2)/Bmod_nod
       b_nod(:, 3) = Bfl(:, 3)/Bmod_nod

       xyg_f = MATMUL(refElPol%N1D, Xf)
       xyDer = MATMUL(refElPol%Nxi1D, Xf)
       uefg = MATMUL(refElPol%N1D, uef)
       qfg = MATMUL(refElPol%N1D, qf)
       ufg_tr = MATMUL(refElPol%N1D, uf_tr)
       bg = MATMUL(refElPol%N1D, b_nod)

#ifndef KEQUATION
       CALL setLocalDiff(xyg_f, ufg_tr, diff_iso_fac, diff_ani_fac)
#else
       q_cylfl = phys%q_cyl(Mesh%T(iel, nod))
       q_cylg = MATMUL(refElPol%N1D, q_cylfl)
       CALL setLocalDiff(xyg_f, ufg_tr, diff_iso_fac, diff_ani_fac, q_cylg)
#endif

       DO g = 1, refElPol%Ngauss1d
          xyDerNorm = NORM2(xyDer(g, :))
          dline = refElPol%gauss_weights1D(g)*xyDerNorm
          IF (switch%axisym) dline = dline*xyg_f(g, 1)
          t_g = xyDer(g, :)/xyDerNorm
          n_g(1) = t_g(2)
          n_g(2) = -t_g(1)
          bn = DOT_PRODUCT(bg(g, 1:2), n_g)
          Qpr = RESHAPE(qfg(g, :), (/2, neq_l/))

          ! hdg_BC.f90:2256: recycled_parallel_source = uefg(2)*bn
          par_flux = uefg(g, 2)*bn
          ! hdg_BC.f90:2262-2264: recycled_diffusion_source =
          !   -Re*(diffiso(1,1)*(Qpr(:,1).n) - diffani(1,1)*bn*(Qpr(:,1).b))
          diff_flux = -(diff_iso_fac(1, 1, g)*(Qpr(1, 1)*n_g(1) + Qpr(2, 1)*n_g(2)) &
               &      - diff_ani_fac(1, 1, g)*(Qpr(1, 1)*bn*bg(g, 1) + Qpr(2, 1)*bn*bg(g, 2)))

          wall_influx(iel) = wall_influx(iel) &
               & + (recycling_coeff*(par_flux + diff_flux) + puff_coeff)*dline
          wall_pump(iel) = wall_pump(iel) + cryopump_coeff*dline
       END DO
    END DO
  END SUBROUTINE neutral_coupling_wall_moments
#endif

  SUBROUTINE venus_time_integration_coefficients(ktis)
    ! Replica of setTimeIntegrationCoefficients (hdg_ComputeJacobian.f90:2532),
    ! unreachable from here: it is an internal procedure of
    ! HDG_computeJacobian. The literals are copied VERBATIM, including their
    ! single precision (11./6. and not 11.d0/6.d0): "same coefficients as the
    ! plasma" must hold bit for bit, not to within 1e-8.
    ! Orders 4-6 are refused: the monolithic order-4 table is corrupted
    ! (ktis(4) written twice - 4./3. then -0.25 - and ktis(5) never).
    ! Replicating it would copy a bug, departing from it would break parity.
    REAL*8, INTENT(OUT) :: ktis(:)
    INTEGER :: it_ord

    ktis = 0.

    IF (time%ik .LT. time%tis) THEN
       it_ord = time%ik
    ELSE
       it_ord = time%tis
    END IF

    SELECT CASE (it_ord)
    CASE (1)
       ktis(1) = 1.
       ktis(2) = 1.
    CASE (2)
       ktis(1) = 1.5
       ktis(2) = 2
       ktis(3) = -0.5
    CASE (3)
       ktis(1) = 11./6.
       ktis(2) = 3.
       ktis(3) = -1.5
       ktis(4) = 1./3.
    CASE DEFAULT
       PRINT *, 'ERROR: neutral BDF order ', it_ord, ' not supported (tis <= 3).'
       PRINT *, '       (The monolithic order-4 coefficient table is itself corrupted:'
       PRINT *, '        setTimeIntegrationCoefficients writes ktis(4) twice.)'
       STOP
    END SELECT
  END SUBROUTINE venus_time_integration_coefficients

  SUBROUTINE neutral_coupling_picard_residual(u_new, residual)
    REAL*8, INTENT(IN)  :: u_new(:)
    REAL*8, INTENT(OUT) :: residual
    REAL*8 :: diff_norm, old_norm, rel_err
    INTEGER :: i

    residual = 0.d0
    IF (.NOT. coupling_active) RETURN
    IF (.NOT. ALLOCATED(u_plasma_old)) THEN
       residual = HUGE(1.d0)
       RETURN
    END IF

    diff_norm = 0.d0
    old_norm = 0.d0
    DO i = 1, SIZE(u_new)
       diff_norm = diff_norm + (u_new(i) - u_plasma_old(i))**2
       old_norm = old_norm + u_plasma_old(i)**2
    END DO
    IF (old_norm > 1.d-12) THEN
       rel_err = SQRT(diff_norm / old_norm)
    ELSE
       rel_err = SQRT(diff_norm)
    END IF
    residual = rel_err
    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [PICARD RESIDUAL] L2 relative error: ", E12.5)') rel_err
    END IF
  END SUBROUTINE neutral_coupling_picard_residual

  SUBROUTINE neutral_coupling_picard_commit(u_new)
    REAL*8, INTENT(IN) :: u_new(:)
    IF (.NOT. coupling_active) RETURN
    IF (.NOT. ALLOCATED(u_plasma_old)) THEN
       ALLOCATE(u_plasma_old(SIZE(u_new)))
    END IF
    u_plasma_old = u_new
  END SUBROUTINE neutral_coupling_picard_commit

  LOGICAL FUNCTION neutral_coupling_is_converged(u_new, residual)
    REAL*8, INTENT(IN)  :: u_new(:)
    REAL*8, INTENT(OUT) :: residual

    neutral_coupling_is_converged = .TRUE.
    residual = 0.d0
    IF (.NOT. coupling_active) RETURN

    CALL neutral_coupling_picard_residual(u_new, residual)
    IF (residual > numer%picard_tol) neutral_coupling_is_converged = .FALSE.
    CALL neutral_coupling_picard_commit(u_new)
  END FUNCTION neutral_coupling_is_converged

  SUBROUTINE save_coupled_neutral_density(file_id)
    USE HDF5
    USE HDF5_io_module
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: group_id
    INTEGER :: ierr

    IF (.NOT. coupling_active) RETURN
    IF (.NOT. ALLOCATED(v_mesh%cells)) RETURN

    ! Group coupled_neutrals: the neutral FV state is U = [n_n, Gamma_x,
    ! Gamma_y] (venus_cell_t%U). All THREE components are saved, so that the
    ! neutral momentum/flow field resolved by the isothermal Euler model can
    ! be visualized too. The velocity is rebuilt in post: u = Gamma/n_n.
    CALL HDF5_group_create('coupled_neutrals', file_id, group_id, ierr)
    CALL save_cell_component_nodal(group_id, 1, 'rhon')     ! neutral density n_n
    CALL save_cell_component_nodal(group_id, 2, 'Gammax')   ! neutral R momentum
    CALL save_cell_component_nodal(group_id, 3, 'Gammay')   ! neutral Z momentum
    CALL HDF5_group_close(group_id, ierr)
  END SUBROUTINE save_coupled_neutral_density

  ! Projects component `comp` of the neutral FV state (v_mesh%cells%U) onto the
  ! high-order nodes of the HDG mesh and writes it into `group_id` as `name`.
  ! Nodal average (P1), then P1 interpolation on the element nodes: the same
  ! treatment for every component.
  SUBROUTINE save_cell_component_nodal(group_id, comp, name)
    USE HDF5
    USE HDF5_io_module
    INTEGER(HID_T), INTENT(IN) :: group_id
    INTEGER, INTENT(IN) :: comp
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: iel, i, gnode, idx
    REAL*8 :: xi, eta
    REAL*8, ALLOCATABLE :: v_glob(:), node_weight(:), v_array(:)

    ALLOCATE(v_glob(Mesh%Nnodes), node_weight(Mesh%Nnodes))
    v_glob = 0.d0
    node_weight = 0.d0
    DO iel = 1, Mesh%Nelems
       DO i = 1, 3
          gnode = Mesh%T(iel, i)
          v_glob(gnode) = v_glob(gnode) + v_mesh%cells(iel)%U(comp)
          node_weight(gnode) = node_weight(gnode) + 1.d0
       END DO
    END DO
    DO i = 1, Mesh%Nnodes
       IF (node_weight(i) > 0.5d0) v_glob(i) = v_glob(i) / node_weight(i)
    END DO

    ALLOCATE(v_array(Mesh%Nelems * Mesh%Nnodesperelem))
    DO iel = 1, Mesh%Nelems
       DO i = 1, Mesh%Nnodesperelem
          idx = (iel - 1) * Mesh%Nnodesperelem + i
          xi = refElPol%coord2d(i, 1)
          eta = refElPol%coord2d(i, 2)
          v_array(idx) = -0.5d0 * (xi + eta) * v_glob(Mesh%T(iel, 1)) + &
                          0.5d0 * (1.d0 + xi) * v_glob(Mesh%T(iel, 2)) + &
                          0.5d0 * (1.d0 + eta) * v_glob(Mesh%T(iel, 3))
       END DO
    END DO

    CALL HDF5_array1D_saving(group_id, v_array, SIZE(v_array), name)
    DEALLOCATE(v_glob, node_weight, v_array)
  END SUBROUTINE save_cell_component_nodal

  ! ------------------------------------------------------------------
  ! Manufactured-solution (MMS) check of the neutral FV scheme.
  !
  ! Purpose: measure the DISCRETIZATION error of the two-point flux (TPFA)
  ! alone, on the real mesh, with no plasma, no Picard loop, no relaxation and
  ! no reference code. An analytic solution n*(R,Z) is injected, and the
  ! deviation n_h - n* is the error of the scheme by itself.
  !
  ! The monolithic model has no equivalent check for its neutrals - they live
  ! in an HDG discretization that is verified as a whole. The Venus FV scheme
  ! is a new discretization; it must prove its accuracy on its own.
  !
  ! Expected behaviour: TPFA is only consistent on a K-orthogonal mesh. On an
  ! arbitrary triangular mesh the error does NOT vanish under refinement - it
  ! plateaus at the level of mesh non-orthogonality. The number to read is the
  ! absolute error and its trend under refinement, not a clean convergence
  ! order.
  !
  ! Operator under test = EXACTLY the production one
  ! (neutral_diffusion_solve_step), with constant D:
  !   coeff_f = D * L_f * R_f / d_centers        (interior faces)
  ! Boundary: Dirichlet = exact n* at the face midpoint (ghost term
  ! D*R_f*L_f/dist(center,face) on the diagonal + ...*n*(face) on the rhs).
  ! That makes the system SPD and non-singular, and tests the INTERIOR
  ! operator, which is the object of the test. The quasi-Neumann production BC
  ! is a source BC (wall recycling), not what is tested here.
  !
  ! Source consistent with the quadrature of the code: the discrete operator
  ! approximates  INT_K [-div(D R grad n)] dA  (PLANE measure, R as a
  ! coefficient), hence  rhs(K) = sum_g w_g detJ_g * [ -D (R_g lap n* +
  ! dn*/dR) ].
  !
  ! Trigger: environment variable VENUS_MMS=1 (inactive otherwise). The test
  ! runs at initialization, prints, and STOPs - it does not perturb runs.
  ! ------------------------------------------------------------------
  SUBROUTINE venus_run_mms_selftest()
#if defined(NEUTRAL) || defined(VENUS)
    CHARACTER(LEN=8) :: env_val
    INTEGER :: env_len, env_stat
    INTEGER :: iel, ifa, g, Npel, cell_L, cell_R, iface_loc, neighbor, iter, imode
    REAL*8  :: Rmin, Rmax, Zmin, Zmax, LR, LZ, kR, kZ, A0, A1, Dconst, nwave, pi
    REAL*8  :: val, dvdR, lap, r_face, dx, dy, d_centers, d_bf, T_bf, fsrc, Rg
    REAL*8  :: dvolu, sum_offdiag, rel_res, field_norm, dnorm
    REAL*8  :: l2num, l2den, linf_num, linf_den, e, nstar_c, resid, bnorm, Aii_n
    LOGICAL :: patch
    REAL*8, ALLOCATABLE :: diag(:), coeff(:), rhs(:), nvec(:), nprev(:), cvol(:)
    REAL*8 :: Xel(Mesh%Nnodesperelem, 2)
    REAL*8 :: xyg(refElPol%Ngauss2d, 2)
    REAL*8 :: J11(refElPol%Ngauss2d), J12(refElPol%Ngauss2d)
    REAL*8 :: J21(refElPol%Ngauss2d), J22(refElPol%Ngauss2d)
    REAL*8 :: detJ(refElPol%Ngauss2d)
    INTEGER, PARAMETER :: mms_gs_max = 500000
    REAL*8,  PARAMETER :: mms_gs_tol = 1.d-11

    CALL GET_ENVIRONMENT_VARIABLE('VENUS_MMS', env_val, env_len, env_stat)
    IF (env_stat /= 0) RETURN                       ! variable not set: inactive
    IF (TRIM(ADJUSTL(env_val)) /= '1') RETURN

    pi = 4.d0*ATAN(1.d0)
    Dconst = 1.d0
    A0 = 1.d0
    A1 = 0.5d0                                        ! n* in [0.5, 1.5] > 0
    nwave = 2.d0                                      ! ~2 wavelengths per domain
    Npel = Mesh%Nnodesperelem

    ! Bounding box from the cell centers (robust with respect to units).
    Rmin =  HUGE(1.d0); Rmax = -HUGE(1.d0)
    Zmin =  HUGE(1.d0); Zmax = -HUGE(1.d0)
    DO iel = 1, v_mesh%ncells
       Rmin = MIN(Rmin, v_mesh%cells(iel)%center(1))
       Rmax = MAX(Rmax, v_mesh%cells(iel)%center(1))
       Zmin = MIN(Zmin, v_mesh%cells(iel)%center(2))
       Zmax = MAX(Zmax, v_mesh%cells(iel)%center(2))
    END DO
    LR = MAX(Rmax - Rmin, 1.d-30)
    LZ = MAX(Zmax - Zmin, 1.d-30)
    kR = 2.d0*pi*nwave/LR
    kZ = 2.d0*pi*nwave/LZ

    ALLOCATE(diag(v_mesh%ncells), rhs(v_mesh%ncells), nvec(v_mesh%ncells))
    ALLOCATE(nprev(v_mesh%ncells), cvol(v_mesh%ncells), coeff(v_mesh%nfaces))

    ! ----- Geometry (field independent): cvol, coeff_face, diag -----
    ! cvol(K) = INT_K R dA (axisym) - the measure of the error norm.
    cvol = 0.d0
    DO iel = 1, Mesh%Nelems
       Xel = Mesh%X(Mesh%T(iel, :), :)
       xyg = MATMUL(refElPol%N2D, Xel)
       J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
       J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
       J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
       J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
       detJ = J11*J22 - J21*J12
       DO g = 1, refElPol%Ngauss2d
          dvolu = refElPol%gauss_weights2D(g) * detJ(g)
          IF (switch%axisym) dvolu = dvolu * xyg(g, 1)
          cvol(iel) = cvol(iel) + dvolu
       END DO
    END DO

    coeff = 0.d0
    diag = 0.d0
    DO ifa = 1, v_mesh%nfaces
       cell_L = v_mesh%faces(ifa)%cell_L
       cell_R = v_mesh%faces(ifa)%cell_R
       r_face = 1.d0
       IF (switch%axisym) r_face = v_mesh%faces(ifa)%midpoint(1)
       IF (cell_R > 0) THEN                           ! interior face: TPFA
          dx = v_mesh%cells(cell_L)%center(1) - v_mesh%cells(cell_R)%center(1)
          dy = v_mesh%cells(cell_L)%center(2) - v_mesh%cells(cell_R)%center(2)
          d_centers = MAX(SQRT(dx**2 + dy**2), 1.d-30)
          coeff(ifa) = Dconst * v_mesh%faces(ifa)%length * r_face / d_centers
          diag(cell_L) = diag(cell_L) + coeff(ifa)
          diag(cell_R) = diag(cell_R) + coeff(ifa)
       ELSE                                           ! boundary face: ghost Dirichlet
          dx = v_mesh%cells(cell_L)%center(1) - v_mesh%faces(ifa)%midpoint(1)
          dy = v_mesh%cells(cell_L)%center(2) - v_mesh%faces(ifa)%midpoint(2)
          d_bf = MAX(SQRT(dx**2 + dy**2), 1.d-30)
          T_bf = Dconst * r_face * v_mesh%faces(ifa)%length / d_bf
          diag(cell_L) = diag(cell_L) + T_bf
       END IF
    END DO

    IF (MPIvar%glob_id == 0) THEN
       WRITE (6, '("   [VENUS MMS] Manufactured solution - isotropic diffusion, D = 1")')
       WRITE (6, '("   [VENUS MMS] mesh: ", I7, " cells, ", I8, " faces (", I7, " on boundary)")') &
            & v_mesh%ncells, v_mesh%nfaces, v_mesh%nfaces_ext
       WRITE (6, '("   [VENUS MMS] R in [", E12.5, ",", E12.5, "]  Z in [", E12.5, ",", E12.5, "]")') &
            & Rmin, Rmax, Zmin, Zmax
    END IF

    ! ----- Two modes: 1 = sinusoid (the measurement), 2 = constant patch -----
    DO imode = 1, 2
       patch = (imode == 2)

       ! rhs = source integral (quadrature of the code) + boundary Dirichlet
       rhs = 0.d0
       DO iel = 1, Mesh%Nelems
          Xel = Mesh%X(Mesh%T(iel, :), :)
          xyg = MATMUL(refElPol%N2D, Xel)
          J11 = MATMUL(refElPol%Nxi2D, Xel(:, 1))
          J12 = MATMUL(refElPol%Nxi2D, Xel(:, 2))
          J21 = MATMUL(refElPol%Neta2D, Xel(:, 1))
          J22 = MATMUL(refElPol%Neta2D, Xel(:, 2))
          detJ = J11*J22 - J21*J12
          DO g = 1, refElPol%Ngauss2d
             CALL mms_field(xyg(g, 1), xyg(g, 2), val, dvdR, lap)
             IF (switch%axisym) THEN
                Rg = xyg(g, 1)
                fsrc = -Dconst * (Rg*lap + dvdR)      ! -div(D R grad n*)
             ELSE
                fsrc = -Dconst * lap                  ! -div(D grad n*)
             END IF
             rhs(iel) = rhs(iel) + refElPol%gauss_weights2D(g) * detJ(g) * fsrc
          END DO
       END DO
       DO ifa = 1, v_mesh%nfaces
          IF (v_mesh%faces(ifa)%cell_R > 0) CYCLE
          cell_L = v_mesh%faces(ifa)%cell_L
          r_face = 1.d0
          IF (switch%axisym) r_face = v_mesh%faces(ifa)%midpoint(1)
          dx = v_mesh%cells(cell_L)%center(1) - v_mesh%faces(ifa)%midpoint(1)
          dy = v_mesh%cells(cell_L)%center(2) - v_mesh%faces(ifa)%midpoint(2)
          d_bf = MAX(SQRT(dx**2 + dy**2), 1.d-30)
          T_bf = Dconst * r_face * v_mesh%faces(ifa)%length / d_bf
          CALL mms_field(v_mesh%faces(ifa)%midpoint(1), v_mesh%faces(ifa)%midpoint(2), val, dvdR, lap)
          rhs(cell_L) = rhs(cell_L) + T_bf * val
       END DO

       ! Starting point = n* at the cell centers: GS then only has to polish
       ! the discretization gap (a few %), it does not invent the solution.
       DO iel = 1, v_mesh%ncells
          CALL mms_field(v_mesh%cells(iel)%center(1), v_mesh%cells(iel)%center(2), val, dvdR, lap)
          nvec(iel) = val
       END DO

       rel_res = 1.d0
       DO iter = 1, mms_gs_max
          nprev = nvec
          DO iel = 1, v_mesh%ncells
             sum_offdiag = 0.d0
             DO iface_loc = 1, 3
                ifa = v_mesh%cell_faces(iface_loc, iel)
                cell_L = v_mesh%faces(ifa)%cell_L
                cell_R = v_mesh%faces(ifa)%cell_R
                IF (cell_L == iel) THEN
                   neighbor = cell_R
                ELSE
                   neighbor = cell_L
                END IF
                IF (neighbor > 0) sum_offdiag = sum_offdiag + coeff(ifa) * nvec(neighbor)
             END DO
             nvec(iel) = (rhs(iel) + sum_offdiag) / diag(iel)
          END DO
          field_norm = SUM(nvec**2)
          dnorm = SUM((nvec - nprev)**2)
          IF (field_norm > 1.d-30) THEN
             rel_res = SQRT(dnorm / field_norm)
          ELSE
             rel_res = 0.d0
          END IF
          IF (rel_res < mms_gs_tol) EXIT
       END DO

       ! TRUE residual ||A n - b|| / ||b||: certifies that the linear solver
       ! does not pollute the measurement of the discretization error.
       resid = 0.d0
       bnorm = 0.d0
       DO iel = 1, v_mesh%ncells
          sum_offdiag = 0.d0
          DO iface_loc = 1, 3
             ifa = v_mesh%cell_faces(iface_loc, iel)
             cell_L = v_mesh%faces(ifa)%cell_L
             cell_R = v_mesh%faces(ifa)%cell_R
             IF (cell_L == iel) THEN
                neighbor = cell_R
             ELSE
                neighbor = cell_L
             END IF
             IF (neighbor > 0) sum_offdiag = sum_offdiag + coeff(ifa) * nvec(neighbor)
          END DO
          Aii_n = diag(iel)*nvec(iel) - sum_offdiag
          resid = resid + (rhs(iel) - Aii_n)**2
          bnorm = bnorm + rhs(iel)**2
       END DO
       resid = SQRT(resid / MAX(bnorm, 1.d-300))

       ! Error against the exact solution, weighted by the axisymmetric volume.
       l2num = 0.d0; l2den = 0.d0; linf_num = 0.d0; linf_den = 0.d0
       DO iel = 1, v_mesh%ncells
          CALL mms_field(v_mesh%cells(iel)%center(1), v_mesh%cells(iel)%center(2), nstar_c, dvdR, lap)
          e = nvec(iel) - nstar_c
          l2num = l2num + cvol(iel) * e**2
          l2den = l2den + cvol(iel) * nstar_c**2
          linf_num = MAX(linf_num, ABS(e))
          linf_den = MAX(linf_den, ABS(nstar_c))
       END DO

       IF (MPIvar%glob_id == 0) THEN
          IF (.NOT. patch) THEN
             WRITE (6, '("   [VENUS MMS] --- sinusoidal field (kR=", E11.4, " kZ=", E11.4, ") ---")') kR, kZ
          ELSE
             WRITE (6, '("   [VENUS MMS] --- patch test (constant field, expected ~1e-14) ---")')
          END IF
          WRITE (6, '("   [VENUS MMS]   GS: ", I7, " sweeps, linear residual ", E12.5)') iter, resid
          WRITE (6, '("   [VENUS MMS]   relative L2 error (volume weighted): ", E12.5)') SQRT(l2num/MAX(l2den, 1.d-300))
          WRITE (6, '("   [VENUS MMS]   relative Linf error                : ", E12.5)') linf_num/MAX(linf_den, 1.d-300)
          IF (resid > 1.d-8) &
               & WRITE (6, '("   [VENUS MMS]   WARNING: GS not converged - the error above is polluted")')
          IF (patch .AND. SQRT(l2num/MAX(l2den, 1.d-300)) > 1.d-10) &
               & WRITE (6, '("   [VENUS MMS]   PATCH FAILURE: the constant field is not reproduced - assembly bug")')
       END IF
    END DO

    DEALLOCATE(diag, rhs, nvec, nprev, cvol, coeff)
    IF (MPIvar%glob_id == 0) WRITE (6, '("   [VENUS MMS] done - stopping (VENUS_MMS).")')
    STOP

  CONTAINS

    ! Manufactured solution and its analytic derivatives. patch (host) switches
    ! to the constant field of the patch test.
    SUBROUTINE mms_field(R, Z, fval, fdR, flap)
      REAL*8, INTENT(IN)  :: R, Z
      REAL*8, INTENT(OUT) :: fval, fdR, flap
      REAL*8 :: aR, aZ, sR, cR, cZ
      IF (patch) THEN
         fval = 1.d0; fdR = 0.d0; flap = 0.d0
         RETURN
      END IF
      aR = kR*(R - Rmin); aZ = kZ*(Z - Zmin)
      sR = SIN(aR); cR = COS(aR); cZ = COS(aZ)
      fval = A0 + A1*sR*cZ
      fdR  = A1*kR*cR*cZ
      flap = -A1*(kR**2 + kZ**2)*sR*cZ
    END SUBROUTINE mms_field
#else
    PRINT *, 'ERROR: venus_run_mms_selftest requires a build with neutral physics.'
    STOP
#endif
  END SUBROUTINE venus_run_mms_selftest

END MODULE neutral_coupling
