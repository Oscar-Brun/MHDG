!*****************************************
! project: MHDG
! file: venus_types.f90
! date: 24/06/2026
! Define data structures for Venus FV solver
!*****************************************
MODULE venus_types
  IMPLICIT NONE

  ! Structure representing a Finite Volume cell (a triangle from HDG)
  TYPE :: venus_cell_t
    REAL*8  :: center(2)       ! Cell center of gravity (xc, yc)
    REAL*8  :: area            ! Cell area (triangle surface)
    REAL*8  :: U(3)            ! Conservative variables: [rho_n, Gamma_nx, Gamma_ny]
    REAL*8  :: U_old(3)        ! Conservative variables at previous time step (for temporal scheme)
    REAL*8  :: grad_U(2,3)     ! Reconstructed gradients for MUSCL order 2: [d/dx, d/dy] for each variable
    REAL*8  :: Ti              ! Plasma ion temperature in this cell (from plasma)
    REAL*8  :: ne              ! Plasma electron density in this cell (from plasma)
    REAL*8  :: te              ! Plasma electron temperature in this cell (from plasma)
    REAL*8  :: u_par           ! Plasma parallel velocity in this cell (from plasma)
    REAL*8  :: S(3)            ! Cell source terms: [S_rho, S_gammax, S_gammay]
    REAL*8  :: Jac(3,3)        ! Cell local Jacobian of the source terms dS/dU
    REAL*8  :: dist_wall       ! Distance to the nearest wall face
    REAL*8  :: normal_wall(2)  ! Normal vector pointing from the nearest wall face to this cell
  END TYPE venus_cell_t

  ! Structure representing a Finite Volume face (an edge of a triangle)
  TYPE :: venus_face_t
    INTEGER :: cell_L          ! Left cell index (always >= 1)
    INTEGER :: cell_R          ! Right cell index (>= 1, or -1 if boundary face)
    REAL*8  :: normal(2)       ! Unit normal vector (pointing from cell_L to cell_R)
    REAL*8  :: length          ! Length of the face (edge length)
    REAL*8  :: midpoint(2)     ! Midpoint coordinates (xf, yf)
    INTEGER :: boundary_type   ! Type of boundary condition if cell_R == -1 (recycling, pump, puff, reflection)
    INTEGER :: bc_flag         ! Raw boundary flag index from the mesh
  END TYPE venus_face_t

  ! Structure representing the entire Finite Volume mesh
  TYPE :: venus_mesh_t
    INTEGER :: ncells          ! Number of cells (= number of elements in HDG)
    INTEGER :: nfaces          ! Total number of faces
    INTEGER :: nfaces_int      ! Number of interior faces
    INTEGER :: nfaces_ext      ! Number of exterior/boundary faces
    
    TYPE(venus_cell_t), ALLOCATABLE :: cells(:)
    TYPE(venus_face_t), ALLOCATABLE :: faces(:)
    
    INTEGER, ALLOCATABLE :: cell_faces(:,:) ! For each cell, indices of its 3 faces (size: 3, ncells)
  END TYPE venus_mesh_t

END MODULE venus_types
