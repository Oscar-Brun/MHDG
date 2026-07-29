!***********************************************************************
! project: MHDG
! file: venus_mesh.f90
! description: Mesh builder for the Venus 2D Finite Volume solver.
!***********************************************************************
MODULE venus_mesh
  USE globals
  USE venus_types
  IMPLICIT NONE
  SAVE

CONTAINS

  ! This routine draws the big map of cells and faces from the high-order HDG mesh.
  SUBROUTINE venus_mesh_build(v_mesh)
    TYPE(venus_mesh_t), INTENT(OUT) :: v_mesh
    INTEGER :: iel, ifa, ie, ii, ifa_loc, nodeA, nodeB, bc
    REAL*8, DIMENSION(2) :: XA, XB, V, normal_raw
    REAL*8 :: norm_val, sum_area

    PRINT *, '--- Venus: Building Finite Volume Mesh ---'

    ! 1. Allocate arrays and set counts
    v_mesh%ncells = Mesh%Nelems
    v_mesh%nfaces = Mesh%Nintfaces + Mesh%Nextfaces
    v_mesh%nfaces_int = Mesh%Nintfaces
    v_mesh%nfaces_ext = Mesh%Nextfaces

    ALLOCATE(v_mesh%cells(v_mesh%ncells))
    ALLOCATE(v_mesh%faces(v_mesh%nfaces))
    ALLOCATE(v_mesh%cell_faces(3, v_mesh%ncells))

    ! 2. Initialize cells
    DO iel = 1, v_mesh%ncells
       v_mesh%cells(iel)%area = Mesh%elemSize(iel)
       
       ! Center of the triangle is the average of its 3 corners
       v_mesh%cells(iel)%center(1) = (Mesh%X(Mesh%Tlin(iel,1), 1) + &
                                      Mesh%X(Mesh%Tlin(iel,2), 1) + &
                                      Mesh%X(Mesh%Tlin(iel,3), 1)) / 3.d0
       v_mesh%cells(iel)%center(2) = (Mesh%X(Mesh%Tlin(iel,1), 2) + &
                                      Mesh%X(Mesh%Tlin(iel,2), 2) + &
                                      Mesh%X(Mesh%Tlin(iel,3), 2)) / 3.d0
       
       ! Associate the 3 faces forming this triangle
       v_mesh%cell_faces(1:3, iel) = Mesh%F(iel, 1:3)
       
       ! Initialize variables
       v_mesh%cells(iel)%U = 0.d0
       v_mesh%cells(iel)%U_old = 0.d0
       v_mesh%cells(iel)%grad_U = 0.d0
       v_mesh%cells(iel)%Ti = 0.d0
       v_mesh%cells(iel)%ne = 0.d0
       v_mesh%cells(iel)%te = 0.d0
       v_mesh%cells(iel)%u_par = 0.d0
       v_mesh%cells(iel)%S = 0.d0
       v_mesh%cells(iel)%Jac = 0.d0
    END DO

    ! 3. Process interior faces
    DO ii = 1, Mesh%Nintfaces
       ifa = ii
       iel = Mesh%intFaces(ii, 1)     ! Left cell
       ifa_loc = Mesh%intFaces(ii, 2) ! Local face index on the left cell (1, 2, or 3)
       
       v_mesh%faces(ifa)%cell_L = iel
       v_mesh%faces(ifa)%cell_R = Mesh%intFaces(ii, 3) ! Right cell
       v_mesh%faces(ifa)%boundary_type = 0            ! Internal face
       v_mesh%faces(ifa)%bc_flag = 0
       
       ! Get the two end-nodes of this face from the left cell geometry
       IF (ifa_loc == 1) THEN
          nodeA = Mesh%Tlin(iel, 1)
          nodeB = Mesh%Tlin(iel, 2)
       ELSEIF (ifa_loc == 2) THEN
          nodeA = Mesh%Tlin(iel, 2)
          nodeB = Mesh%Tlin(iel, 3)
       ELSE
          nodeA = Mesh%Tlin(iel, 3)
          nodeB = Mesh%Tlin(iel, 1)
       END IF
       
       XA = Mesh%X(nodeA, :)
       XB = Mesh%X(nodeB, :)
       
       ! Compute center and length
       v_mesh%faces(ifa)%midpoint = 0.5d0 * (XA + XB)
       v_mesh%faces(ifa)%length = SQRT((XA(1) - XB(1))**2 + (XA(2) - XB(2))**2)
       
       ! Compute normal pointing from Left to Right cell
       V = v_mesh%faces(ifa)%midpoint - v_mesh%cells(iel)%center
       normal_raw(1) = XB(2) - XA(2)
       normal_raw(2) = -(XB(1) - XA(1))
       IF (DOT_PRODUCT(normal_raw, V) < 0.d0) normal_raw = -normal_raw
       norm_val = SQRT(normal_raw(1)**2 + normal_raw(2)**2)
       v_mesh%faces(ifa)%normal = normal_raw / norm_val
    END DO

    ! 4. Process boundary faces
    DO ie = 1, Mesh%Nextfaces
       ifa = ie + Mesh%Nintfaces
       iel = Mesh%extFaces(ie, 1)     ! Left cell (adjacent cell)
       ifa_loc = Mesh%extFaces(ie, 2) ! Local face index on the cell (1, 2, or 3)
       bc = phys%bcflags(Mesh%boundaryFlag(ie)) ! Physical boundary flag
       
       v_mesh%faces(ifa)%cell_L = iel
       v_mesh%faces(ifa)%cell_R = -1   ! External boundary
       v_mesh%faces(ifa)%boundary_type = bc
       v_mesh%faces(ifa)%bc_flag = Mesh%boundaryFlag(ie)
       
       ! Get the two end-nodes of this face
       IF (ifa_loc == 1) THEN
          nodeA = Mesh%Tlin(iel, 1)
          nodeB = Mesh%Tlin(iel, 2)
       ELSEIF (ifa_loc == 2) THEN
          nodeA = Mesh%Tlin(iel, 2)
          nodeB = Mesh%Tlin(iel, 3)
       ELSE
          nodeA = Mesh%Tlin(iel, 3)
          nodeB = Mesh%Tlin(iel, 1)
       END IF
       
       XA = Mesh%X(nodeA, :)
       XB = Mesh%X(nodeB, :)
       
       ! Compute center and length
       v_mesh%faces(ifa)%midpoint = 0.5d0 * (XA + XB)
       v_mesh%faces(ifa)%length = SQRT((XA(1) - XB(1))**2 + (XA(2) - XB(2))**2)
       
       ! Compute normal pointing outwards (from Cell to outside)
       V = v_mesh%faces(ifa)%midpoint - v_mesh%cells(iel)%center
       normal_raw(1) = XB(2) - XA(2)
       normal_raw(2) = -(XB(1) - XA(1))
       IF (DOT_PRODUCT(normal_raw, V) < 0.d0) normal_raw = -normal_raw
       norm_val = SQRT(normal_raw(1)**2 + normal_raw(2)**2)
       v_mesh%faces(ifa)%normal = normal_raw / norm_val
    END DO

    ! 5. Compute distance to nearest wall face for each cell
    DO iel = 1, v_mesh%ncells
       v_mesh%cells(iel)%dist_wall = 1.d10
       v_mesh%cells(iel)%normal_wall = 0.d0
       DO ie = 1, v_mesh%nfaces_ext
          ifa = ie + v_mesh%nfaces_int
          ! Calculate distance from cell center to boundary face midpoint
          V = v_mesh%cells(iel)%center - v_mesh%faces(ifa)%midpoint
          norm_val = SQRT(V(1)**2 + V(2)**2)
          IF (norm_val < v_mesh%cells(iel)%dist_wall) THEN
             v_mesh%cells(iel)%dist_wall = norm_val
             ! Inward unit normal (pointing into the domain)
             v_mesh%cells(iel)%normal_wall = -v_mesh%faces(ifa)%normal
          END IF
       END DO
    END DO

    ! Validation check: sum of cell areas must match total mesh area
    sum_area = 0.d0
    DO iel = 1, v_mesh%ncells
       sum_area = sum_area + v_mesh%cells(iel)%area
    END DO
    PRINT *, '   -> Number of cells: ', v_mesh%ncells
    PRINT *, '   -> Number of faces: ', v_mesh%nfaces
    PRINT *, '   -> Total domain area check: ', sum_area
    PRINT *, '--- Venus: FV Mesh successfully built ---'

  END SUBROUTINE venus_mesh_build

END MODULE venus_mesh
