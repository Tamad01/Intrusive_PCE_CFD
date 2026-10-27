program tcheby_1d_pce
  !=====================================================================
  ! Rotating Rayleigh-Benard / PCE INTRUSIF (P=1, Hermite)
  ! 2 PARAMETRES INCERTAINS (Q=2) -> 3 MODES :
  !   Mode 0 = moyenne
  !   Mode 1 = sensibilite a alpha   (couplage alpha*w)   ~ xi1
  !   Mode 2 = sensibilite a la CL   (temperature paroi)  ~ xi2
  !=====================================================================

  use mpi
  use decomp_2d
  use ifport

  use m_mesh_base
  use m_numerics
  use m_operator_base
  use m_operator_tcheby
  use m_operator_fourier_dft
  use m_fourier_transform
  use m_tensor_product
  use m_hdf5_ifce

  use m_solver_diag_cart_hhi
  use m_solver_diag_cart_hhi_cplx
  use m_navier_stokes_cart
  use m_snapshots
  use m_quadrature

  implicit none

  integer, parameter :: OX = 1, OY = 2, OZ = 3
  integer, parameter :: SCHEME_EXPLICIT_CLASSICAL  = 1
  integer, parameter :: SCHEME_IMPLICIT_CLASSICAL  = 2
  integer, parameter :: SCHEME_IMPLICIT_OLSHANSKII = 3
  integer, parameter :: SCHEME_FROM_LEGACY_INPUT   = 0
  integer, parameter :: line_max = 100

  ! Amplitude d'incertitude sur la temperature de paroi chaude (mode 2).
  ! DOIT etre un parameter (constante) car utilise dans une fonction PURE (udf_bc_hot).
  ! Pour changer la valeur, modifier ici et recompiler.
  real(dp), parameter :: BC_HOT_AMP_PARAM = 0.05_dp

  type(t_mesh_base) :: msh(3)

  type(t_operator_fourier_dft) :: opx, opy
  type(t_operator_tcheby)      :: opz

  ! === PCE : Triplement des solveurs (Mode 0, Mode 1=alpha, Mode 2=paroi) ===
  type(t_solver_diag_cart_hhi)      :: eqn_u0, eqn_v0, eqn_w0, eqn_t0, eqn_p0
  type(t_solver_diag_cart_hhi)      :: eqn_u1, eqn_v1, eqn_w1, eqn_t1, eqn_p1
  type(t_solver_diag_cart_hhi)      :: eqn_u2, eqn_v2, eqn_w2, eqn_t2, eqn_p2
  type(t_solver_diag_cart_hhi_cplx) :: eqn_uv0, eqn_uv1, eqn_uv2
  type(t_quadrature) :: quad

  type(DECOMP_INFO) :: ph

  ! Workspaces partages
  real(dp), allocatable :: w01(:,:,:), w02(:,:,:), w03(:,:,:)
  real(dp), allocatable :: w04(:,:,:), w05(:,:,:), w06(:,:,:)
  real(dp), allocatable :: w07(:,:,:), w08(:,:,:), w09(:,:,:)

  real(dp), allocatable :: dxp0(:,:,:), dyp0(:,:,:), dzp0(:,:,:)
  real(dp), allocatable :: dxp1(:,:,:), dyp1(:,:,:), dzp1(:,:,:)
  real(dp), allocatable :: dxp2(:,:,:), dyp2(:,:,:), dzp2(:,:,:)

  ! === MODE 0 (Moyenne) ===
  real(dp), allocatable :: u0(:,:,:), v0(:,:,:), w0(:,:,:)
  real(dp), allocatable :: pres0(:,:,:), temp0(:,:,:), phi0(:,:,:)
  real(dp), allocatable :: rhs_u0(:,:,:), rhs_v0(:,:,:), rhs_w0(:,:,:)
  real(dp), allocatable :: rhs_t0(:,:,:), rhs_phi0(:,:,:)
  real(dp), allocatable :: nl_u0(:,:,:), nl_v0(:,:,:), nl_w0(:,:,:), nl_t0(:,:,:)
  real(dp), allocatable :: nl_u0_old(:,:,:), nl_v0_old(:,:,:), nl_w0_old(:,:,:), nl_t0_old(:,:,:)
  real(dp), allocatable :: lin_u0(:,:,:), lin_v0(:,:,:), lin_w0(:,:,:), lin_t0(:,:,:)
  real(dp), allocatable :: u0_tilde(:,:,:), v0_tilde(:,:,:), w0_tilde(:,:,:)
  real(dp), allocatable :: u0_old(:,:,:), v0_old(:,:,:)
  real(dp), allocatable :: u0_save(:,:,:), v0_save(:,:,:)

  ! === MODE 1 (sensibilite alpha) ===
  real(dp), allocatable :: u1(:,:,:), v1(:,:,:), w1(:,:,:)
  real(dp), allocatable :: pres1(:,:,:), temp1(:,:,:), phi1(:,:,:)
  real(dp), allocatable :: rhs_u1(:,:,:), rhs_v1(:,:,:), rhs_w1(:,:,:)
  real(dp), allocatable :: rhs_t1(:,:,:), rhs_phi1(:,:,:)
  real(dp), allocatable :: nl_u1(:,:,:), nl_v1(:,:,:), nl_w1(:,:,:), nl_t1(:,:,:)
  real(dp), allocatable :: nl_u1_old(:,:,:), nl_v1_old(:,:,:), nl_w1_old(:,:,:), nl_t1_old(:,:,:)
  real(dp), allocatable :: lin_u1(:,:,:), lin_v1(:,:,:), lin_w1(:,:,:), lin_t1(:,:,:)
  real(dp), allocatable :: u1_tilde(:,:,:), v1_tilde(:,:,:), w1_tilde(:,:,:)
  real(dp), allocatable :: u1_old(:,:,:), v1_old(:,:,:)
  real(dp), allocatable :: u1_save(:,:,:), v1_save(:,:,:)

  ! === MODE 2 (sensibilite CL paroi chaude) ===
  real(dp), allocatable :: u2(:,:,:), v2(:,:,:), w2(:,:,:)
  real(dp), allocatable :: pres2(:,:,:), temp2(:,:,:), phi2(:,:,:)
  real(dp), allocatable :: rhs_u2(:,:,:), rhs_v2(:,:,:), rhs_w2(:,:,:)
  real(dp), allocatable :: rhs_t2(:,:,:), rhs_phi2(:,:,:)
  real(dp), allocatable :: nl_u2(:,:,:), nl_v2(:,:,:), nl_w2(:,:,:), nl_t2(:,:,:)
  real(dp), allocatable :: nl_u2_old(:,:,:), nl_v2_old(:,:,:), nl_w2_old(:,:,:), nl_t2_old(:,:,:)
  real(dp), allocatable :: lin_u2(:,:,:), lin_v2(:,:,:), lin_w2(:,:,:), lin_t2(:,:,:)
  real(dp), allocatable :: u2_tilde(:,:,:), v2_tilde(:,:,:), w2_tilde(:,:,:)
  real(dp), allocatable :: u2_old(:,:,:), v2_old(:,:,:)
  real(dp), allocatable :: u2_save(:,:,:), v2_save(:,:,:)

  ! === Temporaires termes croises (3 modes) ===
  ! base_k = NL(mode k seul) ; sum01 = NL(0+1) ; sum02 = NL(0+2)
  real(dp), allocatable :: u_a(:,:,:), v_a(:,:,:), w_a(:,:,:), t_a(:,:,:)
  real(dp), allocatable :: base_t0(:,:,:), base_t1(:,:,:), base_t2(:,:,:)
  real(dp), allocatable :: sum01_t(:,:,:), sum02_t(:,:,:)
  real(dp), allocatable :: base_u0(:,:,:), base_u1(:,:,:), base_u2(:,:,:)
  real(dp), allocatable :: sum01_u(:,:,:), sum02_u(:,:,:)
  real(dp), allocatable :: base_v0(:,:,:), base_v1(:,:,:), base_v2(:,:,:)
  real(dp), allocatable :: sum01_v(:,:,:), sum02_v(:,:,:)
  real(dp), allocatable :: base_w0(:,:,:), base_w1(:,:,:), base_w2(:,:,:)
  real(dp), allocatable :: sum01_w(:,:,:), sum02_w(:,:,:)

  real(dp), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)

  real(dp) :: geostrophic_error
  real(dp) :: ageostrophic_error
  real(dp) :: splitting_error

  logical :: first_step_ab = .true.

  integer :: nx, ny, nz
  integer :: nb_cpu_y, nb_cpu_z
  real(dp) :: x_min, x_max, y_min, y_max, z_min, z_max
  real(dp) :: k_mode = -1.0_dp
  namelist /parameters_cube/ nx, ny, nz, x_max, y_max, z_max, nb_cpu_y, nb_cpu_z, k_mode

  real(dp) :: pr, ra, ek
  namelist /parameters_physical/ pr, ra, ek

  real(dp) :: tmax, dt, cfl
  integer  :: dn_dump = 1000
  integer  :: scheme_type = SCHEME_FROM_LEGACY_INPUT
  logical  :: resume = .false., explicit = .false.
  character(len=1024) :: root_dir
  real(dp) :: dt_min     = 1.0e-8_dp
  real(dp) :: dt_max     = 1.0e0_dp
  real(dp) :: cfl_target = 0.1_dp
  namelist /parameters_timescheme/ tmax, dt, dt_min, dt_max, cfl, cfl_target, dn_dump, root_dir, resume, &
       explicit, scheme_type

  logical :: resume_snap = .false.
  integer :: dn_snap = 1000
  character(len=1024) :: snap_dir
  namelist /parameters_diagnostics/ dn_snap, snap_dir, resume_snap

  integer :: ierr, rank, iostat
  integer :: i, j, k, it_time
  integer :: n(3), cpu_grid(2)
  integer :: snap_id

  real(dp) :: xmin(3), xmax(3)

  ! === PCE : Variables incertitude ===
  real(dp) :: nu, omega, omega_tilde
  real(dp) :: alpha0, alpha1
  real(dp) :: kappa_diff
  real(dp) :: t_cur, snap_time
  real(dp) :: wall_start, step_walltime
  real(dp) :: div_max, kinetic_w, std_w
  real(dp) :: var_alpha, var_paroi, sobol_alpha, sobol_paroi

  character(len=1024) :: input_file
  character(len=1024) :: file_dump, base_snap, filename
  character(len=6)    :: snap_num

  call mpi_init(ierr)
  call h5open_f(ierr)
  call mpi_comm_rank(mpi_comm_world, rank, ierr)

  call command_line_read_input(input_file)
  call read_input_file(input_file)
  call broadcast_input_parameters()

  call setup_physical_parameters()
  call setup_decomposition_and_operators()
  call setup_equations()
  call allocate_fields()
  call initialize_fields()
  call initialize_snapshot_output()

  t_cur = 0.0_dp
  snap_time = 0.0_dp

  do it_time = 1, 100000
     call advance_one_step_abcn()
     if (t_cur > tmax) then
        call write_final_state_and_stop()
     end if
  end do

  call dump_state('WRITE')
  call mpi_finalize(ierr)

contains

  !===================================================================
  ! Input handling
  !===================================================================

  subroutine read_input_file(file_name)
    character(len=*), intent(in) :: file_name
    integer :: ierr_loc

    if (rank == 0) then
       open(unit=input_unit, file=trim(file_name), status='old', action='read')
       read(input_unit, nml=parameters_cube,        iostat=iostat)
       if (iostat /= 0) write(*,'(a,i0)') 'WARNING parameters_cube, iostat = ', iostat
       read(input_unit, nml=parameters_physical,    iostat=iostat)
       if (iostat /= 0) write(*,'(a,i0)') 'WARNING parameters_physical, iostat = ', iostat
       read(input_unit, nml=parameters_timescheme,  iostat=iostat)
       if (iostat /= 0) write(*,'(a,i0)') 'WARNING parameters_timescheme, iostat = ', iostat
       read(input_unit, nml=parameters_diagnostics, iostat=iostat)
       if (iostat /= 0) write(*,'(a,i0)') 'WARNING parameters_diagnostics, iostat = ', iostat
       close(input_unit)

       call configure_scheme_from_input()
       call validate_scheme_type()

       write(*, nml=parameters_timescheme, iostat=iostat)

       ierr_loc = system('mkdir -p ' // trim(root_dir) // '/dump')
       ierr_loc = system('mkdir -p ' // trim(root_dir) // trim(snap_dir))
    end if
  end subroutine read_input_file

  subroutine broadcast_input_parameters()
    x_min = 0.0_dp
    y_min = 0.0_dp
    z_min = 0.0_dp

    call MPI_BCAST(nx,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ny,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nz,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(x_min,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(x_max,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(y_min,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(y_max,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(z_min,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(z_max,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(k_mode,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (k_mode > 0.0_dp) then
       x_max = 5 * 2.0_dp * pi / k_mode
       y_max = 5 * 2.0_dp * pi / k_mode
    end if

    call MPI_BCAST(nb_cpu_y, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nb_cpu_z, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(pr,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ra,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ek,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call MPI_BCAST(resume,   1,    MPI_LOGICAL,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(root_dir, 1024, MPI_CHARACTER,        0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dn_dump,  1,    MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt,       1,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_min,   1,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dt_max,   1,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cfl,      1,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cfl_target, 1,  MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(scheme_type, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(explicit, 1,    MPI_LOGICAL,          0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(tmax,     1,    MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    explicit = (scheme_type == SCHEME_EXPLICIT_CLASSICAL)

    call MPI_BCAST(dn_snap,     1,    MPI_INTEGER,   0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(snap_dir,    1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(resume_snap, 1,    MPI_LOGICAL,   0, MPI_COMM_WORLD, ierr)
  end subroutine broadcast_input_parameters

  subroutine configure_scheme_from_input()
    if (scheme_type == SCHEME_FROM_LEGACY_INPUT) then
       if (explicit) then
          scheme_type = SCHEME_EXPLICIT_CLASSICAL
       else
          scheme_type = SCHEME_IMPLICIT_OLSHANSKII
       end if
    end if
    explicit = (scheme_type == SCHEME_EXPLICIT_CLASSICAL)
  end subroutine configure_scheme_from_input

  subroutine validate_scheme_type()
    if (scheme_type < SCHEME_EXPLICIT_CLASSICAL .or. &
        scheme_type > SCHEME_IMPLICIT_OLSHANSKII) then
       write(*,'(a,i0)') 'Invalid scheme_type = ', scheme_type
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine validate_scheme_type

  !===================================================================
  ! Setup
  !===================================================================

  subroutine setup_physical_parameters()
    n = [nx, ny, nz]
    cpu_grid = [nb_cpu_y, nb_cpu_z]

    xmin = [x_min, y_min, z_min]
    xmax = [x_max, y_max, z_max]

    file_dump = trim(root_dir) // 'dump/dump.h5'
    base_snap = trim(root_dir) // trim(snap_dir) // 'snap_'

    nu    = sqrt(pr / ra)
    omega = sqrt(pr / ra) / ek

    ! === PCE : incertitudes ===
    alpha0     = 1.0_dp             ! couplage moyen
    alpha1     = 0.05_dp * alpha0   ! incertitude sur alpha (5%) -> mode 1 (xi1)
    kappa_diff = 1.0_dp / sqrt(pr * ra)
    ! incertitude paroi chaude (mode 2) : voir BC_HOT_AMP_PARAM en tete de programme

    call update_omega_tilde()
  end subroutine setup_physical_parameters

  subroutine update_omega_tilde()
    if (scheme_type == SCHEME_EXPLICIT_CLASSICAL) then
       omega_tilde = 0.0_dp
    else
       omega_tilde = omega * dt
    end if
  end subroutine update_omega_tilde

  subroutine setup_decomposition_and_operators()
    call DECOMP_2D_INIT( &
         NX = n(1) + 1, NY = n(2) + 1, NZ = n(3) + 1, &
         P_ROW = cpu_grid(1), P_COL = cpu_grid(2))

    call GET_DECOMP_INFO(ph)

    call START_TENSOR_PRODUCT( &
         ph%xst, ph%xen, n(1) + 1, &
         ph%yst, ph%yen, n(2) + 1, &
         ph%zst, ph%zen, n(3) + 1)

    call START_FOURIER_TRANSFORMS( &
         ph%xst, ph%xen, n(1) + 1, &
         ph%yst, ph%yen, n(2) + 1, &
         ph%zst, ph%zen, n(3) + 1)

    call M_TENSOR_PRODUCT_UT( &
         ph%xst, ph%xen, nx, &
         ph%yst, ph%yen, ny, &
         ph%zst, ph%zen, nz)

    call msh(1)%initialize(xmin(1), xmax(1), n(1), .true.)
    call msh(2)%initialize(xmin(2), xmax(2), n(2), .true.)
    call msh(3)%initialize(xmin(3), xmax(3), n(3))

    call START_OP_BASE()
    call opx%INIT_OPERATOR_FOURIER_DFT(mesh=msh(OX), axis=OX)
    call opy%INIT_OPERATOR_FOURIER_DFT(mesh=msh(OY), axis=OY)
    call opz%INIT_OPERATOR_TCHEBY(mesh=msh(OZ), axis=OZ)
  end subroutine setup_decomposition_and_operators

  subroutine setup_equations()
    real(dp) :: dirichlet(2), neumann(2)
    real(dp) :: nu_momentum(3), nu_energy(3), nu_poisson(3)
    real(dp) :: kappa

    dirichlet = [1.0_dp, 0.0_dp]
    neumann   = [0.0_dp, 1.0_dp]

    nu_momentum = -0.5_dp * nu * [1.0_dp, 1.0_dp, 1.0_dp]

    kappa = 1.0_dp / sqrt(pr * ra)
    nu_energy = -0.5_dp * kappa * [1.0_dp, 1.0_dp, 1.0_dp]
    nu_poisson = [1.0_dp, 1.0_dp, 1.0_dp]

    ! --- MODE 0 ---
    call setup_scalar_hhi_solver(eqn_u0, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_v0, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_w0, nu_momentum, dirichlet, dirichlet, udf_null, udf_null)
    call eqn_uv0%SET_PARAMS(nu=nu_momentum, sigma=0.0_dp, omega=0.0_dp)
    call eqn_uv0%SET_BCS(axis=OZ, bcs_minus=neumann, bcs_plus=neumann)
    call eqn_uv0%INITIALISE(msh, opx, opy, opz, ph)
    call eqn_uv0%SET_BVS(mesh=msh, axis=OZ, udf_minus=udf_null, udf_plus=udf_null)
    call setup_scalar_hhi_solver(eqn_p0, nu_poisson, neumann, neumann, udf_null, udf_null)

    ! --- MODE 1 (sensibilite alpha) ---
    call setup_scalar_hhi_solver(eqn_u1, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_v1, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_w1, nu_momentum, dirichlet, dirichlet, udf_null, udf_null)
    call eqn_uv1%SET_PARAMS(nu=nu_momentum, sigma=0.0_dp, omega=0.0_dp)
    call eqn_uv1%SET_BCS(axis=OZ, bcs_minus=neumann, bcs_plus=neumann)
    call eqn_uv1%INITIALISE(msh, opx, opy, opz, ph)
    call eqn_uv1%SET_BVS(mesh=msh, axis=OZ, udf_minus=udf_null, udf_plus=udf_null)
    call setup_scalar_hhi_solver(eqn_p1, nu_poisson, neumann, neumann, udf_null, udf_null)

    ! --- MODE 2 (sensibilite CL paroi chaude) ---
    call setup_scalar_hhi_solver(eqn_u2, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_v2, nu_momentum, neumann, neumann, udf_null, udf_null)
    call setup_scalar_hhi_solver(eqn_w2, nu_momentum, dirichlet, dirichlet, udf_null, udf_null)
    call eqn_uv2%SET_PARAMS(nu=nu_momentum, sigma=0.0_dp, omega=0.0_dp)
    call eqn_uv2%SET_BCS(axis=OZ, bcs_minus=neumann, bcs_plus=neumann)
    call eqn_uv2%INITIALISE(msh, opx, opy, opz, ph)
    call eqn_uv2%SET_BVS(mesh=msh, axis=OZ, udf_minus=udf_null, udf_plus=udf_null)
    call setup_scalar_hhi_solver(eqn_p2, nu_poisson, neumann, neumann, udf_null, udf_null)

    ! --- Solveurs temperature : meme diffusion deterministe pour les 3 modes ---
    !     Mode 0 et 1 : CL homogene (udf_null)
    !     Mode 2      : CL paroi chaude (bas) INCERTAINE = udf_bc_hot
    call setup_scalar_hhi_solver(eqn_t0, nu_energy, dirichlet, dirichlet, udf_null,   udf_null)
    call setup_scalar_hhi_solver(eqn_t1, nu_energy, dirichlet, dirichlet, udf_null,   udf_null)
    call setup_scalar_hhi_solver(eqn_t2, nu_energy, dirichlet, dirichlet, udf_bc_hot, udf_null)

    ! --- QUADRATURE ---
    call INIT_quadrature_hhi(quad, xmin, xmax, &
         ph%xst, ph%xen, n(1), &
         ph%yst, ph%yen, n(2), &
         ph%zst, ph%zen, n(3))
  end subroutine setup_equations

  function quad_mean(fi) result(mean_val)
    implicit none
    real(dp), allocatable, intent(inout) :: fi(:,:,:)
    real(dp) :: mean_val
    call get_quadrature_hhi(quad, fi, w09, &
         ph%xst, ph%xen, n(1), &
         ph%yst, ph%yen, n(2), &
         ph%zst, ph%zen, n(3))
    mean_val = w09(ph%xst(1), ph%xst(2), ph%xst(3))
  end function quad_mean

  subroutine compute_laplacian(fi, lap_fi)
    real(dp), allocatable, intent(inout) :: fi(:,:,:)
    real(dp), allocatable, intent(inout) :: lap_fi(:,:,:)
    call GRAD(opx, opy, opz, fi, w06, w07, w08)
    call DIV(opx, opy, opz, w06, w07, w08, lap_fi, w01, w02, w03)
  end subroutine compute_laplacian

  subroutine setup_scalar_hhi_solver(eqn, nu_coef, bcs_minus, bcs_plus, udf_minus, udf_plus)
    type(t_solver_diag_cart_hhi), intent(inout) :: eqn
    real(dp), intent(in) :: nu_coef(3)
    real(dp), intent(in) :: bcs_minus(2), bcs_plus(2)
    interface
       real(dp) pure function udf_minus(t_cur, x, y, z)
         import :: dp
         real(dp), intent(in) :: t_cur, x, y, z
       end function udf_minus
       real(dp) pure function udf_plus(t_cur, x, y, z)
         import :: dp
         real(dp), intent(in) :: t_cur, x, y, z
       end function udf_plus
    end interface
    call eqn%SET_PARAMS(nu=nu_coef, sigma=0.0_dp)
    call eqn%SET_BCS(axis=OZ, bcs_minus=bcs_minus, bcs_plus=bcs_plus)
    call eqn%INITIALISE(msh, opx, opy, opz, ph)
    call eqn%SET_BVS(mesh=msh, axis=OZ, udf_minus=udf_minus, udf_plus=udf_plus)
  end subroutine setup_scalar_hhi_solver

  !===================================================================
  ! Allocation
  !===================================================================

  subroutine alloc_one(a)
    real(dp), allocatable, intent(inout) :: a(:,:,:)
    call alloc_x(a, opt_global=.true.); a = 0.0_dp
  end subroutine alloc_one

  subroutine allocate_fields()
    ! Workspaces
    call alloc_one(w01); call alloc_one(w02); call alloc_one(w03)
    call alloc_one(w04); call alloc_one(w05); call alloc_one(w06)
    call alloc_one(w07); call alloc_one(w08); call alloc_one(w09)
    call alloc_one(x);   call alloc_one(y);   call alloc_one(z)

    ! --- MODE 0 ---
    call alloc_one(u0); call alloc_one(v0); call alloc_one(w0)
    call alloc_one(pres0); call alloc_one(temp0); call alloc_one(phi0)
    call alloc_one(dxp0); call alloc_one(dyp0); call alloc_one(dzp0)
    call alloc_one(rhs_u0); call alloc_one(rhs_v0); call alloc_one(rhs_w0)
    call alloc_one(rhs_t0); call alloc_one(rhs_phi0)
    call alloc_one(nl_u0); call alloc_one(nl_v0); call alloc_one(nl_w0); call alloc_one(nl_t0)
    call alloc_one(nl_u0_old); call alloc_one(nl_v0_old); call alloc_one(nl_w0_old); call alloc_one(nl_t0_old)
    call alloc_one(lin_u0); call alloc_one(lin_v0); call alloc_one(lin_w0); call alloc_one(lin_t0)
    call alloc_one(u0_tilde); call alloc_one(v0_tilde); call alloc_one(w0_tilde)
    call alloc_one(u0_old); call alloc_one(v0_old)
    call alloc_one(u0_save); call alloc_one(v0_save)

    ! --- MODE 1 ---
    call alloc_one(u1); call alloc_one(v1); call alloc_one(w1)
    call alloc_one(pres1); call alloc_one(temp1); call alloc_one(phi1)
    call alloc_one(dxp1); call alloc_one(dyp1); call alloc_one(dzp1)
    call alloc_one(rhs_u1); call alloc_one(rhs_v1); call alloc_one(rhs_w1)
    call alloc_one(rhs_t1); call alloc_one(rhs_phi1)
    call alloc_one(nl_u1); call alloc_one(nl_v1); call alloc_one(nl_w1); call alloc_one(nl_t1)
    call alloc_one(nl_u1_old); call alloc_one(nl_v1_old); call alloc_one(nl_w1_old); call alloc_one(nl_t1_old)
    call alloc_one(lin_u1); call alloc_one(lin_v1); call alloc_one(lin_w1); call alloc_one(lin_t1)
    call alloc_one(u1_tilde); call alloc_one(v1_tilde); call alloc_one(w1_tilde)
    call alloc_one(u1_old); call alloc_one(v1_old)
    call alloc_one(u1_save); call alloc_one(v1_save)

    ! --- MODE 2 ---
    call alloc_one(u2); call alloc_one(v2); call alloc_one(w2)
    call alloc_one(pres2); call alloc_one(temp2); call alloc_one(phi2)
    call alloc_one(dxp2); call alloc_one(dyp2); call alloc_one(dzp2)
    call alloc_one(rhs_u2); call alloc_one(rhs_v2); call alloc_one(rhs_w2)
    call alloc_one(rhs_t2); call alloc_one(rhs_phi2)
    call alloc_one(nl_u2); call alloc_one(nl_v2); call alloc_one(nl_w2); call alloc_one(nl_t2)
    call alloc_one(nl_u2_old); call alloc_one(nl_v2_old); call alloc_one(nl_w2_old); call alloc_one(nl_t2_old)
    call alloc_one(lin_u2); call alloc_one(lin_v2); call alloc_one(lin_w2); call alloc_one(lin_t2)
    call alloc_one(u2_tilde); call alloc_one(v2_tilde); call alloc_one(w2_tilde)
    call alloc_one(u2_old); call alloc_one(v2_old)
    call alloc_one(u2_save); call alloc_one(v2_save)

    ! --- TEMPORAIRES termes croises ---
    call alloc_one(u_a); call alloc_one(v_a); call alloc_one(w_a); call alloc_one(t_a)
    call alloc_one(base_t0); call alloc_one(base_t1); call alloc_one(base_t2)
    call alloc_one(sum01_t); call alloc_one(sum02_t)
    call alloc_one(base_u0); call alloc_one(base_u1); call alloc_one(base_u2)
    call alloc_one(sum01_u); call alloc_one(sum02_u)
    call alloc_one(base_v0); call alloc_one(base_v1); call alloc_one(base_v2)
    call alloc_one(sum01_v); call alloc_one(sum02_v)
    call alloc_one(base_w0); call alloc_one(base_w1); call alloc_one(base_w2)
    call alloc_one(sum01_w); call alloc_one(sum02_w)
  end subroutine allocate_fields

  subroutine initialize_fields()
    call build_grid_arrays()
    if (resume) then
       call read_restart()
    else
       call set_initial_condition()
    end if
    ! AB2 history (Coriolis explicite) pour les 3 modes
    u0_old = u0 ; v0_old = v0 ; u0_save = u0 ; v0_save = v0
    u1_old = u1 ; v1_old = v1 ; u1_save = u1 ; v1_save = v1
    u2_old = u2 ; v2_old = v2 ; u2_save = u2 ; v2_save = v2
  end subroutine initialize_fields

  subroutine build_grid_arrays()
    forall(i = ph%xst(1):ph%xen(1), &
           j = ph%xst(2):ph%xen(2), &
           k = ph%xst(3):ph%xen(3))
       x(i,j,k) = msh(1)%x(i)
       y(i,j,k) = msh(2)%x(j)
       z(i,j,k) = msh(3)%x(k)
    end forall
  end subroutine build_grid_arrays

  subroutine read_restart()
    real(dp) :: dumped_dt
    call DUMP_HDF5_BASIC(file_dump, 'READ', t_cur, dumped_dt, msh, u0, v0, w0, pres0, temp0, w09)
  end subroutine read_restart

  subroutine set_initial_condition()
    integer :: is(3), ie(3), ii, jj, kk
    real(dp) :: eps, kx_mode, ky_mode

    t_cur = 0.0_dp
    call dump_state('NEW')
    eps = 1.0e-6_dp

    u0 = 0.0_dp; v0 = 0.0_dp; w0 = 0.0_dp; pres0 = 0.0_dp; phi0 = 0.0_dp; temp0 = 0.0_dp
    u1 = 0.0_dp; v1 = 0.0_dp; w1 = 0.0_dp; pres1 = 0.0_dp; phi1 = 0.0_dp; temp1 = 0.0_dp
    u2 = 0.0_dp; v2 = 0.0_dp; w2 = 0.0_dp; pres2 = 0.0_dp; phi2 = 0.0_dp; temp2 = 0.0_dp

    ! Perturbation initiale sur le mode 0 uniquement
    if (k_mode > 0.0_dp) then
       is = get_is()
       ie = get_ie()
       kx_mode = k_mode
       ky_mode = 0.0_dp
       do kk = is(3), ie(3)
          do jj = is(2), ie(2)
             do ii = is(1), ie(1)
                temp0(ii,jj,kk) = eps * cos(kx_mode*x(ii,jj,kk)) &
                                      * cos(ky_mode*y(ii,jj,kk)) &
                                      * sin(pi*z(ii,jj,kk))
             end do
          end do
       end do
    else
       call random_seed()
       call random_number(w01)
       temp0 = eps * (w01 - 1.0_dp)
    end if
  end subroutine set_initial_condition

  subroutine initialize_snapshot_output()
    if (.not. resume_snap) then
       call EXPORT_GRID(trim(base_snap) // 'grid.h5', x, y, z)
    end if
    filename = trim(base_snap) // 'grid.h5'
    call update_id_and_time(trim(filename), snap_time, snap_id)
    write(snap_num, '(I6.6)') snap_id
  end subroutine initialize_snapshot_output

  !===================================================================
  ! Time stepping: AB2--Crank--Nicolson
  !===================================================================

  subroutine update_dt_cfl(dt_inout, cfl_current, cfl_target_in, dt_min_in, dt_max_in)
    real(dp), intent(inout) :: dt_inout
    real(dp), intent(in)    :: cfl_current, cfl_target_in, dt_min_in, dt_max_in
    real(dp), parameter :: eps_cfl    = 1.0e-14_dp
    real(dp), parameter :: grow_max   = 1.05_dp
    real(dp), parameter :: shrink_max = 0.5_dp
    real(dp) :: dt_candidate
    if (cfl_current > eps_cfl) then
       dt_candidate = dt_inout * cfl_target_in / cfl_current
    else
       dt_candidate = grow_max * dt_inout
    end if
    dt_candidate = min(dt_candidate, grow_max * dt_inout)
    dt_candidate = max(dt_candidate, shrink_max * dt_inout)
    dt_candidate = min(dt_candidate, dt_max_in)
    dt_candidate = max(dt_candidate, dt_min_in)
    dt_inout = dt_candidate
  end subroutine update_dt_cfl

  subroutine advance_one_step_abcn()
    call GetCFL(msh(1), msh(2), msh(3), u0, v0, w0, dt, cfl)
    call update_dt_cfl(dt, cfl, cfl_target, dt_min, dt_max)
    call update_omega_tilde()
    t_cur = t_cur + dt
    snap_time = snap_time + dt
    wall_start = MPI_Wtime()
    call advance_abcn_step()
    call write_diagnostics_and_outputs()
  end subroutine advance_one_step_abcn

  subroutine advance_abcn_step()
    real(dp) :: sigma
    real(dp) :: nu_poisson(3)
    real(dp) :: coef_h, coef_z
    real(dp) :: omega_cn
    real(dp) :: div_before, div_after
    integer :: isb(3), ieb(3), isi(3), iei(3)

    sigma    = 1.0_dp / dt
    omega_cn = 0.5_dp * omega

    !===============================================================
    ! TEMPERATURE : diffusion deterministe, couplage alpha*w incertain (xi1)
    !===============================================================

    ! 1. Advection : on calcule NL pour mode0, mode1, mode2, (0+1), (0+2)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u0, v0, w0, temp0, base_t0, &
         w01, w02, w03, w04, w05, w06, w07, w08, w09)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u1, v1, w1, temp1, base_t1, &
         w01, w02, w03, w04, w05, w06, w07, w08, w09)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u2, v2, w2, temp2, base_t2, &
         w01, w02, w03, w04, w05, w06, w07, w08, w09)
    ! somme 0+1
    u_a = u0 + u1 ; v_a = v0 + v1 ; w_a = w0 + w1 ; t_a = temp0 + temp1
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u_a, v_a, w_a, t_a, sum01_t, &
         w01, w02, w03, w04, w05, w06, w07, w08, w09)
    ! somme 0+2
    u_a = u0 + u2 ; v_a = v0 + v2 ; w_a = w0 + w2 ; t_a = temp0 + temp2
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u_a, v_a, w_a, t_a, sum02_t, &
         w01, w02, w03, w04, w05, w06, w07, w08, w09)

    ! 2. Repartition Galerkin P=1, Q=2 (formules verifiees symboliquement)
    nl_t0 = base_t0 + base_t1 + base_t2
    nl_t1 = sum01_t - base_t0 - base_t1
    nl_t2 = sum02_t - base_t0 - base_t2

    ! 3. Couplage flottabilite alpha*w (incertitude alpha = xi1)
    !    mode0 : alpha0*w0 + alpha1*w1
    !    mode1 : alpha0*w1 + alpha1*w0
    !    mode2 : alpha0*w2            (alpha1 ne couple pas avec mode2 en P=1)
    nl_t0 = nl_t0 - (alpha0*w0 + alpha1*w1)
    nl_t1 = nl_t1 - (alpha0*w1 + alpha1*w0)
    nl_t2 = nl_t2 - (alpha0*w2)

    ! 4. Diffusion thermique deterministe (Crank-Nicolson)
    call compute_laplacian(temp0, lin_t0)
    call compute_laplacian(temp1, lin_t1)
    call compute_laplacian(temp2, lin_t2)

    if (first_step_ab) then
       rhs_t0 = sigma*temp0 + 0.5_dp*kappa_diff*lin_t0 - nl_t0
       rhs_t1 = sigma*temp1 + 0.5_dp*kappa_diff*lin_t1 - nl_t1
       rhs_t2 = sigma*temp2 + 0.5_dp*kappa_diff*lin_t2 - nl_t2
    else
       rhs_t0 = sigma*temp0 + 0.5_dp*kappa_diff*lin_t0 - (1.5_dp*nl_t0 - 0.5_dp*nl_t0_old)
       rhs_t1 = sigma*temp1 + 0.5_dp*kappa_diff*lin_t1 - (1.5_dp*nl_t1 - 0.5_dp*nl_t1_old)
       rhs_t2 = sigma*temp2 + 0.5_dp*kappa_diff*lin_t2 - (1.5_dp*nl_t2 - 0.5_dp*nl_t2_old)
    end if

    call eqn_t0%SOLVE(temp0, rhs_t0, sigma, ph)
    call eqn_t1%SOLVE(temp1, rhs_t1, sigma, ph)
    call eqn_t2%SOLVE(temp2, rhs_t2, sigma, ph)

    !===============================================================
    ! QUANTITE DE MOUVEMENT
    !===============================================================

    call GRAD(opx, opy, opz, pres0, dxp0, dyp0, dzp0)
    call GRAD(opx, opy, opz, pres1, dxp1, dyp1, dzp1)
    call GRAD(opx, opy, opz, pres2, dxp2, dyp2, dzp2)

    ! 1. Advection : NL pour mode0, mode1, mode2, (0+1), (0+2)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u0, v0, w0, base_u0, base_v0, base_w0, &
         w04, w05, w06, w07, w08, w09, w01, w02)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u1, v1, w1, base_u1, base_v1, base_w1, &
         w04, w05, w06, w07, w08, w09, w01, w02)
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u2, v2, w2, base_u2, base_v2, base_w2, &
         w04, w05, w06, w07, w08, w09, w01, w02)
    ! somme 0+1
    u_a = u0 + u1 ; v_a = v0 + v1 ; w_a = w0 + w1
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u_a, v_a, w_a, sum01_u, sum01_v, sum01_w, &
         w04, w05, w06, w07, w08, w09, w01, w02)
    ! somme 0+2
    u_a = u0 + u2 ; v_a = v0 + v2 ; w_a = w0 + w2
    call COMPUTE_NON_LINEAR_TERMS( &
         opx, opy, opz, u_a, v_a, w_a, sum02_u, sum02_v, sum02_w, &
         w04, w05, w06, w07, w08, w09, w01, w02)

    ! 2. Repartition Galerkin
    nl_u0 = base_u0 + base_u1 + base_u2
    nl_v0 = base_v0 + base_v1 + base_v2
    nl_w0 = base_w0 + base_w1 + base_w2

    nl_u1 = sum01_u - base_u0 - base_u1
    nl_v1 = sum01_v - base_v0 - base_v1
    nl_w1 = sum01_w - base_w0 - base_w1

    nl_u2 = sum02_u - base_u0 - base_u2
    nl_v2 = sum02_v - base_v0 - base_v2
    nl_w2 = sum02_w - base_w0 - base_w2

    ! 3. Couplage Boussinesq (deterministe, mode par mode)
    nl_w0 = nl_w0 - temp0
    nl_w1 = nl_w1 - temp1
    nl_w2 = nl_w2 - temp2

    ! 4. Laplaciens cinematiques
    call compute_laplacian(u0, lin_u0); call compute_laplacian(v0, lin_v0); call compute_laplacian(w0, lin_w0)
    call compute_laplacian(u1, lin_u1); call compute_laplacian(v1, lin_v1); call compute_laplacian(w1, lin_w1)
    call compute_laplacian(u2, lin_u2); call compute_laplacian(v2, lin_v2); call compute_laplacian(w2, lin_w2)

    ! 5. Prediction AB2 (3 modes)
    if (first_step_ab) then
       rhs_u0 = sigma*u0 + 0.5_dp*nu*lin_u0 - nl_u0 - dxp0
       rhs_v0 = sigma*v0 + 0.5_dp*nu*lin_v0 - nl_v0 - dyp0
       rhs_w0 = sigma*w0 + 0.5_dp*nu*lin_w0 - nl_w0 - dzp0
       rhs_u1 = sigma*u1 + 0.5_dp*nu*lin_u1 - nl_u1 - dxp1
       rhs_v1 = sigma*v1 + 0.5_dp*nu*lin_v1 - nl_v1 - dyp1
       rhs_w1 = sigma*w1 + 0.5_dp*nu*lin_w1 - nl_w1 - dzp1
       rhs_u2 = sigma*u2 + 0.5_dp*nu*lin_u2 - nl_u2 - dxp2
       rhs_v2 = sigma*v2 + 0.5_dp*nu*lin_v2 - nl_v2 - dyp2
       rhs_w2 = sigma*w2 + 0.5_dp*nu*lin_w2 - nl_w2 - dzp2
    else
       rhs_u0 = sigma*u0 + 0.5_dp*nu*lin_u0 - (1.5_dp*nl_u0 - 0.5_dp*nl_u0_old) - dxp0
       rhs_v0 = sigma*v0 + 0.5_dp*nu*lin_v0 - (1.5_dp*nl_v0 - 0.5_dp*nl_v0_old) - dyp0
       rhs_w0 = sigma*w0 + 0.5_dp*nu*lin_w0 - (1.5_dp*nl_w0 - 0.5_dp*nl_w0_old) - dzp0
       rhs_u1 = sigma*u1 + 0.5_dp*nu*lin_u1 - (1.5_dp*nl_u1 - 0.5_dp*nl_u1_old) - dxp1
       rhs_v1 = sigma*v1 + 0.5_dp*nu*lin_v1 - (1.5_dp*nl_v1 - 0.5_dp*nl_v1_old) - dyp1
       rhs_w1 = sigma*w1 + 0.5_dp*nu*lin_w1 - (1.5_dp*nl_w1 - 0.5_dp*nl_w1_old) - dzp1
       rhs_u2 = sigma*u2 + 0.5_dp*nu*lin_u2 - (1.5_dp*nl_u2 - 0.5_dp*nl_u2_old) - dxp2
       rhs_v2 = sigma*v2 + 0.5_dp*nu*lin_v2 - (1.5_dp*nl_v2 - 0.5_dp*nl_v2_old) - dyp2
       rhs_w2 = sigma*w2 + 0.5_dp*nu*lin_w2 - (1.5_dp*nl_w2 - 0.5_dp*nl_w2_old) - dzp2
    end if

    ! Sauvegarde des vitesses (Coriolis explicite)
    u0_save = u0 ; v0_save = v0
    u1_save = u1 ; v1_save = v1
    u2_save = u2 ; v2_save = v2

    ! 6. Coriolis et viscosite
    select case (scheme_type)
    case (SCHEME_EXPLICIT_CLASSICAL)
       if (first_step_ab) then
          rhs_u0 = rhs_u0 - omega*v0 ; rhs_v0 = rhs_v0 + omega*u0
          rhs_u1 = rhs_u1 - omega*v1 ; rhs_v1 = rhs_v1 + omega*u1
          rhs_u2 = rhs_u2 - omega*v2 ; rhs_v2 = rhs_v2 + omega*u2
       else
          rhs_u0 = rhs_u0 - omega*(1.5_dp*v0 - 0.5_dp*v0_old)
          rhs_v0 = rhs_v0 + omega*(1.5_dp*u0 - 0.5_dp*u0_old)
          rhs_u1 = rhs_u1 - omega*(1.5_dp*v1 - 0.5_dp*v1_old)
          rhs_v1 = rhs_v1 + omega*(1.5_dp*u1 - 0.5_dp*u1_old)
          rhs_u2 = rhs_u2 - omega*(1.5_dp*v2 - 0.5_dp*v2_old)
          rhs_v2 = rhs_v2 + omega*(1.5_dp*u2 - 0.5_dp*u2_old)
       end if
       call eqn_u0%SOLVE(u0, rhs_u0, sigma, ph) ; call eqn_v0%SOLVE(v0, rhs_v0, sigma, ph)
       call eqn_u1%SOLVE(u1, rhs_u1, sigma, ph) ; call eqn_v1%SOLVE(v1, rhs_v1, sigma, ph)
       call eqn_u2%SOLVE(u2, rhs_u2, sigma, ph) ; call eqn_v2%SOLVE(v2, rhs_v2, sigma, ph)
    case (SCHEME_IMPLICIT_CLASSICAL, SCHEME_IMPLICIT_OLSHANSKII)
       rhs_u0 = rhs_u0 - omega_cn*v0 ; rhs_v0 = rhs_v0 + omega_cn*u0
       call eqn_uv0%SOLVE(u0, v0, rhs_u0, rhs_v0, sigma, -omega_cn, ph)
       rhs_u1 = rhs_u1 - omega_cn*v1 ; rhs_v1 = rhs_v1 + omega_cn*u1
       call eqn_uv1%SOLVE(u1, v1, rhs_u1, rhs_v1, sigma, -omega_cn, ph)
       rhs_u2 = rhs_u2 - omega_cn*v2 ; rhs_v2 = rhs_v2 + omega_cn*u2
       call eqn_uv2%SOLVE(u2, v2, rhs_u2, rhs_v2, sigma, -omega_cn, ph)
    end select

    call eqn_w0%SOLVE(w0, rhs_w0, sigma, ph)
    call eqn_w1%SOLVE(w1, rhs_w1, sigma, ph)
    call eqn_w2%SOLVE(w2, rhs_w2, sigma, ph)

    ! Sauvegarde des predictions
    u0_tilde = u0 ; v0_tilde = v0 ; w0_tilde = w0
    u1_tilde = u1 ; v1_tilde = v1 ; w1_tilde = w1
    u2_tilde = u2 ; v2_tilde = v2 ; w2_tilde = w2

    !===============================================================
    ! PROJECTION DE LA PRESSION (3 modes)
    !===============================================================

    call DIV(opx, opy, opz, u0, v0, w0, rhs_phi0, w01, w02, w03)
    isi = get_is() ; iei = get_ie()
    div_before = maxval(abs(rhs_phi0(isi(1):iei(1), isi(2):iei(2), isi(3):iei(3))))
    call MPI_ALLREDUCE(MPI_IN_PLACE, div_before, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

    call DIV(opx, opy, opz, u1, v1, w1, rhs_phi1, w01, w02, w03)
    call DIV(opx, opy, opz, u2, v2, w2, rhs_phi2, w01, w02, w03)

    select case (scheme_type)
    case (SCHEME_EXPLICIT_CLASSICAL, SCHEME_IMPLICIT_CLASSICAL)
       nu_poisson = [1.0_dp, 1.0_dp, 1.0_dp]
       isb = get_is_b([0, 0, 0]) ; ieb = get_ie_b([0, 0, 0])

       rhs_phi0 = rhs_phi0 / dt
       call eqn_p0%SOLVE_POISSON(phi0, rhs_phi0, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi0, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u0(i,j,k) = u0(i,j,k) - dt*w01(i,j,k)
          v0(i,j,k) = v0(i,j,k) - dt*w02(i,j,k)
          w0(i,j,k) = w0(i,j,k) - dt*w03(i,j,k)
       end forall

       rhs_phi1 = rhs_phi1 / dt
       call eqn_p1%SOLVE_POISSON(phi1, rhs_phi1, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi1, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u1(i,j,k) = u1(i,j,k) - dt*w01(i,j,k)
          v1(i,j,k) = v1(i,j,k) - dt*w02(i,j,k)
          w1(i,j,k) = w1(i,j,k) - dt*w03(i,j,k)
       end forall

       rhs_phi2 = rhs_phi2 / dt
       call eqn_p2%SOLVE_POISSON(phi2, rhs_phi2, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi2, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u2(i,j,k) = u2(i,j,k) - dt*w01(i,j,k)
          v2(i,j,k) = v2(i,j,k) - dt*w02(i,j,k)
          w2(i,j,k) = w2(i,j,k) - dt*w03(i,j,k)
       end forall

    case (SCHEME_IMPLICIT_OLSHANSKII)
       coef_h = sigma / (sigma**2 + omega_cn**2)
       coef_z = 1.0_dp / sigma
       nu_poisson = [coef_h, coef_h, coef_z]
       isb = get_is_b([0, 0, 0]) ; ieb = get_ie_b([0, 0, 0])

       call eqn_p0%SOLVE_POISSON(phi0, rhs_phi0, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi0, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u0(i,j,k) = u0(i,j,k) - (sigma*w01(i,j,k) - omega_cn*w02(i,j,k))/(sigma**2 + omega_cn**2)
          v0(i,j,k) = v0(i,j,k) - (sigma*w02(i,j,k) + omega_cn*w01(i,j,k))/(sigma**2 + omega_cn**2)
          w0(i,j,k) = w0(i,j,k) - w03(i,j,k)/sigma
       end forall

       call eqn_p1%SOLVE_POISSON(phi1, rhs_phi1, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi1, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u1(i,j,k) = u1(i,j,k) - (sigma*w01(i,j,k) - omega_cn*w02(i,j,k))/(sigma**2 + omega_cn**2)
          v1(i,j,k) = v1(i,j,k) - (sigma*w02(i,j,k) + omega_cn*w01(i,j,k))/(sigma**2 + omega_cn**2)
          w1(i,j,k) = w1(i,j,k) - w03(i,j,k)/sigma
       end forall

       call eqn_p2%SOLVE_POISSON(phi2, rhs_phi2, nu_poisson, ph)
       call GRAD(opx, opy, opz, phi2, w01, w02, w03)
       forall(i = isb(1):ieb(1), j = isb(2):ieb(2), k = isb(3):ieb(3))
          u2(i,j,k) = u2(i,j,k) - (sigma*w01(i,j,k) - omega_cn*w02(i,j,k))/(sigma**2 + omega_cn**2)
          v2(i,j,k) = v2(i,j,k) - (sigma*w02(i,j,k) + omega_cn*w01(i,j,k))/(sigma**2 + omega_cn**2)
          w2(i,j,k) = w2(i,j,k) - w03(i,j,k)/sigma
       end forall
    end select

    select case (scheme_type)
    case (SCHEME_EXPLICIT_CLASSICAL)
       call eqn_u0%APPLY_BCS(u0, ph) ; call eqn_v0%APPLY_BCS(v0, ph)
       call eqn_u1%APPLY_BCS(u1, ph) ; call eqn_v1%APPLY_BCS(v1, ph)
       call eqn_u2%APPLY_BCS(u2, ph) ; call eqn_v2%APPLY_BCS(v2, ph)
    case (SCHEME_IMPLICIT_CLASSICAL, SCHEME_IMPLICIT_OLSHANSKII)
       call eqn_uv0%apply_bcs(u0, v0, ph)
       call eqn_uv1%apply_bcs(u1, v1, ph)
       call eqn_uv2%apply_bcs(u2, v2, ph)
    end select

    ! Pressions totales
    pres0 = pres0 + phi0
    pres1 = pres1 + phi1
    pres2 = pres2 + phi2

    call compute_divergence_error()
    div_after = div_max

    ! Historique NL pour AB2
    nl_u0_old = nl_u0 ; nl_v0_old = nl_v0 ; nl_w0_old = nl_w0 ; nl_t0_old = nl_t0
    nl_u1_old = nl_u1 ; nl_v1_old = nl_v1 ; nl_w1_old = nl_w1 ; nl_t1_old = nl_t1
    nl_u2_old = nl_u2 ; nl_v2_old = nl_v2 ; nl_w2_old = nl_w2 ; nl_t2_old = nl_t2

    ! Historique vitesses pour AB2 (Coriolis explicite)
    u0_old = u0_save ; v0_old = v0_save
    u1_old = u1_save ; v1_old = v1_save
    u2_old = u2_save ; v2_old = v2_save

    first_step_ab = .false.
  end subroutine advance_abcn_step

  !===================================================================
  ! Diagnostics and outputs
  !===================================================================

  subroutine write_diagnostics_and_outputs()
    real(dp), parameter :: eps_growth = 1.0e-300_dp
    real(dp), save :: kinetic_w_old = -1.0_dp
    real(dp), save :: growth_rate_w = 0.0_dp

    call compute_divergence_error()
    call compute_vertical_kinetic_energy()
    call compute_geostrophic_error()

    call GetCFL(msh(1), msh(2), msh(3), u0, v0, w0, dt, cfl)

    if (kinetic_w_old > eps_growth .and. kinetic_w > eps_growth) then
       growth_rate_w = 0.5_dp * log(kinetic_w / kinetic_w_old) / dt
    else
       growth_rate_w = 0.0_dp
    end if
    kinetic_w_old = kinetic_w

    step_walltime = MPI_Wtime() - wall_start
    call MPI_ALLREDUCE(MPI_IN_PLACE, step_walltime, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (rank == 0) then
       print '(i9,13(1x,e15.8))', it_time, t_cur, dt, cfl, kinetic_w, std_w, growth_rate_w, &
                                  div_max, step_walltime, &
                                  geostrophic_error, ageostrophic_error, splitting_error, &
                                  sobol_alpha, sobol_paroi
       ! fort.20 : t, kinetic_w, std_w, growth, var_alpha, var_paroi, sobol_alpha, sobol_paroi
       write(20, '(8(1x,e15.8))') t_cur, kinetic_w, std_w, growth_rate_w, &
                                  var_alpha, var_paroi, sobol_alpha, sobol_paroi
    end if

    if (mod(it_time, dn_dump) == 0) then
       call dump_state('WRITE')
    end if

    if (mod(it_time, dn_snap) == 0) then
       ! Mode 0 (moyenne)
       call reserve_next_snapshot_filename(filename)
       call EXPORT_snapshot(trim(filename), t_cur, u0, v0, w0, pres0, temp0)
       ! Mode 1 (sensibilite alpha)
       call reserve_next_snapshot_filename(filename)
       call EXPORT_snapshot(trim(filename), t_cur, u1, v1, w1, pres1, temp1)
       ! Mode 2 (sensibilite paroi)
       call reserve_next_snapshot_filename(filename)
       call EXPORT_snapshot(trim(filename), t_cur, u2, v2, w2, pres2, temp2)
    end if
  end subroutine write_diagnostics_and_outputs

  subroutine compute_divergence_error()
    integer :: is_loc2(3), ie_loc2(3)
    call DIV(opx, opy, opz, u0, v0, w0, w09, w01, w02, w03)
    is_loc2 = get_is()
    ie_loc2 = get_ie()
    div_max = maxval(abs(w09(is_loc2(1):ie_loc2(1), is_loc2(2):ie_loc2(2), is_loc2(3):ie_loc2(3))))
    call MPI_ALLREDUCE(MPI_IN_PLACE, div_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  end subroutine compute_divergence_error

subroutine compute_vertical_kinetic_energy()
    real(dp) :: e00, e11, e22, e01, e02
    real(dp), parameter :: tiny_v = 1.0e-300_dp
    ! Energie verticale <w^2> et decomposition de la variance par parametre.
    ! w = w0 + w1 xi1 + w2 xi2  (xi ~ N(0,1) independants)
    !
    ! On calcule les moments spatiaux :
    !   e00 = <w0^2>   e11 = <w1^2>   e22 = <w2^2>
    !   e01 = <w0*w1>  e02 = <w0*w2>  (covariances mode0-mode1 et mode0-mode2)
    w01 = w0*w0 ; e00 = quad_mean(w01)
    w01 = w1*w1 ; e11 = quad_mean(w01)
    w01 = w2*w2 ; e22 = quad_mean(w01)
    w01 = w0*w1 ; e01 = quad_mean(w01)
    w01 = w0*w2 ; e02 = quad_mean(w01)

    ! Esperance de l'energie : E[<w^2>] = <w0^2> + <w1^2> + <w2^2>
    kinetic_w = e00 + e11 + e22

    ! Variance de l'ENERGIE <w^2> (formule a 3 modes, coherente avec le cas 2-modes) :
    !   Var = 4*e01^2 + 2*e11^2  (termes du mode 1, alpha)
    !       + 4*e02^2 + 2*e22^2  (termes du mode 2, paroi)
    ! L'ecart-type total est la racine de cette variance.
    std_w = sqrt(max(4.0_dp*e01**2 + 2.0_dp*e11**2 &
                   + 4.0_dp*e02**2 + 2.0_dp*e22**2, 0.0_dp))

    ! Decomposition de la variance par parametre (pour les indices de Sobol).
    ! On attribue a chaque parametre sa contribution propre a la variance de l'energie :
    !   var_alpha = 4*e01^2 + 2*e11^2   (du a alpha, xi1)
    !   var_paroi = 4*e02^2 + 2*e22^2   (du a la paroi, xi2)
    var_alpha = 4.0_dp*e01**2 + 2.0_dp*e11**2
    var_paroi = 4.0_dp*e02**2 + 2.0_dp*e22**2

    ! Indices de Sobol (premier ordre)
    sobol_alpha = var_alpha / max(var_alpha + var_paroi, tiny_v)
    sobol_paroi = var_paroi / max(var_alpha + var_paroi, tiny_v)
  end subroutine compute_vertical_kinetic_energy


  subroutine compute_geostrophic_error()
    real(dp) :: mean_num, mean_den
    real(dp) :: mean_ag, mean_u
    real(dp) :: mean_rec, mean_grad_phi
    real(dp) :: sigma_loc, omega_rec

    call GRAD(opx, opy, opz, pres0, w01, w02, w03)
    w04 = (w01 + omega*v0)**2 + (w02 - omega*u0)**2
    w05 = (omega*v0)**2 + (omega*u0)**2 + w01**2 + w02**2
    mean_num = quad_mean(w04)
    mean_den = quad_mean(w05)
    geostrophic_error = sqrt(mean_num)

    if (abs(omega) > 1.0e-300_dp) then
       w04 = (u0 - w02/omega)**2 + (v0 + w01/omega)**2
       w05 = u0**2 + v0**2
       mean_ag = quad_mean(w04)
       mean_u  = quad_mean(w05)
       ageostrophic_error = sqrt(mean_ag) / max(sqrt(mean_u), 1.0e-300_dp)
    else
       ageostrophic_error = 0.0_dp
    end if

    sigma_loc = 1.0_dp / dt
    if (scheme_type == SCHEME_EXPLICIT_CLASSICAL) then
       omega_rec = 0.0_dp
    else
       omega_rec = 0.5_dp * omega
    end if

    call GRAD(opx, opy, opz, phi0, w01, w02, w03)
    w04 = (sigma_loc*(u0 - u0_tilde) + omega_rec*(v0 - v0_tilde) + w01)**2 + &
          (sigma_loc*(v0 - v0_tilde) - omega_rec*(u0 - u0_tilde) + w02)**2 + &
          (sigma_loc*(w0 - w0_tilde)                             + w03)**2
    w05 = w01**2 + w02**2 + w03**2
    mean_rec      = quad_mean(w04)
    mean_grad_phi = quad_mean(w05)
    splitting_error = sqrt(mean_rec) / max(sqrt(mean_grad_phi), 1.0e-300_dp)
  end subroutine compute_geostrophic_error

  !===================================================================
  ! I/O helpers
  !===================================================================

  subroutine reserve_next_snapshot_filename(file_name)
    character(len=*), intent(out) :: file_name
    file_name = trim(base_snap) // 'grid.h5'
    call update_id_and_time(trim(file_name), snap_time, snap_id)
    write(snap_num, '(I6.6)') snap_id
    file_name = trim(base_snap) // snap_num // '.h5'
  end subroutine reserve_next_snapshot_filename

  subroutine dump_state(mode)
    character(len=*), intent(in) :: mode
    call DUMP_HDF5_BASIC(file_dump, mode, t_cur, dt, msh, u0, v0, w0, pres0, temp0, w09)
  end subroutine dump_state

  subroutine write_final_state_and_stop()
    call reserve_next_snapshot_filename(filename)
    call EXPORT_snapshot(trim(filename), t_cur, u0, v0, w0, pres0, temp0)
    call reserve_next_snapshot_filename(filename)
    call EXPORT_snapshot(trim(filename), t_cur, u1, v1, w1, pres1, temp1)
    call reserve_next_snapshot_filename(filename)
    call EXPORT_snapshot(trim(filename), t_cur, u2, v2, w2, pres2, temp2)
    call dump_state('WRITE')
    call mpi_finalize(ierr)
    stop
  end subroutine write_final_state_and_stop

  !===================================================================
  ! Local index helpers
  !===================================================================

  function get_is_b(halo) result(res)
    integer, intent(in) :: halo(3)
    integer :: res(3)
    res = ph%xst
    if (ph%xst(1) == 1) res(1) = res(1) + halo(1)
    if (ph%xst(2) == 1) res(2) = res(2) + halo(2)
    if (ph%xst(3) == 1) res(3) = res(3) + halo(3)
  end function get_is_b

  function get_ie_b(halo) result(res)
    integer, intent(in) :: halo(3)
    integer :: res(3)
    res = ph%xen
    if (ph%xen(1) == n(1) + 1) res(1) = res(1) - halo(1)
    if (ph%xen(2) == n(2) + 1) res(2) = res(2) - halo(2)
    if (ph%xen(3) == n(3) + 1) res(3) = res(3) - halo(3)
  end function get_ie_b

  function get_is() result(res)
    integer :: res(3)
    res = ph%xst
    if (ph%xst(1) == 1) res(1) = res(1) + 1
    if (ph%xst(2) == 1) res(2) = res(2) + 1
    if (ph%xst(3) == 1) res(3) = res(3) + 1
  end function get_is

  function get_ie() result(res)
    integer :: res(3)
    res = ph%xen
    if (ph%xen(1) == n(1) + 1) res(1) = res(1) - 1
    if (ph%xen(2) == n(2) + 1) res(2) = res(2) - 1
    if (ph%xen(3) == n(3) + 1) res(3) = res(3) - 1
  end function get_ie

  !===================================================================
  ! CFL routines
  !===================================================================

  subroutine GetCFL(msh_x, msh_y, msh_z, ux, uy, uz, dt_in, cfl_out)
    type(t_mesh_base), intent(in) :: msh_x, msh_y, msh_z
    real(dp), allocatable, intent(in) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(dp), intent(in)  :: dt_in
    real(dp), intent(out) :: cfl_out
    real(dp) :: local_cfl, speed_over_dx, dx, dy, dz
    integer :: ii, jj, kk
    integer :: is_loc2(3), ie_loc2(3)
    is_loc2 = get_is()
    ie_loc2 = get_ie()
    local_cfl = 0.0_dp
    do kk = is_loc2(3), ie_loc2(3)
       do jj = is_loc2(2), ie_loc2(2)
          do ii = is_loc2(1), ie_loc2(1)
             dx = msh_x%x(ii+1) - msh_x%x(ii)
             dy = msh_y%x(jj+1) - msh_y%x(jj)
             dz = msh_z%x(kk+1) - msh_z%x(kk)
             speed_over_dx = abs(ux(ii,jj,kk))/dx + abs(uy(ii,jj,kk))/dy + abs(uz(ii,jj,kk))/dz
             local_cfl = max(local_cfl, speed_over_dx * dt_in)
          end do
       end do
    end do
    cfl_out = local_cfl
    call MPI_ALLREDUCE(MPI_IN_PLACE, cfl_out, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  end subroutine GetCFL

  !===================================================================
  ! User-defined boundary/value functions
  !===================================================================

  real(dp) pure function udf_null(time_arg, x_arg, y_arg, z_arg)
    real(dp), intent(in) :: time_arg, x_arg, y_arg, z_arg
    udf_null = 0.0_dp
  end function udf_null

  ! CL incertaine sur la paroi chaude (bas) pour le MODE 2.
  ! La valeur retournee est directement la valeur Dirichlet imposee a la paroi.
  ! Utilise la CONSTANTE BC_HOT_AMP_PARAM (obligatoire car fonction pure).
  real(dp) pure function udf_bc_hot(time_arg, x_arg, y_arg, z_arg)
    real(dp), intent(in) :: time_arg, x_arg, y_arg, z_arg
    udf_bc_hot = BC_HOT_AMP_PARAM
  end function udf_bc_hot

  real(dp) pure function udf_unit(time_arg, x_arg, y_arg, z_arg)
    real(dp), intent(in) :: time_arg, x_arg, y_arg, z_arg
    udf_unit = 1.0_dp
  end function udf_unit

  real(dp) pure function udf_linear(time_arg, x_arg, y_arg, z_arg)
    real(dp), intent(in) :: time_arg, x_arg, y_arg, z_arg
    udf_linear = 1.0_dp - x_arg
  end function udf_linear

end program tcheby_1d_pce