!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine beam_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &                  ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: I,J

    real (prec)  ::  xg(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  u(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  ul(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E,G,A,k_prime(3),Izz,Iyy,Jp,kappa,ly,hz,qv,qw,leng             ! Material properties
    real (prec)  ::  T(3,3),Q(12,12)
    real (prec)  ::  D1,D2,D3,D4,D5,D6,D7,D8,temp
    real (prec)  ::  phi_y, phi_z, phil_y,phil_z

    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:


    fail = .false.

    ! x = [u1 v1 w1 ; θx1 θy1 θz1 ; u2 v2 w2 ; θx2 θy2 θz2]^T
    
    xg = reshape(element_coords+dof_total+dof_increment,(/3,length_coord_array/3/))
    u = reshape(dof_increment,(/3,length_coord_array/3/))
write(*,*) u
    E = element_properties(1)
    G = element_properties(2)
    A = element_properties(3)
    k_prime(1:3) = element_properties(4:6)
    Izz = element_properties(7)
    Iyy = element_properties(8)
    Jp = element_properties(9)
    kappa = element_properties(10)
    ly = element_properties(11)
    hz = element_properties(12)
    qv = element_properties(13)
    qw = element_properties(14)
    leng = dsqrt(((xg(1,1)-xg(1,3))**2 + (xg(2,1)-xg(2,3))**2 + (xg(3,1)-xg(3,3))**2)**2)

    ! coordinate transformation
    T(1,1) = (xg(1,3) - xg(1,1))/leng
    T(1,2) = (xg(2,3) - xg(2,1))/leng
    T(1,3) = (xg(3,3) - xg(3,1))/leng

    temp = dsqrt((T(1,2)*k_prime(3) - T(1,3)*k_prime(2))**2+(T(1,3)*k_prime(1) - &
    T(1,1)*k_prime(3))**2+(T(1,1)*k_prime(2) - T(1,2)*k_prime(1))**2)
    T(2,1) = -(T(1,2)*k_prime(3) - T(1,3)*k_prime(2))/temp
    T(2,2) = -(T(1,3)*k_prime(1) - T(1,1)*k_prime(3))/temp
    T(2,3) = -(T(1,1)*k_prime(2) - T(1,2)*k_prime(1))/temp

    temp = dsqrt((T(1,2)*T(2,3) - T(1,3)*T(2,2))**2+(T(1,3)*T(2,1) - &
    T(1,1)*T(2,3))**2+(T(1,1)*T(2,2) - T(1,2)*T(2,1))**2)
    T(3,1) = (T(1,2)*T(2,3) - T(1,3)*T(2,2))/temp
    T(3,2) = (T(1,3)*T(2,1) - T(1,1)*T(2,3))/temp
    T(3,3) = (T(1,1)*T(2,2) - T(1,2)*T(2,1))/temp

!    write(6,"(3D12.5,/,3D12.5,/,3D12.5)")  ((T(i,j),j=1,3),i=1,3)

    Q = 0.d0
    Q(1:3,1:3) = T
    Q(4:6,4:6) = T
    Q(7:9,7:9) = T
    Q(10:12,10:12) = T

    phi_y = 12*E*Iyy/(kappa*G*A*leng*leng)
    phi_z = 12*E*Izz/(kappa*G*A*leng*leng)
    phil_y = 1/(1 + phi_y)
    phil_z = 1/(1 + phi_z)

    element_residual = 0.d0
    element_stiffness = 0.d0

    D1 = E*A/leng
    D2 = 12.d0*phil_z*E*Izz/leng**3
    D3 = 12.d0*phil_y*E*Iyy/leng**3
    D4 = G*Jp/leng
    D5 = (4+phi_y)*phil_y*E*Iyy/leng
    D6 = (4+phi_z)*phil_z*E*Izz/leng
    D7 =  (2-phi_y)*phil_y*E*Iyy/leng
    D8 =  (2-phi_z)*phil_z*E*Izz/leng

    element_stiffness(1,1) = D1
    element_stiffness(7,1) = -D1
    element_stiffness(2,2) = D2
    element_stiffness(6,2) = D2/2.d0*leng
    element_stiffness(8,2) = -D2
    element_stiffness(12,2) = D2/2.d0*leng
    element_stiffness(3,3) = D3
    element_stiffness(5,3) = -D3/2.d0*leng
    element_stiffness(9,3) = -D3
    element_stiffness(11,3) = -D3/2.d0*leng
    element_stiffness(4,4) = D4
    element_stiffness(10,4) = -D4
    element_stiffness(5,5) = D5
    element_stiffness(9,5) = D3/2.d0*leng
    element_stiffness(11,5) = D7
    element_stiffness(6,6) = D6
    element_stiffness(8,6) = -D2/2.d0*leng
    element_stiffness(12,6) = D8
    element_stiffness(7,7) = D1
    element_stiffness(8,8) = D2
    element_stiffness(12,8) = -D2/2.d0*leng
    element_stiffness(9,9) = D3
    element_stiffness(11,9) = D3/2.d0*leng
    element_stiffness(10,10) = D4
    element_stiffness(11,11) = D5
    element_stiffness(12,12) = D6
    do i = 1,12
        do j = i,12
            if(j > i) then
                element_stiffness(i,j) = element_stiffness(j,i)
            endif
        enddo
    enddo
  !    R= l/12 [0 6qv0 6qw0 0 lqw0 lqv0 0 −6qv0 6qw0 0 −lqw0 lqv0]T
    element_residual(2) = 6.d0*qv
    element_residual(3) = 6.d0*qw
    element_residual(5) = -leng*qw
    element_residual(6) = leng*qv
    element_residual(8) = 6.d0*qv
    element_residual(9) = 6.d0*qw
    element_residual(11) = leng*qw
    element_residual(12) = -leng*qv
    element_residual = element_residual/12.d0*leng

!    write(6,'(A16,12D12.5)') 'element_residual',element_residual

    element_residual = matmul(transpose(Q),element_residual)
    element_stiffness = matmul(matmul(transpose(Q),element_stiffness),Q)

        ! x = [u1 v1 w1 ; θx1 θy1 θz1 ; u2 v2 w2 ; θx2 θy2 θz2]^T
    ul = reshape(matmul(Q,(dof_increment)),(/3,length_coord_array/3/))
    ! stress
    updated_state_variables = 0.d0
    updated_state_variables(1) = abs(E*(ul(1,1)-ul(1,3))/leng)
    updated_state_variables(2) = abs(G*(ul(1,2)-ul(1,4))/leng)*max(ly,hZ)
    updated_state_variables(3) = abs(E*ly*(ul(3,2)-ul(3,4))/leng)
    updated_state_variables(5) = abs(E*hz*(ul(2,2)-ul(2,4))/leng)
    updated_state_variables(4) = abs(G*phil_z*phi_z*(2*ul(2,1)+ul(3,2)*leng-2*ul(2,3)+ul(3,4)*leng)/2/leng)
    updated_state_variables(6) = abs(G*phil_y*phi_y*(2*ul(3,1)+ul(2,2)*leng-2*ul(3,3)+ul(2,4)*leng)/2/leng)

    !write(*,*)  updated_state_variables

    return
end subroutine beam_3dbasic



!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_beam_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k,I,J

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    real (prec)  ::  B_diff(20,3)                      !
    real (prec)  ::  el_vol                            ! element volume

  
    return
end subroutine fieldvars_beam_3dbasic

