!     Subroutines for basic 3D  hyperelastic elements



!==========================SUBROUTINE el_hyperelast_3d ==============================
subroutine el_hyperelast_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_Eulerian_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only:  dNbardy => vol_avg_shape_function_derivatives_y_3D
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
    integer      :: n_points,kint,I,J

    real (prec)  ::  B_strain(3,3),inv_B(3,3),B_kk,Bvec(6),inv_Bvec(6),J_b           ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  K_stress(6),stress(6),m_stress(3,3),tau(6),tau_kk               ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                                                          ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array),B_bar(6,length_dof_array)                 ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant                                         ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                                       ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  K1, mu1                                                         ! Material properties
    real (prec)  ::  el_vol,B_diff(20,3)                                             ! element volume
    real (prec)  ::  F(3,3),inv_F(3,3),J_f,eta,F_bar(3,3)                            ! deformation gradient
    real (prec)  ::  G(6,9), B_star(9,3*n_nodes)                                     !
    real (prec)  ::  u(3,length_dof_array/3)                                         ! displacement and displacement_increment
    real (prec)  ::  delta1(3,3),delta2(6,6),identity(6)                             ! Constant matrix
    real (prec)  ::  IinvB(6,6),BinvB(6,6),II(6,6)                                   !
    real (prec)  ::  P(3*n_nodes,3*n_nodes),Q(3*n_nodes,3*n_nodes),Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  dNdy_akbi(3*n_nodes,3*n_nodes),dNdy_aibk(3*n_nodes,3*n_nodes),dNbardy_aibk(3*n_nodes,3*n_nodes)
    real (prec)  ::  dNdyvec(3*n_nodes),dNbardyvec(3*n_nodes),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Svec(3*n_nodes),Smat(3*n_nodes,3*n_nodes),S(3,n_nodes)


    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
    u = reshape(dof_total+dof_increment,(/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    ! constant matrix

    delta1 = 0.d0
    delta1(1,1) = 1.d0
    delta1(2,2) = 1.d0
    delta1(3,3) = 1.d0
    delta2 = 0.d0
    delta2(1:3,1:3) = delta1
    delta2(4,4) = 1.d0/2.d0
    delta2(5,5) = 1.d0/2.d0
    delta2(6,6) = 1.d0/2.d0
    identity = 0.d0
    identity(1:3) = 1.d0


  ! Initialize any variables that are incremented in the loops below

   element_residual = 0.d0
   element_stiffness = 0.d0
   el_vol = 0.d0
   dNbardy = 0.d0
   P = 0.d0
   eta = 0.d0
   !    Initialize integration points and weights
    call initialize_integration_points(n_points, n_nodes, xi, w)
   !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

   !  Compute the deformation gradients by differentiating displacements
       F = delta1 + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
   !    Compute inverse of the deformation gradient and J
       call invert_small(F,inv_F,J_f)
   !    Convert the shape function derivatives to derivatives with respect to deformed coordinates
       dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F)
   !    Add contribution to dNbardy and J_bar from current integration point
       dNbardy(1:n_nodes,1:3) = dNbardy(1:n_nodes,1:3) + dNdy(1:n_nodes,1:3)*J_f*w(kint)*determinant
    !   Add contribution to P_1 from current integration point
       dNdyvec(1:3*n_nodes) = reshape(transpose(dNdy(1:n_nodes,1:3)),(/3*n_nodes/))
       dNdy_aibk =spread(dNdyvec,dim=2,ncopies=3*n_nodes)*spread(dNdyvec,dim=1,ncopies=3*n_nodes)
       do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
       end do
            dNdy_akbi = Pmat*transpose(Pmat)
            P = P + (dNdy_aibk-dNdy_akbi)*J_f*w(kint)*determinant
  !   Add contribution to element volume from current integration point
            eta = eta + J_f*w(kint)*determinant
            el_vol = el_vol + w(kint)*determinant
    enddo
            eta = eta/el_vol
            dNbardy = dNbardy/el_vol/eta
    !   expression for P
       dNbardyvec = reshape(transpose(dNbardy(1:n_nodes,1:3)),(/3*n_nodes/))
       dNbardy_aibk =spread(dNbardyvec,dim=2,ncopies=3*n_nodes)*spread(dNbardyvec,dim=1,ncopies=3*n_nodes)
       P = P/el_vol/eta - dNbardy_aibk
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
   !  Compute the deformation gradients by differentiating displacements
       F = delta1 + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
   !    Compute inverse of the deformation gradient and J
       call invert_small(F,inv_F,J_f)
   !    Convert the shape function derivatives to derivatives with respect to deformed coordinates
       dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F)
        F_bar = (eta/J_f)**(1.d0/3.d0) * F
   !    Compute the Kirchhoff stress and material tangent stiffness
    mu1 = element_properties(1)
    K1 = element_properties(2)
   !    Compute Left Cauchy Green deformation tensor
    B_strain = matmul(F_bar,transpose(F_bar))
    B_kk = B_strain(1,1)+B_strain(2,2)+B_strain(3,3)
   !   Compute the Kirchhoff stress
    m_stress = mu1*(B_strain- delta1*B_kk/3.d0)/(eta**(5.d0/3.d0))+K1*(eta-1.d0)*delta1
    stress(1) = m_stress(1,1)
    stress(2) = m_stress(2,2)
    stress(3) = m_stress(3,3)
    stress(4) = m_stress(1,2)
    stress(5) = m_stress(1,3)
    stress(6) = m_stress(2,3)
    tau = stress*eta
    tau_kk = sum(tau(1:3))
   ! Compute the Tangent Stiffness D
    call invert_small(B_strain,inv_B,J_b)
    Bvec(1) = B_strain(1,1)
    Bvec(2) = B_strain(2,2)
    Bvec(3) = B_strain(3,3)
    Bvec(4) = B_strain(1,2)
    Bvec(5) = B_strain(1,3)
    Bvec(6) = B_strain(2,3)
    inv_Bvec(1) = inv_B(1,1)
    inv_Bvec(2) = inv_B(2,2)
    inv_Bvec(3) = inv_B(3,3)
    inv_Bvec(4) = inv_B(1,2)
    inv_Bvec(5) = inv_B(1,3)
    inv_Bvec(6) = inv_B(2,3)
    IinvB = spread(identity,dim=2,ncopies=6)*spread(inv_Bvec,dim=1,ncopies=6)
    BinvB = spread(Bvec,dim=2,ncopies=6)    *spread(inv_Bvec,dim=1,ncopies=6)
    II    = spread(identity,dim=2,ncopies=6)*spread(identity,dim=1,ncopies=6)
    D = mu1/(eta**(2.d0/3.d0))*delta2 + mu1/(eta**(2.d0/3.d0))/3*(B_kk/3.d0*IinvB-II-BinvB)+K1*eta*(eta-1.d0/2.d0)*IinvB
   !    Assemble B_bar
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        B_bar = B
        B_diff = (dNbardy - dNdy)/3.d0
            do I = 1,3
            B_bar(I,1:3*n_nodes-2:3)= B_bar(I,1:3*n_nodes-2:3) + B_diff(1:n_nodes,1)
            B_bar(I,2:3*n_nodes-1:3)= B_bar(I,2:3*n_nodes-1:3) + B_diff(1:n_nodes,2)
            B_bar(I,3:3*n_nodes :3) = B_bar(I,3:3*n_nodes :3) + B_diff(1:n_nodes,3)
            enddo
     !    Assemble matrix G
    G = 0.d0
    G(1,1) = Bvec(1)
    G(2,2) = Bvec(2)
    G(3,3) = Bvec(3)
    G(4,1) = Bvec(4)
    G(4,2) = Bvec(4)
    G(5,1) = Bvec(5)
    G(5,3) = Bvec(5)
    G(6,2) = Bvec(6)
    G(6,3) = Bvec(6)
    G(1,4) = Bvec(4)
    G(2,5) = Bvec(4)
    G(4,4) = Bvec(2)
    G(4,5) = Bvec(1)
    G(5,4) = Bvec(6)
    G(6,5) = Bvec(5)
    G(1,6) = Bvec(5)
    G(3,7) = Bvec(5)
    G(4,6) = Bvec(6)
    G(5,6) = Bvec(3)
    G(5,7) = Bvec(1)
    G(6,7) = Bvec(4)
    G(2,8) = Bvec(6)
    G(3,9) = Bvec(6)
    G(4,8) = Bvec(5)
    G(5,9) = Bvec(4)
    G(6,8) = Bvec(3)
    G(6,9) = Bvec(2)
    G = 2*G
     !    Assemble matrix B_star
        B_star = 0.d0
        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        B_diff = (dNbardy - dNdy)/3.d0
            do I = 1,3
            B_star(I,1:3*n_nodes-2:3)= B_star(I,1:3*n_nodes-2:3) + B_diff(1:n_nodes,1)
            B_star(I,2:3*n_nodes-1:3)= B_star(I,2:3*n_nodes-1:3) + B_diff(1:n_nodes,2)
            B_star(I,3:3*n_nodes :3) = B_star(I,3:3*n_nodes :3) + B_diff(1:n_nodes,3)
            enddo
      ! Compute Sigma and Q
        S = reshape(matmul(transpose(B),tau),(/3,length_dof_array/3/))
      do i = 1,n_nodes
        Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
        Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
        Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
        Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
      end do
        Sigma = Pmat*transpose(Smat)
       do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
       end do
            Q = Pmat*transpose(Pmat)
      ! Add the contributions from the current integration point to the integrals
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - &
        matmul(transpose(B_bar(1:6,1:3*n_nodes)),tau)*w(kint)*determinant
         !    write(6,*) tau(1:6)
        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B_bar(1:6,1:3*n_nodes)),matmul(matmul(D,G),B_star))*w(kint)*determinant&
            - Sigma*w(kint)*determinant + tau_kk/3.d0*(P+Q)*w(kint)*determinant
    end do
    return
end subroutine el_hyperelast_3d


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_hyperelast_3d_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
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
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

    if (n_properties == 2) then

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu)
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44

    stress = matmul(D,strain+dstrain)

    else if (n_properties == 4) then

    call  hypoelastic_tangent_3D(strain+dstrain,n_properties,element_properties,stress,D)

    end if


        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_hyperelast_3d_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperelast_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_Eulerian_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only:  dNbardy => vol_avg_shape_function_derivatives_y_3D
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

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    logical      :: strcmp

    ! Local Variables
    integer      :: n_points,kint,I,J,k

    real (prec)  ::  B_strain(3,3),inv_B(3,3),B_kk,Bvec(6),inv_Bvec(6),J_b           ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  K_stress(6),stress(6),m_stress(3,3),tau(6),tau_kk               ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                                                          ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array),B_bar(6,length_dof_array)                 ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant                                         ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)                                       ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  K1, mu1                                                         ! Material properties
    real (prec)  ::  el_vol,B_diff(20,3)                                             ! element volume
    real (prec)  ::  F(3,3),inv_F(3,3),J_f,eta,F_bar(3,3)                            ! deformation gradient
    real (prec)  ::  G(6,9), B_star(9,3*n_nodes)                                     !
    real (prec)  ::  u(3,length_dof_array/3)                                         ! displacement and displacement_increment
    real (prec)  ::  delta1(3,3),delta2(6,6),identity(6)                             ! Constant matrix
    real (prec)  ::  IinvB(6,6),BinvB(6,6),II(6,6)                                   !
    real (prec)  ::  P(3*n_nodes,3*n_nodes),Q(3*n_nodes,3*n_nodes),Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  dNdy_akbi(3*n_nodes,3*n_nodes),dNdy_aibk(3*n_nodes,3*n_nodes),dNbardy_aibk(3*n_nodes,3*n_nodes)
    real (prec)  ::  dNdyvec(3*n_nodes),dNbardyvec(3*n_nodes),Pvec(3*n_nodes),Pmat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Svec(3*n_nodes),Smat(3*n_nodes,3*n_nodes),S(3,n_nodes)
    real (prec)  ::  smises,sdev(6)                                                           ! Mises stress


    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio


    x = reshape(element_coords,(/3,length_coord_array/3/))
    u = reshape(dof_total+dof_increment,(/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    ! constant matrix

    delta1 = 0.d0
    delta1(1,1) = 1.d0
    delta1(2,2) = 1.d0
    delta1(3,3) = 1.d0
    delta2 = 0.d0
    delta2(1:3,1:3) = delta1
    delta2(4,4) = 1.d0/2.d0
    delta2(5,5) = 1.d0/2.d0
    delta2(6,6) = 1.d0/2.d0
    identity = 0.d0
    identity(1:3) = 1.d0


  ! Initialize any variables that are incremented in the loops below

   el_vol = 0.d0
   dNbardy = 0.d0
   P = 0.d0
   eta = 0.d0
   !    Initialize integration points and weights
    call initialize_integration_points(n_points, n_nodes, xi, w)
    nodal_fieldvariables = 0.d0

   !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

   !  Compute the deformation gradients by differentiating displacements
       F = delta1 + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
   !    Compute inverse of the deformation gradient and J
       call invert_small(F,inv_F,J_f)
   !    Convert the shape function derivatives to derivatives with respect to deformed coordinates
       dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F)
   !    Add contribution to dNbardy and J_bar from current integration point
       dNbardy(1:n_nodes,1:3) = dNbardy(1:n_nodes,1:3) + dNdy(1:n_nodes,1:3)*J_f*w(kint)*determinant
    !   Add contribution to P_1 from current integration point
       dNdyvec(1:3*n_nodes) = reshape(transpose(dNdy(1:n_nodes,1:3)),(/3*n_nodes/))
       dNdy_aibk =spread(dNdyvec,dim=2,ncopies=3*n_nodes)*spread(dNdyvec,dim=1,ncopies=3*n_nodes)
       do i = 1,n_nodes
            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
       end do
            dNdy_akbi = Pmat*transpose(Pmat)
            P = P + (dNdy_aibk-dNdy_akbi)*J_f*w(kint)*determinant
  !   Add contribution to element volume from current integration point
            eta = eta + J_f*w(kint)*determinant
            el_vol = el_vol + w(kint)*determinant
    enddo
            eta = eta/el_vol
            dNbardy = dNbardy/el_vol/eta
    !   expression for P
       dNbardyvec = reshape(transpose(dNbardy(1:n_nodes,1:3)),(/3*n_nodes/))
       dNbardy_aibk =spread(dNbardyvec,dim=2,ncopies=3*n_nodes)*spread(dNbardyvec,dim=1,ncopies=3*n_nodes)
            P = P/el_vol/eta - dNbardy_aibk


    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
   !  Compute the deformation gradients by differentiating displacements
       F = delta1 + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
   !    Compute inverse of the deformation gradient and J
       call invert_small(F,inv_F,J_f)
   !    Convert the shape function derivatives to derivatives with respect to deformed coordinates
       dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_F)
        F_bar = (eta/J_f)**(1.d0/3.d0) * F
   !    Compute the Kirchhoff stress and material tangent stiffness
    mu1 = element_properties(1)
    K1 = element_properties(2)
   !    Compute Left Cauchy Green deformation tensor
    B_strain = matmul(F_bar,transpose(F_bar))
    B_kk = B_strain(1,1)+B_strain(2,2)+B_strain(3,3)
   !   Compute the Kirchhoff stress
    m_stress = mu1*(B_strain- delta1*B_kk/3.d0)/(eta**(5.d0/3.d0))+K1*(eta-1)*delta1
    stress(1) = m_stress(1,1)
    stress(2) = m_stress(2,2)
    stress(3) = m_stress(3,3)
    stress(4) = m_stress(1,2)
    stress(5) = m_stress(1,3)
    stress(6) = m_stress(2,3)
    tau = stress*J_f
    tau_kk = sum(tau(1:3))
   ! Compute the Tangent Stiffness D
!    call invert_small(B_strain,inv_B,J_b)
!    Bvec(1) = B_strain(1,1)
!    Bvec(2) = B_strain(2,2)
!    Bvec(3) = B_strain(3,3)
!    Bvec(4) = B_strain(1,2)
!    Bvec(5) = B_strain(1,3)
!    Bvec(6) = B_strain(2,3)
!    inv_Bvec(1) = inv_B(1,1)
!    inv_Bvec(2) = inv_B(2,2)
!    inv_Bvec(3) = inv_B(3,3)
!    inv_Bvec(4) = inv_B(1,2)
!    inv_Bvec(5) = inv_B(1,3)
!    inv_Bvec(6) = inv_B(2,3)
!    IinvB = spread(identity,dim=2,ncopies=6)*spread(inv_Bvec,dim=1,ncopies=6)
!    BinvB = spread(Bvec,dim=2,ncopies=6)    *spread(inv_Bvec,dim=1,ncopies=6)
!    II    = spread(identity,dim=2,ncopies=6)*spread(identity,dim=1,ncopies=6)
!    D = mu1/(J_f**(2.d0/3.d0))*delta2 + mu1/(J_f**(2.d0/3.d0))/3*(B_kk/3.d0*IinvB-II-BinvB)+K1*J_f*(J_f-1.d0/2.d0)*IinvB
!   !    Assemble B_bar
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!
!        B_bar = B
!        B_diff = (dNbardy - dNdy)/3.d0
!            do I = 1,3
!            B_bar(I,1:3*n_nodes-2:3)= B_bar(I,1:3*n_nodes-2:3) + B_diff(1:n_nodes,1)
!            B_bar(I,2:3*n_nodes-1:3)= B_bar(I,2:3*n_nodes-1:3) + B_diff(1:n_nodes,2)
!            B_bar(I,3:3*n_nodes :3) = B_bar(I,3:3*n_nodes :3) + B_diff(1:n_nodes,3)
!            enddo
!     !    Assemble matrix G
!    G = 0.d0
!    G(1,1) = Bvec(1)
!    G(2,2) = Bvec(2)
!    G(3,3) = Bvec(3)
!    G(4,1) = Bvec(4)
!    G(4,2) = Bvec(4)
!    G(5,1) = Bvec(5)
!    G(5,3) = Bvec(5)
!    G(6,2) = Bvec(6)
!    G(6,3) = Bvec(6)
!    G(1,4) = Bvec(4)
!    G(2,5) = Bvec(4)
!    G(4,4) = Bvec(3)
!    G(4,5) = Bvec(1)
!    G(5,4) = Bvec(6)
!    G(6,5) = Bvec(5)
!    G(1,6) = Bvec(5)
!    G(3,7) = Bvec(5)
!    G(4,6) = Bvec(6)
!    G(5,6) = Bvec(3)
!    G(5,7) = Bvec(1)
!    G(6,7) = Bvec(4)
!    G(2,8) = Bvec(6)
!    G(3,9) = Bvec(5)
!    G(4,8) = Bvec(5)
!    G(5,9) = Bvec(4)
!    G(6,8) = Bvec(3)
!    G(6,9) = Bvec(2)
!    G = 2*G
!     !    Assemble matrix B_star
!        B_star = 0.d0
!        B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
!        B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
!        B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
!        B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
!        B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
!        B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
!        B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
!        B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
!        B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)
!
!        B_diff = (dNbardy - dNdy)/3.d0
!            do I = 1,3
!            B_star(I,1:3*n_nodes-2:3)= B_star(I,1:3*n_nodes-2:3) + B_diff(1:n_nodes,1)
!            B_star(I,2:3*n_nodes-1:3)= B_star(I,2:3*n_nodes-1:3) + B_diff(1:n_nodes,2)
!            B_star(I,3:3*n_nodes :3) = B_star(I,3:3*n_nodes :3) + B_diff(1:n_nodes,3)
!            enddo
!      ! Compute Sigma and Q
!        S = reshape(matmul(transpose(B),tau),(/3,length_dof_array/3/))
!      do i = 1,n_nodes
!        Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!        Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
!        Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!        Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
!      end do
!        Sigma = Pmat*transpose(Smat)
!       do i = 1,n_nodes
!            Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
!            Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
!       end do
!            Q = Pmat*transpose(Pmat)


        sdev = tau
        sdev(1:3) = sdev(1:3)-tau_kk/3.d0
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + tau(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
    end do
  
    return
end subroutine fieldvars_hyperelast_3d

