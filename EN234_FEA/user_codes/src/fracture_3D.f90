subroutine fracture_3d_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
        !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
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

    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint,I,J,Nindex,Nsize,Ndelete,jj

    real (prec)  ::  x(3,length_coord_array/3)                                       !  Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  u(3,length_dof_array/3)                                         !  Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  du(3,length_dof_array/3)                                        !  Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  delta(3,3)                                                      ! Constant matrix
    real (prec)  ::  dxidx(3,3), determinant                                         ! Jacobian inverse and determinant
    real (prec)  ::  F(3,3),dF(3,3),midF(3,3),inv_midF(3,3),J_f                      ! deformation gradient
    real (prec)  ::  dL(3,3),dL_kk,deta,eta,el_vol,dLbar(3,3)
    real (prec)  ::  dstrain(3,3),dspin(3,3),inv_Rot(3,3),J_Rot,dRot(3,3),stress1(6)
    real (prec)  ::  B(6,length_dof_array),B_bar(6,length_dof_array),B_diff(20,3)    ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  D(6,6),dstress(1:6),dstrainl(1:6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)

    x = reshape(element_coords,(/3,length_coord_array/3/))
    du = reshape(dof_increment,(/3,length_dof_array/3/))
    u = reshape(dof_total,(/3,length_dof_array/3/))
    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27



    call initialize_integration_points(n_points, n_nodes, xi, w)

    ! constant matrix

    delta = 0.d0
    delta(1,1) = 1.d0
    delta(2,2) = 1.d0
    delta(3,3) = 1.d0

    element_residual = 0.d0
    Ndelete = 0.d0
    deta = 0.d0
    eta = 0.d0
    el_vol = 0.d0
    dNbardy =0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

   !  Compute the deformation gradients by differentiating displacements
       F = delta + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
       dF = matmul(du(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
       midF = F + dF/2.d0

 !      write(*,*) midf
          call invert_small(midF,inv_midF,J_f)
   !    Convert the shape function derivatives to derivatives with respect to deformed coordinates
       dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_midF)
   !    Add contribution to dNbardy and J_bar from current integration point
       dNbardy(1:n_nodes,1:3) = dNbardy(1:n_nodes,1:3) + dNdy(1:n_nodes,1:3)*J_f*w(kint)*determinant
!    Add contribution to the volume averaged Jacobian and the average volumetric strain increment
        dL = matmul(dF,inv_midF)
        dL_kk = dL(1,1) + dL(2,2) + dL(3,3)
            deta = deta + J_f*dL_kk*w(kint)*determinant
            eta = eta + J_f*w(kint)*determinant
            el_vol = el_vol + w(kint)*determinant
    enddo
            eta = eta/el_vol
            deta = deta/eta/el_vol
            dNbardy = dNbardy/eta/el_vol

    !     --  Loop over integration points a second time
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
   !  Compute the deformation gradients by differentiating displacements
       F = delta + matmul(u(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
       dF = matmul(du(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))
       midF = F + dF/2.d0
       call invert_small(midF,inv_midF,J_f)
      dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),inv_midF)
  ! Compute the corrected velocity gradient increment
        dL = matmul(dF,inv_midF)
        dL_kk = dL(1,1) + dL(2,2) + dL(3,3)
        dLbar = dL + (deta - dL_kk)/3.d0*delta
   !  Compute the strain and spin increments
      dstrain = (dLbar + transpose(dLbar))/2.d0
      dspin = (dLbar - transpose(dLbar))/2.d0
      call invert_small((delta - dspin/2.d0),inv_Rot,J_Rot)
      dRot = matmul(inv_Rot,(delta + dspin/2.d0))

      Nsize = 8
      Nindex = (kint-1)*Nsize + 1
   ! write(*,*)  dF
      call gurson(element_properties,n_properties,Nsize,initial_state_variables(Nindex:(Nindex+Nsize-1)), &
                updated_state_variables(Nindex:(Nindex+Nsize-1)),dstrain,dRot,stress1,Ndelete)

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

         element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_bar),stress1)*w(kint)*determinant


    end do

    if (Ndelete ==  n_points) then
    element_deleted = .true.
    else
    element_deleted = .false.
    endif
    return
end subroutine fracture_3d_dynamic

!------------------------------------------------------------------------------------------------------------------------------------

subroutine fieldvars_fracture_3d(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    logical      :: strcmp

    integer      :: n_points,kint,k,Nindex,Nsize

    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  ematrix, Vf
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  p, smises                          ! Pressure and Mises stress

    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)

      Nsize = 8
      Nindex = (kint-1)*Nsize + 1

    stress = updated_state_variables(Nindex:Nindex+5)
    ematrix = updated_state_variables(Nindex+6)
    Vf = updated_state_variables(Nindex+7)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt(dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)

        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do

    end do

    return
end subroutine fieldvars_fracture_3d


