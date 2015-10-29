!     Subroutines to calculate the stress vector and the matrial tangent stiffness matrix



!==========================SUBROUTINE hypoelastic_tangent ==============================
subroutine hypoelastic_tangent_3D(strain,n_properties,element_properties,&       ! Input variables
                               stress,D)      ! Output variables                    ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: n_properties                                           ! # properties for the element
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file

    ! Local Variables
    integer      :: I,J

    real (prec)  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  stress_0, strain_0,mp_n, K       ! Material properties
    real (prec)  ::  delta1(6,6),delta2(6,6)
    real (prec)  ::  Et,Es,strain_e,stress_e
!    real (prec)  ::  sdev(6)
    real (prec)  ::  strain_kk
    real (prec)  ::  e(6),e_dyadic_e(6,6)

    !     El props are:


    stress_0 = element_properties(1)
    strain_0 = element_properties(2)
    mp_n   = element_properties(3)
	K = element_properties(4)
	
    D = 0.d0
    ! delta 1 and delta 2

    delta1 = 0.d0
    delta1(1,1) = 2.d0
    delta1(2,2) = 2.d0
    delta1(3,3) = 2.d0
    delta1(4,4) = 1.d0
    delta1(5,5) = 1.d0
    delta1(6,6) = 1.d0
    delta2 = 0.d0
    delta2(1:3,1:3) = 1.d0

    strain_kk = sum(strain(1:3))/3.d0
    e = strain
    e(1:3) = e(1:3) - strain_kk
    e(4:6) = e(4:6)/2.d0
    e_dyadic_e = spread(e,dim=2,ncopies=6)*spread(e,dim=1,ncopies=6)

    strain_e =dsqrt(dot_product(e(1:3),e(1:3)) + 2.d0*dot_product(e(4:6),e(4:6)))*dsqrt(2.d0/3.d0)
    if (strain_e .lt. strain_0) then
        stress_e = stress_0 *( dsqrt((1+mp_n**2)/(mp_n-1)**2 - (mp_n/(mp_n-1)-strain_e/strain_0)**2) - 1/(mp_n-1))
        Es = stress_e/strain_e
        Et = (stress_0/strain_0)/dsqrt((1+mp_n**2)/(mp_n-1)**2 -(mp_n/(mp_n-1)-strain_e/strain_0)**2)&
        *(mp_n/(mp_n-1)-strain_e/strain_0)
        if (all(abs(strain) == 0.d0)) then
        Es = Et
        end if
    else
        stress_e = stress_0 *(strain_e/strain_0)**(1/mp_n)
        Es = stress_e/strain_e
        Et = stress_0/strain_0/mp_n * (strain_e/strain_0)**(1/mp_n-1)
        if (all(abs(strain) ==0.d0)) then
        Es = Et
        end if
    endif

    D = (4.d0/9.d0/strain_e**2)*(Et-Es)*e_dyadic_e+Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2


    stress(1:3) = 2.d0/3.d0 *stress_e * e(1:3) /strain_e +3.d0*K*strain_kk
    stress(4:6) = 2.d0/3.d0 *stress_e * e(4:6) /strain_e


    if (all(abs(strain) == 0.d0)) then
    D = Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2
    stress = 0.d0
    end if

    return
end subroutine hypoelastic_tangent_3D


!==========================SUBROUTINE hypoelastic_tangent ==============================
subroutine hypoelastic_tangent_2D(strain,n_properties,element_properties,&       ! Input variables
                               stress,D)      ! Output variables                    ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: n_properties                                           ! # properties for the element
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file

    ! Local Variables
    integer      :: I,J

    real (prec)  ::  strain(3)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  stress_0, strain_0,mp_n, K       ! Material properties
    real (prec)  ::  delta1(3,3),delta2(3,3)
    real (prec)  ::  Et,Es,strain_e,stress_e
!    real (prec)  ::  sdev(6)
    real (prec)  ::  strain_kk
    real (prec)  ::  e(3),e_dyadic_e(3,3)

    !     El props are:


    stress_0 = element_properties(1)
    strain_0 = element_properties(2)
    mp_n   = element_properties(3)
    K = element_properties(4)

    D = 0.d0
    ! delta 1 and delta 2

    delta1 = 0.d0
    delta1(1,1) = 2.d0
    delta1(2,2) = 2.d0
    delta1(3,3) = 1.d0
    delta2 = 0.d0
    delta2(1:2,1:2) = 1.d0

    strain_kk = sum(strain(1:2))/3.d0
    e = strain
    e(1:2) = e(1:2) - strain_kk
    e(3) = e(3)/2.d0
    e_dyadic_e = spread(e,dim=2,ncopies=3)*spread(e,dim=1,ncopies=3)

    strain_e =dsqrt(dot_product(e(1:2),e(1:2)) + 2.d0*e(3)*e(3))*dsqrt(2.d0/3.d0)
    if (strain_e .lt. strain_0) then
        stress_e = stress_0 *( dsqrt((1+mp_n**2)/(mp_n-1)**2 - (mp_n/(mp_n-1)-strain_e/strain_0)**2) - 1/(mp_n-1))
        Es = stress_e/strain_e
        Et = (stress_0/strain_0)/dsqrt((1+mp_n**2)/(mp_n-1)**2 -(mp_n/(mp_n-1)-strain_e/strain_0)**2)&
        *(mp_n/(mp_n-1)-strain_e/strain_0)
        if (all(abs(strain) == 0.d0)) then
        Es = Et
        end if
    else
        stress_e = stress_0 *(strain_e/strain_0)**(1/mp_n)
        Es = stress_e/strain_e
        Et = stress_0/strain_0/mp_n * (strain_e/strain_0)**(1/mp_n-1)
        if (all(abs(strain) ==0.d0)) then
        Es = Et
        end if
    endif

    D = (4.d0/9.d0/strain_e**2)*(Et-Es)*e_dyadic_e+Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2


    stress(1:2) = 2.d0/3.d0 *stress_e * e(1:2) /strain_e +3.d0*K*strain_kk
    stress(3) = 2.d0/3.d0 *stress_e * e(3) /strain_e


    if (all(abs(strain) == 0.d0)) then
    D = Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2
    stress = 0.d0
    end if

    return
end subroutine hypoelastic_tangent_2D

!==========================SUBROUTINE hypoelastic_tangent ==============================
subroutine hypoelastic_tangent_2D_stress(strain,n_properties,element_properties,&       ! Input variables
                               stress,D)      ! Output variables                    ! Output variables
    use Types
    use ParamIO
    use Printparameters
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    implicit none

    integer, intent( in )         :: n_properties                                           ! # properties for the element
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file

    ! Local Variables
    integer      :: I,J,nit,maxit

    real (prec)  ::  strain(3)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(3)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(3,3)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  stress_0, strain_0,mp_n, K       ! Material properties
    real (prec)  ::  delta1(3,3),delta2(3,3)
    real (prec)  ::  Et,Es,strain_e,stress_e
!    real (prec)  ::  sdev(6)
    real (prec)  ::  strain_kk
    real (prec)  ::  e(3),e_dyadic_e(3,3)
    real (prec)  ::  e_33,strain_33,stress_33,dstrain_33,dstress_strain_33
    real (prec)  ::  tol,relax,err1



    !     El props are:



    stress_0 = element_properties(1)
    strain_0 = element_properties(2)
    mp_n   = element_properties(3)
    K = element_properties(4)

    D = 0.d0
    ! delta 1 and delta 2

    delta1 = 0.d0
    delta1(1,1) = 2.d0
    delta1(2,2) = 2.d0
    delta1(3,3) = 1.d0
    delta2 = 0.d0
    delta2(1:2,1:2) = 1.d0

    !Newton Raphson iteration
      tol = 0.0001d0
      maxit = 50
      relax = 0.7d0
      err1 =  1.d0
      nit = 0
      strain_33 = 0.d0
    if (all(abs(strain) == 0.d0)) then
    strain_33 = 0.d0

    else

    do while((err1 .gt. tol).and.(nit .lt. maxit))
    nit = nit + 1

    strain_kk = (sum(strain(1:2))+strain_33)/3.d0
    e = strain
    e(1:2) = e(1:2) - strain_kk
    e(3) = e(3)/2.d0
    e_33 = strain_33 - strain_kk
    e_dyadic_e = spread(e,dim=2,ncopies=3)*spread(e,dim=1,ncopies=3)

    strain_e =dsqrt(dot_product(e(1:2),e(1:2)) + 2.d0*e(3)*e(3)+e_33*e_33)*dsqrt(2.d0/3.d0)
    if (strain_e .lt. strain_0) then
        stress_e = stress_0 *( dsqrt((1+mp_n**2)/(mp_n-1)**2 - (mp_n/(mp_n-1)-strain_e/strain_0)**2) - 1/(mp_n-1))
        Es = stress_e/strain_e
        Et = (stress_0/strain_0)/dsqrt((1+mp_n**2)/(mp_n-1)**2 -(mp_n/(mp_n-1)-strain_e/strain_0)**2)&
        *(mp_n/(mp_n-1)-strain_e/strain_0)
        if (all(abs(strain) == 0.d0)) then
        Es = Et
        end if
    else
        stress_e = stress_0 *(strain_e/strain_0)**(1/mp_n)
        Es = stress_e/strain_e
        Et = stress_0/strain_0/mp_n * (strain_e/strain_0)**(1/mp_n-1)
        if (all(abs(strain) ==0.d0)) then
        Es = Et
        end if
    endif
    stress_33 = 2.d0/3.d0 *stress_e * strain_33 /strain_e +3.d0*K*strain_kk
    dstress_strain_33 = (4.d0/9.d0/strain_e**2)*(Et-Es)*strain_33**2+4.d0*Es/9.d0+K
    dstrain_33 = -stress_33/dstress_strain_33


    strain_33 = strain_33 + relax *dstrain_33
! Check convergence
    err1 = dstrain_33*dstrain_33
    err1 = dsqrt(err1)

    write(IOW,*)' Newton Raphson for strain_33 completed successfully '
    write(IOW,'(1x,D12.5,I5)') err1,nit

    end do



    D = (4.d0/9.d0/strain_e**2)*(Et-Es)*e_dyadic_e+Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2


    stress(1:2) = 2.d0/3.d0 *stress_e * e(1:2) /strain_e +3.d0*K*strain_kk
    stress(3) = 2.d0/3.d0 *stress_e * e(3) /strain_e

    end if

    if (all(abs(strain) == 0.d0)) then
    D = Es/3.d0*delta1+(K-2.d0*Es/9.d0)*delta2
    stress = 0.d0
    end if

    return
end subroutine hypoelastic_tangent_2D_stress

