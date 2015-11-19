    subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables, &
                 updated_state_variables,dstrain,dRot,stress1,Ndelete)
    use Types
    use ParamIO
    use Globals, only : TIME,DTIME
    use Element_Utilities, only : invert_small

    implicit none

    integer, intent( in )     :: n_properties
    integer, intent( in )     :: n_state_variables

    real (prec), intent( in )   :: element_properties(n_properties)
    real (prec), intent( in )   :: initial_state_variables(n_state_variables)
    real (prec), intent( in )   :: dstrain(3,3)
    real (prec), intent( in )   :: dRot(3,3)

    real (prec), intent( out )  :: stress1(6)
    real (prec), intent( out )  :: updated_state_variables(n_state_variables)

    integer      ::  Ndelete,maxit,nit

    real (prec)  ::  E,nu,Y,E0_dot,m,q_1,q_2,q_3,f_N,e_N,s_N,f_c,f_F,f_Fbar                                 ! Material properties
    real (prec)  ::  stress0(6),ematrix,Vf,stress0_m(3,3),stress1_m(3,3)
    real (prec)  ::  dev_strain(3,3),dev_stress(3,3),stress_kk,strain_kk
    real (prec)  ::  delta1(3,3),S_star(3,3),p_star,Se_star,f_star
    real (prec)  ::  phi,phi_sqrt,dee,dev,C_Se,C_p,Se,p,A,F_1,F_2
    real (prec)  ::  dphi_dSe,dphi_dp,dA_ddee,dA_ddev,dphi_dSe_ddee,dphi_dSe_ddev
    real (prec)  ::  dphi_dp_ddee,dphi_dp_ddev,dF_1ddee,dF_1ddev,dF_2ddee,dF_2ddev
    real (prec)  ::  K(2,2),R(2,1),inv_K(2,2),J_K,dd(2,1),wnorm
    real (prec)  ::  tol,relax,err1
    real (prec)  ::  dematrix,Vf1,ematrix1

    delta1 = 0.d0
    delta1(1,1) = 1.d0
    delta1(2,2) = 1.d0
    delta1(3,3) = 1.d0

    E = element_properties(1)
    nu = element_properties(2)
    Y = element_properties(3)
    E0_dot = element_properties(4)
    m = element_properties(5)
    q_1 = element_properties(6)
    q_2 = element_properties(7)
    q_3 = element_properties(8)
    f_N = element_properties(9)
    e_N = element_properties(10)
    s_N = element_properties(11)
    f_c = element_properties(12)
    f_F = element_properties(13)
    f_Fbar = (q_1 + dsqrt(q_1**2 - q_3))/q_3


    stress0 = initial_state_variables(1:6)      ! Stress at start of increment
    ematrix = initial_state_variables(7)
    Vf = initial_state_variables(8)

!   Compute the elastic predictors for the deviatoric and hydrostatic stress

    stress0_m(1,1) = stress0(1)
    stress0_m(2,2) = stress0(2)
    stress0_m(3,3) = stress0(3)
    stress0_m(1,2) = stress0(4)
    stress0_m(2,1) = stress0(4)
    stress0_m(1,3) = stress0(5)
    stress0_m(3,1) = stress0(5)
    stress0_m(2,3) = stress0(6)
    stress0_m(3,2) = stress0(6)

    stress_kk = sum(stress0(1:3))
    strain_kk = dstrain(1,1) + dstrain(2,2) + dstrain(3,3)
    dev_stress = stress0_m - stress_kk*delta1/3.d0
    dev_strain = dstrain - strain_kk*delta1/3.d0

   ! S_star = dev_stress + E/(1+nu)*dev_strain + matmul(matmul(dRot,dev_stress),transpose(dRot))
    S_star = E/(1+nu)*dev_strain + matmul(matmul(dRot,dev_stress),transpose(dRot))
    p_star = stress_kk/3.d0 + E/3.d0/(1.d0-2.d0*nu)*strain_kk
    Se_star = dsqrt(3.d0/2.d0*sum(S_star*S_star))
    if (Vf .lt. f_c) then
        f_star = Vf
    else if (Vf .lt. f_F) then
        f_star = f_c + (f_Fbar -f_c)/(f_F - f_c)*(Vf - f_c)
    else
        f_star = f_Fbar
    end if
    C_Se = - 3.d0/2.d0*E/(1+nu)
    C_p = - E/3.d0/(1.d0-2.d0*nu)


    phi_sqrt = Se_star**2/Y**2+2*q_1*f_star*cosh(3.d0/2.d0*q_2*p_star/Y)-(1+q_3*f_star**2)


 !    write(*,*) Vf,f_c,f_F,Ndelete
 !  write(*,*) stress0,dstrain
   !write(*,*) stress0
    if (phi_sqrt .lt. 1.d-10) then
        stress1_m = S_star + p_star*delta1
        dematrix = 0.d0

    else

        dee = 0.d0
        dev = 0.d0

        tol = 1.d-8
        maxit = 20
        relax = 1.d0
        err1 = 1.d0
        nit = 0


        do while ((err1 .gt. tol) .and. (nit .lt. maxit))          ! Newton Raphson loop
        nit = nit + 1

        !calculate phi
        Se = Se_star + C_Se*dee
        p = p_star + C_p*dev

        phi_sqrt = Se**2/Y**2+2.d0*q_1*f_star*cosh(3.d0/2.d0*q_2*p/Y) - (1.d0+q_3*f_star**2)
        phi = dsqrt(phi_sqrt)
        !calculate F
        dphi_dSe = Se/phi/Y/Y
        dphi_dp = 3.d0/2.d0*q_1*q_2*f_star/phi/Y*sinh(3.d0/2.d0*q_2*p/Y)
        A = dsqrt((dphi_dSe)**2+2.d0/9.d0*(dphi_dp)**2)
        F_1 = A*(dee/DTIME/E0_dot)-dphi_dSe*phi**m
        F_2 = A*(dev/DTIME/E0_dot)-dphi_dp*phi**m
        !calculate the derivatives of F

        dA_ddee = -A/phi*dphi_dSe*C_Se
        dA_ddev = -A/phi*dphi_dp *C_p

        dphi_dSe_ddee = (1.d0/phi/Y**2 - Se/phi**2/Y**2*dphi_dSe)*C_Se;
        dphi_dSe_ddev = -1.d0/phi*dphi_dSe*dphi_dp*C_p
        dphi_dp_ddee = -1.d0/phi*dphi_dSe*dphi_dp*C_Se
        dphi_dp_ddev = (-1.d0/phi*dphi_dp + 3.d0/2.d0*q_1*q_2*f_star/phi/Y*cosh(3.d0/2.d0*q_2*p/Y)*3.d0/2.d0*q_2/Y)*C_p

        dF_1ddee = dA_ddee*(dee/DTIME/E0_dot)+A/DTIME/E0_dot-dphi_dSe_ddee*phi**m-dphi_dSe*(m*phi**(m-1)*dphi_dSe*C_Se)
        dF_1ddev = dA_ddev*(dee/DTIME/E0_dot)            -dphi_dSe_ddev*phi**m-dphi_dSe   *(m*phi**(m-1)*dphi_dp *C_p)
        dF_2ddee = dA_ddee*(dev/DTIME/E0_dot)            -dphi_dp_ddee *phi**m-dphi_dp    *(m*phi**(m-1)*dphi_dSe*C_Se)
        dF_2ddev = dA_ddev*(dev/DTIME/E0_dot)+A/DTIME/E0_dot-dphi_dp_ddev*phi**m -dphi_dp *(m*phi**(m-1)*dphi_dp *C_p)

        K(1,1) = dF_1ddee
        K(1,2) = dF_1ddev
        K(2,1) = dF_2ddee
        K(2,2) = dF_2ddev

        R(1,1) = -F_1
        R(2,1) = -F_2

        call invert_small(K,inv_K,J_K)
        dd = matmul(inv_K,R)

        dee = dee + relax*dd(1,1)
        dev = dev + relax*dd(2,1)

        !   Check convergence

        wnorm = dee*dee + dev*dev
        err1 = dd(1,1)*dd(1,1) + dd(2,1)*dd(2,1)
        err1 = dsqrt(err1/wnorm)
        end do

        write(IOW,'(//A25,I5,A26)') '    Inner NR Step number ',nit,'         For Inner NR     '
        write(IOW,'(A36,D12.5)')    '    Correction                      ',err1
        write(IOW,'(A36,D12.5)')    '    dee                             ',dee
        write(IOW,'(A36,D12.5)')    '    dev                             ',dev

        Se = Se_star + C_Se*dee
        p = p_star + C_p*dev

        phi_sqrt = Se**2/Y**2+2.d0*q_1*f_star*cosh(3.d0/2.d0*q_2*p/Y) - (1.d0+q_3*f_star**2)
        phi = dsqrt(phi_sqrt)
        !calculate F
        dphi_dSe = Se/phi/Y/Y
        dphi_dp = 3.d0/2.d0*q_1*q_2*f_star/phi/Y*sinh(3.d0/2.d0*q_2*p/Y)
        A = dsqrt((dphi_dSe)**2+2.d0/9.d0*(dphi_dp)**2)
        stress1_m = S_star - dee*E/(1+nu)*3.d0/2.d0*S_star/Se_star + (p_star - E/3.d0/(1.d0-2.d0*nu)*dev)*delta1
        dematrix = DTIME*E0_dot*phi**m/(1.d0-Vf)/A*(dphi_dSe*Se + 1.d0/3.d0*dphi_dp*p)
    end if
        ematrix1 = ematrix + dematrix
        Vf1 = 1.d0 + (Vf-1.d0)*exp(-dev) + f_N*dematrix/s_N/dsqrt(2.d0*PI)*exp(-1.d0/2.d0*((ematrix1-e_N)/s_N)**2)
        stress1(1) = stress1_m(1,1)
        stress1(2) = stress1_m(2,2)
        stress1(3) = stress1_m(3,3)
        stress1(4) = stress1_m(1,2)
        stress1(5) = stress1_m(1,3)
        stress1(6) = stress1_m(2,3)

        updated_state_variables(1:6) = stress1(1:6)
        updated_state_variables(7) = ematrix1
        updated_state_variables(8) = Vf1
    if (Vf1 .ge. f_F) then
        Ndelete = Ndelete +1
    end if

   write(*,*) Vf1, f_F, Ndelete
    return
    end subroutine gurson


