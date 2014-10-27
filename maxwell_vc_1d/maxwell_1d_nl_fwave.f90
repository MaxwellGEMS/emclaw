! ===============================================================================
subroutine rp1(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ===============================================================================
!
!   This version outputs f-waves.

!   Riemann solver for the time dependent nonlinear Maxwell em equations in 1d,
!   in this case eps and mu dependend on time and position,
!   variable coefficients
!       kappa1(q,t,x)*(q1)_t +  (q2)_x  = -(eps)_t*(q1)
!       kappa2(q,t,x)*(q2)_t +  (q1)_x = -(mu)_t*(q2)
!
!   where   q1=q1, q2=q2, eps=f(x,t), and mu=g(x,t)
!           kappa1 = eps + 2*chi2_e*q1 + 3*chi3_e*q1^2
!           kappa2 = mu  + 2*chi2_m*q2 + 3*chi3_m*q2^2
!
!   aux(1,i) = eps(i)
!   aux(2,i) = mu(i)
!   aux(3,i) = eps_t(i)
!   aux(4,i) = mu_t(i)

!   function f(x_i,t_i) gives the permittivity value at the ith cell at step t_i, 
!   g(x_i,t_i) gives the permeability value at the ith cell at step t_i
!    
!   For RIP:   f(x_i,t_i)=f(x_i-v*t_i), and g(x_i,t_i)=g(x_i-v*t_i)
!   the system assumes the em functions to be some constant value + transient part

!   On input,   ql contains the state vector at the left edge of each cell
!               qr contains the state vector at the right edge of each cell

!   On output,  fwave contains the waves as jumps in f,
!               s the speeds,
!
!               amdq = A^- Delta q,
!               apdq = A^+ Delta q,
!

! Note that the ith Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    integer,          intent(in)  :: mx, num_ghost, maxnx, num_aux, num_eqn, num_waves

    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: fwave(num_eqn,num_waves,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) ::    s(num_waves,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: apdq(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: amdq(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer          :: i, m, nl
    double precision :: q1i, q1im, q2i, q2im
    double precision :: dq1, dq2, b1, b2
    double precision :: epsi, epsim, mui, muim
    double precision :: epsti, epstim, muti, mutim
    double precision :: kappa1, kappa2, eo, mo, zo, co, dx
    double precision :: c,z
    double precision :: chi2_e, chi2_m, chi3_e, chi3_m

    common /cparam/  dx, chi2_e, chi2_m, chi3_e, chi3_m, eo, mo, co, zo, nl

!   split the jump in q at each interface into waves
    z = zo
    c = co
    do i = 2-num_ghost, mx+num_ghost
        epsi   = auxl(1,i  )
        epsim  = auxr(1,i-1)
        mui    = auxl(2,i  )
        muim   = auxr(2,i-1)
        epsti  = auxl(3,i  )
        epstim = auxr(3,i-1)
        muti   = auxl(4,i  )
        mutim  = auxr(4,i-1)
        q1i    = ql(1,i  )
        q1im   = qr(1,i-1)
        q2i    = ql(2,i  )
        q2im   = qr(2,i-1)

        dq2 = (q2i - q2im)/eo
        dq1 = (q1i - q1im)/mo

        b1 = (dq1*z - dq2)/(2.d0*z)
        b2 = (dq1*z + dq2)/(2.d0*z)

        kappa1 = 0.5d0*(epsi + epsim)
        kappa2 = 0.5d0*(mui  + muim )

        if (nl.eq.1) then  
            kappa1 = kappa1 + 0.5d0*chi2_e*(q1i + q1im) + (3.d0/8.d0)*chi3_e*((q1i + q1im)**2)
            kappa2 = kappa2 + 0.5d0*chi2_m*(q2i + q2im) + (3.d0/8.d0)*chi3_m*((q2i + q2im)**2)
        end if

        fwave(1,1,i) = b1 *(-z) / kappa1
        fwave(2,1,i) = b1 / kappa2
        s(1,i) = -c

        fwave(1,2,i) = b2 *(z) / kappa1
        fwave(2,2,i) = b2 / kappa2
        s(2,i) = c
    end do

!   compute the leftgoing and rightgoing fluctuations:
!   Note s(1,i) < 0   and   s(2,i) > 0.

    do m=1,num_eqn
        do i = 2-num_ghost, mx+num_ghost
            amdq(m,i) = fwave(m,1,i)
            apdq(m,i) = fwave(m,2,i)
        end do
    end do

    return
end subroutine rp1