! ===============================================================================
subroutine rp1(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ===============================================================================
!
!   This version outputs f-waves.

!   Riemann solver for the time dependent nonlinear Maxwell equations in 1d conservative form,
!   in this case eps and mu dependend on time and position,
!   variable coefficients
!       (q1)_t +  (q2/eta2)_x = 0
!       (q2)_t +  (q1/eta1)_x = 0
!
!   where   q1=eta1*E, q2=eta2*H, eps=f(x,t), and mu=g(x,t)
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

! Note that the ith Riemann problem has left  state qr(:,i-1)
!                                   and right state ql(:,i)
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

    integer          :: i, m
    double precision :: q1i, q1im, q2i, q2im
    double precision :: df1, df2, b1, b2
    double precision :: epsi, epsim, mui, muim
    double precision :: eo, mo, zo, co, dx
    double precision :: ci, cim, zi, zim
    double precision :: chi3_e, chi3_m

    common /cparam/  dx, eo, mo, co, zo, nl

!   split the jump in q at each interface into waves

    do i = 2-num_ghost, mx+num_ghost
        epsi   = auxl(1,i  )
        epsim  = auxr(1,i-1)
        mui    = auxl(2,i  )
        muim   = auxr(2,i-1)
        q1i    = ql(1,i  )
        q1im   = qr(1,i-1)
        q2i    = ql(2,i  )
        q2im   = qr(2,i-1)

        zi  = sqrt(epsi/mui)
        zim = sqrt(epsim/muim)

        ci  = 1.0d0/sqrt(epsi*mui)
        cim = 1.0d0/sqrt(epsim*muim)

        df1 = q2i/mui  - q2im/muim
        df2 = q1i/epsi - q1im/epsim

        if (nl.eq.1) then
            zi  = zi*sqrt((1.0d0 - 3.0d0*chi3_m*(q2i**2))/(1.0d0 - 3.0d0*chi3_e*(q1i**2))) 
            zim = zim*sqrt((1.0d0 - 3.0d0*chi3_m*(q2im**2))/(1.0d0 - 3.0d0*chi3_e*(q1im**2))) 

            ci  = ci*sqrt((1.0d0 - 3.0d0*chi3_m*(q2i**2))*(1.0d0 - 3.0d0*chi3_e*(q1i**2)))
            cim = cim*sqrt((1.0d0 - 3.0d0*chi3_m*(q2im**2))*(1.0d0 - 3.0d0*chi3_e*(q1im**2)))

            df1 = df1 - chi3_m*((q2i**3)/mui  - (q2im**3)/muim)
            df2 = df2 - chi3_e*((q1i**3)/epsi - (q1im**3)/epsim)
        end if

        b1 = (df2*zi  - df1)/(zi + zim)
        b2 = (df2*zim + df1)/(zi + zim)

        fwave(1,1,i) = b1*(-zim)
        fwave(2,1,i) = b1
        s(1,i) = -cim

        fwave(1,2,i) = b2*(zi)
        fwave(2,2,i) = b2
        s(2,i) = ci
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