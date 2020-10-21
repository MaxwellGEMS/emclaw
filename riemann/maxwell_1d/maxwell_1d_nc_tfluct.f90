! ============================================================================
subroutine tfluct1(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,amdq2)
! ============================================================================
!
!   "Internal" Riemann solver for Maxwell's equations in 1d.
!
!   On input, q     contains the cell average state vector
!             ql    contains the state vector at the left edge of each cell
!             qr    contains the state vector at the right edge of each cell
!             auxl  contains the auxiliary vector at the left edge of each cell
!             auxr  contains the state vector at the right edge of each cell
!             maxnx contains the number of physical points (without ghost cells)
!             num_ghost   contains the number of ghost cells
!             num_eqn  contains the number of equations 
!
!   On output, amdq2 contains the decomposition of the flux difference
!              f(qr(i)) - f(ql(i)).
!           _         _  
!   f(q) = |   q2/eo   |
!          |_  q1/mo  _|

    implicit none

    integer,          intent(in)  :: mx, num_ghost, maxnx, num_aux, num_eqn, num_waves

    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: amdq2(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer          :: i, nl, psi
    double precision :: q1i, q1im, q2i, q2im
    double precision :: df1, df2, psi1, psi2
    double precision :: epsi, epsim, mui, muim
    double precision :: epsti, epstim, muti, mutim
    double precision :: kappa1, kappa2, eo, mo, zo, co, dx
    double precision :: chi2_e, chi2_m, chi3_e, chi3_m

    common /cparam/  dx, chi2_e, chi2_m, chi3_e, chi3_m, eo, mo, co, zo, nl, psi

    psi1 = 0.d0
    psi2 = 0.d0

    do i = 1-num_ghost, mx+num_ghost
        q1i   = qr(1,i)
        q1im  = ql(1,i)
        q2i   = qr(2,i)
        q2im  = ql(2,i)

        epsi  = auxr(1,i)
        epsim = auxl(1,i)
        mui   = auxr(2,i)
        muim  = auxl(2,i)

        kappa1 = 0.5d0*(epsi + epsim)
        kappa2 = 0.5d0*(mui  + muim )

        if (nl.eq.1) then  
            kappa1 = kappa1 + 0.5d0*chi2_e*(q1i + q1im) + (3.d0/8.d0)*chi3_e*((q1i + q1im)**2)
            kappa2 = kappa2 + 0.5d0*chi2_m*(q2i + q2im) + (3.d0/8.d0)*chi3_m*((q2i + q2im)**2)
        end if

        df1 = (q2i - q2im)/eo
        df2 = (q1i - q1im)/mo

        if (psi.eq.1) then
            epsti  = auxr(3,i)
            epstim = auxl(3,i)
            muti   = auxr(4,i)
            mutim  = auxl(4,i)

            psi1 = -0.5d0*(epsti*q1i + epstim*q1im)
            psi2 = -0.5d0*(muti*q2i  + mutim*q2im )

            amdq2(1,i) = (df1 - dx*psi1)/kappa1
            amdq2(2,i) = (df2 - dx*psi2)/kappa2
        else
            amdq2(1,i) = df1/kappa1
            amdq2(2,i) = df2/kappa2
        end if

    enddo

    return
end subroutine tfluct1
