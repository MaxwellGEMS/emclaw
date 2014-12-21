! ============================================================================
subroutine tfluct1(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,amdq2)
! ============================================================================
!
!   "Internal" Riemann solver for Maxwell's equations in 1d conservative form.
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
!   f(q) = |   q2/mu   |
!          |_  q1/eps _|

    implicit none

    integer,          intent(in)  :: mx, num_ghost, maxnx, num_aux, num_eqn, num_waves

    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: amdq2(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer          :: i
    double precision :: q1i, q1im, q2i, q2im
    double precision :: df1, df2
    double precision :: epsi, epsim, mui, muim
    double precision :: eo, mo, zo, co

    common /cparam/  eo, mo, co, zo, nl

    do i = 1-num_ghost, mx+num_ghost
        q1i   = qr(1,i)
        q1im  = ql(1,i)
        q2i   = qr(2,i)
        q2im  = ql(2,i)

        epsi  = auxr(1,i)
        epsim = auxl(1,i)
        mui   = auxr(2,i)
        muim  = auxl(2,i)

        df1 = q2i/mui  - q2im/muim
        df2 = q1i/epsi - q1im/epsim

        if (nl.eq.1) then
            df1 = df1 - chi3_m*((q2i**3)/mui  - (q2im**3)/muim)
            df2 = df2 - chi3_e*((q1i**3)/epsi - (q1im**3)/epsim)
        end if

        amdq2(1,i) = df1
        amdq2(2,i) = df2
    enddo

    return
end subroutine tfluct1
