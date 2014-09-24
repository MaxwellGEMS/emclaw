! ==================================================================
subroutine tfluct3(ixyz,maxnx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,amdq2)
! ==================================================================
!
!   "Internal" Riemann solver for Maxwell's equations in 3D. 
!
!   On input, q     contains the cell average state vector
!             ql    contains the state vector at the left edge of each cell
!             qr    contains the state vector at the right edge of each cell
!             auxl  contains the auxiliary vector at the left edge of each cell
!             auxr  contains the state vector at the right edge of each cell
!             maxnx contains the number of physical points (without ghost cells)
!             mbc   contains the number of ghost cells
!             meqn  contains the number of equations
!
!   On output, amdq2 contains the decomposition of the flux difference
!              f(qr(i)) - f(ql)) - deltaxyz*psi(i).
!

    implicit none

    double precision :: auxl(maux,1-mbc:maxnx+mbc)
    double precision :: auxr(maux,1-mbc:maxnx+mbc)
    double precision ::   ql(meqn,1-mbc:maxnx+mbc)
    double precision ::   qr(meqn,1-mbc:maxnx+mbc)
    double precision :: amdq2(meqn,1-mbc:maxnx+mbc)
    double precision :: chi2, chi3, vac(6)
    double precision :: psi(6), kappa(6), dq(6)
    integer :: i, mx, mbc, maxnx, maux, meqn, mwaves, ixyz

    double precision :: eta1i, eta1im, eta2i, eta2im, eta3i, eta3im
    double precision :: eta4i, eta4im, eta5i, eta5im, eta6i, eta6im
    double precision :: etat1i, etat1im, etat2i, etat2im, etat3i, etat3im
    double precision :: etat4i, etat4im, etat5i, etat5im, etat6i, etat6im
    double precision :: q1i, q1im, q2i, q2im, q3i, q3im
    double precision :: q4i, q4im, q5i, q5im, q6i, q6im
    double precision :: eo, mo
    double precision :: df1, df2, df3, df4, df5, df6
    double precision :: dx, dy, dz
    
    common /cparam/  dx, dy, dz, chi2(6), chi3(6), eo, mo

!     split the jump in q at each interface into waves
    vac(1:3) = mo
    vac(4:6) = eo
    do i = 2-mbc, mx+mbc
        eta1i   = auxr(1,i)
        eta1im  = auxl(1,i)
        eta2i   = auxr(2,i)
        eta2im  = auxl(2,i)
        eta3i   = auxr(3,i)
        eta3im  = auxl(3,i)
        eta4i   = auxr(4,i)
        eta4im  = auxl(4,i)
        eta5i   = auxr(5,i)
        eta5im  = auxl(5,i)
        eta6i   = auxr(6,i)
        eta6im  = auxl(6,i)
        etat1i  = auxr(7,i)
        etat1im = auxl(7,i)
        etat2i  = auxr(8,i)
        etat2im = auxl(8,i)
        etat3i  = auxr(9,i)
        etat3im = auxl(9,i)
        etat4i  = auxr(10,i)
        etat4im = auxl(10,i)
        etat5i  = auxr(11,i)
        etat5im = auxl(11,i)
        etat6i  = auxr(12,i)
        etat6im = auxl(12,i)

        q1i     = qr(1,i)
        q1im    = ql(1,i)
        q2i     = qr(2,i)
        q2im    = ql(2,i)
        q3i     = qr(3,i)
        q3im    = ql(3,i)
        q4i     = qr(4,i)
        q4im    = ql(4,i)
        q5i     = qr(5,i)
        q5im    = ql(5,i)
        q6i     = qr(6,i)
        q6im    = ql(6,i)

        psi(1) = -0.5d0*(etat1i*q1i + etat1im*q1im)
        psi(2) = -0.5d0*(etat2i*q2i + etat2im*q2im)
        psi(3) = -0.5d0*(etat3i*q3i + etat3im*q3im)
        psi(4) = -0.5d0*(etat4i*q4i + etat4im*q4im)
        psi(5) = -0.5d0*(etat5i*q5i + etat5im*q5im)
        psi(6) = -0.5d0*(etat6i*q6i + etat6im*q6im)        

!     flux difference
        dq(1) = (q1i - q1im)
        dq(2) = (q2i - q2im)
        dq(3) = (q3i - q3im)
        dq(4) = (q4i - q4im)
        dq(5) = (q5i - q5im)
        dq(6) = (q6i - q6im)

        kappa(1) = 0.5d0*(eta1i + eta1im + 2.d0*chi2(1)*(q1i + q1im) + 3.d0*chi3(1)*((q1i + q1im)**2))
        kappa(2) = 0.5d0*(eta2i + eta2im + 2.d0*chi2(2)*(q2i + q2im) + 3.d0*chi3(2)*((q2i + q2im)**2))
        kappa(3) = 0.5d0*(eta3i + eta3im + 2.d0*chi2(3)*(q3i + q3im) + 3.d0*chi3(3)*((q3i + q3im)**2))
        kappa(4) = 0.5d0*(eta4i + eta4im + 2.d0*chi2(4)*(q4i + q4im) + 3.d0*chi3(4)*((q4i + q4im)**2))
        kappa(5) = 0.5d0*(eta5i + eta5im + 2.d0*chi2(5)*(q5i + q5im) + 3.d0*chi3(5)*((q5i + q5im)**2))
        kappa(6) = 0.5d0*(eta6i + eta6im + 2.d0*chi2(6)*(q6i + q6im) + 3.d0*chi3(6)*((q6i + q6im)**2))


        if (ixyz == 1) then
            df2 = dq(6)/vac(2) - dx*psi(2)
            df3 = dq(5)/vac(3) - dx*psi(3)
            df5 = dq(3)/vac(5) - dx*psi(5)
            df6 = dq(2)/vac(6) - dx*psi(6)
            
            amdq2(1,i) = 0.0d0
            amdq2(2,i) = df2/kappa(2)
            amdq2(3,i) = df3/kappa(3)
            amdq2(4,i) = 0.0d0
            amdq2(5,i) = df5/kappa(5)
            amdq2(6,i) = df6/kappa(6)

        else if (ixyz == 2) then
            df1 = dq(6)/vac(1) - dy*psi(1)
            df3 = dq(4)/vac(3) - dy*psi(3)
            df4 = dq(3)/vac(4) - dy*psi(4)
            df6 = dq(1)/vac(6) - dy*psi(6)

            amdq2(1,i) = df1/kappa(1)
            amdq2(2,i) = 0.0d0
            amdq2(3,i) = df3/kappa(3)
            amdq2(4,i) = df4/kappa(4)
            amdq2(5,i) = 0.0d0
            amdq2(6,i) = df6/kappa(6)

        else if (ixyz == 3) then
            df1 = dq(5)/vac(1) - dz*psi(1)
            df2 = dq(4)/vac(2) - dz*psi(2)
            df4 = dq(2)/vac(4) - dz*psi(4)
            df5 = dq(1)/vac(5) - dz*psi(5)

            amdq2(1,i) = df1/kappa(1)
            amdq2(2,i) = df2/kappa(2)
            amdq2(3,i) = 0.0d0
            amdq2(4,i) = df4/kappa(4)
            amdq2(5,i) = df5/kappa(5)
            amdq2(6,i) = 0.0d0

        endif
    enddo

    return

    end subroutine tfluct3

