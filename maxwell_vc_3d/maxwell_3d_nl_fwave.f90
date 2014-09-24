!     ==================================================================
    subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!     ==================================================================

!   # Riemann solver for Maxwell's  equations in 3d, with time-space 
!   # varying material properties.

!   # in this case eps and mu dependend on time and position,
!   # variable coefficients

!   # Note that although there are 6 eigenvectors per direction, two 
!   # eigenvalues are always zero and so we only need to compute 4 waves.

!   # Solve Riemann problems along one slice of data.
!   # This data is along a slice in the x-direction if ixyz=1
!   #                               the y-direction if ixyz=2.
!   #                               the z-direction if ixyz=3.

!   # On output, fwave contains the waves as jumps in f, and
!   #            s the speeds,
!   #
!   #            amdq = A^- Delta q,
!   #            apdq = A^+ Delta q,
!   #                   the decomposition of the flux difference minus the source term
!   #                       f(qr(i-1)) - f(ql(i)) - psi(()q,x,t)
!   #                   into leftgoing and rightgoing parts respectively.
!   #

!   # Note that the ith Riemann problem has left state qr(:,i-1)
!   #                                    and right state ql(:,i)
!   # From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    double precision :: auxl(maux,1-mbc:maxm+mbc)
    double precision :: auxr(maux,1-mbc:maxm+mbc)
    double precision :: fwave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision ::    s(mwaves,1-mbc:maxm+mbc)
    double precision ::   ql(meqn,1-mbc:maxm+mbc)
    double precision ::   qr(meqn,1-mbc:maxm+mbc)
    double precision :: apdq(meqn,1-mbc:maxm+mbc)
    double precision :: amdq(meqn,1-mbc:maxm+mbc)
    double precision :: chi2, chi3, vac(6)
    double precision :: psi(6), beta(4), kappa(6), dq(6)
    integer :: i, mx, mbc, maxm, maux, meqn, mwaves, m, ixyz

    double precision :: eta1i, eta1im, eta2i, eta2im, eta3i, eta3im
    double precision :: eta4i, eta4im, eta5i, eta5im, eta6i, eta6im
    double precision :: etat1i, etat1im, etat2i, etat2im, etat3i, etat3im
    double precision :: etat4i, etat4im, etat5i, etat5im, etat6i, etat6im
    double precision :: q1i, q1im, q2i, q2im, q3i, q3im
    double precision :: q4i, q4im, q5i, q5im, q6i, q6im
    double precision :: ci, cim, zi, zim
!    double precision :: kappa1, kappa2, kappa3, kappa4, kappa5, kappa6 
    double precision :: eo, mo, zo, co
    double precision :: df1, df2, df3, df4, df5, df6
!    double precision :: psi1, psi2, psi3, psi4, psi5, psi6
    double precision :: dx, dy, dz
!    double precision :: dq1, dq2, dq3, dq4, dq5, dq6
!    double precision :: beta1, beta2, beta3, beta4
    common /cparam/  dx, dy, dz, chi2(6), chi3(6), eo, mo, co, zo

!     # split the jump in q at each interface into waves
    vac(1:3) = mo
    vac(4:6) = eo
    do i = 2-mbc, mx+mbc
        eta1i   = auxl(1,i  )
        eta1im  = auxr(1,i-1)
        eta2i   = auxl(2,i  )
        eta2im  = auxr(2,i-1)
        eta3i   = auxl(3,i  )
        eta3im  = auxr(3,i-1)
        eta4i   = auxl(4,i  )
        eta4im  = auxr(4,i-1)
        eta5i   = auxl(5,i  )
        eta5im  = auxr(5,i-1)
        eta6i   = auxl(6,i  )
        eta6im  = auxr(6,i-1)
        etat1i  = auxl(7,i  )
        etat1im = auxr(7,i-1)
        etat2i  = auxl(8,i  )
        etat2im = auxr(8,i-1)
        etat3i  = auxl(9,i  )
        etat3im = auxr(9,i-1)
        etat4i  = auxl(10,i  )
        etat4im = auxr(10,i-1)
        etat5i  = auxl(11,i  )
        etat5im = auxr(11,i-1)
        etat6i  = auxl(12,i  )
        etat6im = auxr(12,i-1)

        q1i     = ql(1,i)
        q1im    = qr(1,i-1)
        q2i     = ql(2,i)
        q2im    = qr(2,i-1)
        q3i     = ql(3,i)
        q3im    = qr(3,i-1)
        q4i     = ql(4,i)
        q4im    = qr(4,i-1)
        q5i     = ql(5,i)
        q5im    = qr(5,i-1)
        q6i     = ql(6,i)
        q6im    = qr(6,i-1)

!     # calculate velocity c = 1/sqrt(eps*mu) and impedance Z = sqrt(eps/mu)
        ci      = co 
        cim     = co
        zi      = zo
        zim     = zo

        psi(1) = -0.5d0*(etat1i*q1i + etat1im*q1im)
        psi(2) = -0.5d0*(etat2i*q2i + etat2im*q2im)
        psi(3) = -0.5d0*(etat3i*q3i + etat3im*q3im)
        psi(4) = -0.5d0*(etat4i*q4i + etat4im*q4im)
        psi(5) = -0.5d0*(etat5i*q5i + etat5im*q5im)
        psi(6) = -0.5d0*(etat6i*q6i + etat6im*q6im)        

!     # flux difference
        
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
            
            beta(1) = ( df2 + df6*zi )/(zi + zim)
            beta(2) = (-df3 + df5*zi )/(zi + zim)
            beta(3) = (-df2 + df6*zim)/(zi + zim)
            beta(4) = ( df3 + df5*zim)/(zi + zim)
            
            fwave(1,1,i)   = 0.0d0
            fwave(2,1,i)   = beta(1)*(zim)/kappa(2)
            fwave(3:5,1,i) = 0.0d0
            fwave(6,1,i)   = beta(1)/kappa(6)

            fwave(1:2,2,i) = 0.0d0
            fwave(3,2,i)   = beta(2)*(-zim)/kappa(2)
            fwave(4,2,i)   = 0.0d0
            fwave(5,2,i)   = beta(2)/kappa(5)
            fwave(6,2,i)   = 0.0d0

            fwave(1,3,i)   = 0.0d0
            fwave(2,3,i)   = beta(3)*(-zi)/kappa(2)
            fwave(3:5,3,i) = 0.0d0
            fwave(6,3,i)   = beta(3)/kappa(6)

            fwave(1:2,4,i) = 0.0d0
            fwave(3,4,i)   = beta(4)*(zi)/kappa(2)
            fwave(4,4,i)   = 0.0d0
            fwave(5,4,i)   = beta(4)/kappa(5)
            fwave(6,4,i)   = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        else if (ixyz == 2) then
            df1 = dq(6)/vac(1) - dy*psi(1)
            df3 = dq(4)/vac(3) - dy*psi(3)
            df4 = dq(3)/vac(4) - dy*psi(4)
            df6 = dq(1)/vac(6) - dy*psi(6)

            beta(1) = (-df1 + df6*zi )/(zi + zim)
            beta(2) = ( df3 + df4*zi )/(zi + zim)
            beta(3) = ( df1 + df6*zim)/(zi + zim)
            beta(4) = (-df3 + df4*zim)/(zi + zim)

            fwave(1,1,i)   = beta(1)*(-zim)/kappa(1)
            fwave(2:5,1,i) = 0.0d0
            fwave(6,1,i)   = beta(1)/kappa(6)

            fwave(1:2,2,i) = 0.0d0
            fwave(3,2,i)   = beta(2)*(zim)/kappa(3)
            fwave(4,2,i)   = beta(2)/kappa(4)
            fwave(5:6,2,i) = 0.0d0

            fwave(1,3,i)   = beta(3)*(zi)/kappa(1)
            fwave(2:5,3,i) = 0.0d0
            fwave(6,3,i)   = beta(3)/kappa(6)

            fwave(1:2,4,i) = 0.0d0
            fwave(3,4,i)   = beta(4)*(-zi)/kappa(3)
            fwave(4,4,i)   = beta(4)/kappa(4)
            fwave(5:6,4,i) = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        else if (ixyz == 3) then
            df1 = dq(5)/vac(1) - dz*psi(1)
            df2 = dq(4)/vac(2) - dz*psi(2)
            df4 = dq(2)/vac(4) - dz*psi(4)
            df5 = dq(1)/vac(5) - dz*psi(5)

            beta(1) = ( df1 + df5*zi )/(zi + zim)
            beta(2) = (-df2 + df4*zi )/(zi + zim)
            beta(3) = (-df1 + df5*zim)/(zi + zim)
            beta(4) = ( df2 + df4*zim)/(zi + zim)

            fwave(1,1,i)   = beta(1)*(zim)/kappa(1)
            fwave(2:4,1,i) = 0.0d0
            fwave(5,1,i)   = beta(1)/kappa(5)
            fwave(6,1,i)   = 0.0d0

            fwave(1,2,i)   = 0.0d0
            fwave(2,2,i)   = beta(2)*(-zim)/kappa(2)
            fwave(3,2,i)   = 0.0d0
            fwave(4,2,i)   = beta(2)/kappa(4)
            fwave(5:6,2,i) = 0.0d0

            fwave(1,3,i)   = beta(1)*(-zi)/kappa(1)
            fwave(2:4,3,i) = 0.0d0
            fwave(5,3,i)   = beta(1)/kappa(5)
            fwave(6,3,i)   = 0.0d0

            fwave(1,4,i)   = 0.0d0
            fwave(2,4,i)   = beta(2)*(zi)/kappa(2)
            fwave(3,4,i)   = 0.0d0
            fwave(4,4,i)   = beta(2)/kappa(4)
            fwave(5:6,4,i) = 0.0d0

            s(1:2,i) = -cim
            s(3:4,i) = ci

        endif
    
    enddo


!   # compute the leftgoing and rightgoing fluctuations:
!   # Note s(1:2,i) < 0   and   s(3:4,i) > 0.

    do m=1,meqn
        do i = 2-mbc, mx+mbc
            amdq(m,i) = fwave(m,1,i) + fwave(m,2,i)
            apdq(m,i) = fwave(m,3,i) + fwave(m,4,i)
        end do
    end do

    return

    end subroutine rpn3

