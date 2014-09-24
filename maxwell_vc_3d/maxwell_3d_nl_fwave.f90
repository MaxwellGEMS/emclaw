! ===============================================================================
subroutine rpn3(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! ===============================================================================
!   Riemann solver for Maxwell's  equations in 3d, with time-space 
!   varying material properties.

!   in this case eps and mu dependend on time and position,
!   variable coefficients

!   Note that although there are 6 eigenvectors per direction, two 
!   eigenvalues are always zero and so we only need to compute 4 waves.

!   Solve Riemann problems along one slice of data.
!   This data is along a slice in the x-direction if ixyz=1
!                                 the y-direction if ixyz=2.
!                                 the z-direction if ixyz=3.

!   On output, fwave contains the waves as jumps in f, and
!              s the speeds,
!
!              amdq = A^- Delta q,
!              apdq = A^+ Delta q,
!                     the decomposition of the flux difference minus the source term
!                         f(qr(i-1)) - f(ql(i)) - psi(()q,x,t)
!                     into leftgoing and rightgoing parts respectively.
!

!   Note that the ith Riemann problem has left state qr(:,i-1)
!                                      and right state ql(:,i)
!   From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    integer,          intent(in)  :: ixyz, mx, num_ghost, maxnx, num_aux, num_eqn, num_waves

    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: fwave(num_eqn,num_waves,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) ::    s(num_waves,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: apdq(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(out) :: amdq(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer          :: i, m
    double precision :: q1i, q1im, q2i, q2im, q3i, q3im, q4i, q4im, q5i, q5im, q6i, q6im
    double precision :: dq1, dq2, dq3, dq4, dq5, dq6, b1, b2, b3, b4
    double precision :: df1, df2, df3, df4, df5, df6
    double precision :: eta1i, eta1im, eta2i, eta2im, eta3i, eta3im
    double precision :: eta4i, eta4im, eta5i, eta5im, eta6i, eta6im
    double precision :: etat1i, etat1im, etat2i, etat2im, etat3i, etat3im
    double precision :: etat4i, etat4im, etat5i, etat5im, etat6i, etat6im
    double precision :: kappa1, kappa2, kappa3, kappa4, kappa5, kappa6
    double precision :: eo, mo, zo, co, dx, dy, dz
    double precision :: ci, cim, zi, zim, z
    double precision :: chi2(6), chi3(6)

    common /cparam/  dx, dy, dz, chi2, chi3, zo, co, eo, mo

!   pre calculate variables that won't cahnge in loop
    ci  = co
    cim = co
    zi  = zo
    zim = zo
    z   = zi + zim

    fwave(:,:,:) = 0.0d0

!   split the jump in q at each interface into waves
    do i = 2-num_ghost, mx+num_ghost
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

        q1i     = ql(1,i  )
        q1im    = qr(1,i-1)
        q2i     = ql(2,i  )
        q2im    = qr(2,i-1)
        q3i     = ql(3,i  )
        q3im    = qr(3,i-1)
        q4i     = ql(4,i  )
        q4im    = qr(4,i-1)
        q5i     = ql(5,i  )
        q5im    = qr(5,i-1)
        q6i     = ql(6,i  )
        q6im    = qr(6,i-1)

!       flux difference
        dq1 = (q1i - q1im)
        dq2 = (q2i - q2im)
        dq3 = (q3i - q3im)
        dq4 = (q4i - q4im)
        dq5 = (q5i - q5im)
        dq6 = (q6i - q6im)

        kappa1 = 0.5d0*(eta1i + eta1im + 2.d0*chi2(1)*(q1i + q1im) + 3.d0*chi3(1)*((q1i + q1im)**2))
        kappa2 = 0.5d0*(eta2i + eta2im + 2.d0*chi2(2)*(q2i + q2im) + 3.d0*chi3(2)*((q2i + q2im)**2))
        kappa3 = 0.5d0*(eta3i + eta3im + 2.d0*chi2(3)*(q3i + q3im) + 3.d0*chi3(3)*((q3i + q3im)**2))
        kappa4 = 0.5d0*(eta4i + eta4im + 2.d0*chi2(4)*(q4i + q4im) + 3.d0*chi3(4)*((q4i + q4im)**2))
        kappa5 = 0.5d0*(eta5i + eta5im + 2.d0*chi2(5)*(q5i + q5im) + 3.d0*chi3(5)*((q5i + q5im)**2))
        kappa6 = 0.5d0*(eta6i + eta6im + 2.d0*chi2(6)*(q6i + q6im) + 3.d0*chi3(6)*((q6i + q6im)**2))

        if (ixyz == 1) then
            df2 = dq6/eo
            df3 = dq5/eo
            df5 = dq3/mo
            df6 = dq2/mo

            b1 = (-df2 + df6*zi )/z
            b2 = (-df3 - df5*zi )/z
            b3 = ( df3 - df5*zim)/z
            b4 = ( df2 + df6*zim)/z

            fwave(2,1,i) = b1*(-zim)/kappa2
            fwave(6,1,i) = b1/kappa6
            s(1,i) = -cim

            fwave(3,2,i) = b2*(zim)/kappa3
            fwave(5,2,i) = b2/kappa5
            s(2,i) = -cim

            fwave(3,3,i) = b3*(-zi)/kappa3
            fwave(5,3,i) = b3/kappa5
            s(3,i) = ci

            fwave(2,4,i) = b4*(zi)/kappa2
            fwave(6,4,i) = b4/kappa6
            s(4,i) = ci

        else if (ixyz == 2) then
            df1 = dq6/eo
            df3 = dq4/eo
            df4 = dq3/mo
            df6 = dq1/mo

            b1 = (-df1 - df6*zi )/z
            b2 = (-df3 + df4*zi )/z
            b3 = ( df3 + df4*zim)/z
            b4 = ( df1 - df6*zim)/z

            fwave(1,1,i) = b1*(zim)/kappa1
            fwave(6,1,i) = b1/kappa6
            s(1,i) = -cim

            fwave(3,2,i) = b2*(-zim)/kappa3
            fwave(4,2,i) = b2/kappa4
            s(2,i) = -cim

            fwave(3,3,i) = b3*(zi)/kappa3
            fwave(4,3,i) = b3/kappa6
            s(3,i) = ci

            fwave(1,4,i) = b4*(-zi)/kappa1
            fwave(6,4,i) = b4/kappa6
            s(4,i) = ci

        else if (ixyz == 3) then
            df1 = dq5/eo
            df2 = dq4/eo
            df4 = dq2/mo
            df5 = dq1/mo

            b1 = (-df1 + df5*zi )/z
            b2 = (-df2 - df4*zi )/z
            b3 = ( df2 - df4*zim)/z
            b4 = ( df1 + df5*zim)/z

            fwave(1,1,i) = b1*(-zim)/kappa1
            fwave(5,1,i) = b1/kappa5
            s(1,i) = -cim

            fwave(2,2,i) = b2*(zim)/kappa2
            fwave(4,2,i) = b2/kappa4
            s(2,i) = -cim

            fwave(2,3,i) = b3*(-zi)/kappa2
            fwave(4,3,i) = b3/kappa4
            s(3,i) = ci

            fwave(1,4,i) = b4*(zi)/kappa1
            fwave(5,4,i) = b4/kappa5
            s(4,i) = ci

        endif
    enddo

!   compute the leftgoing and rightgoing fluctuations:
!   Note s(1:2,i) < 0   and   s(3:4,i) > 0.

    do m=1,num_eqn
        do i = 2-num_ghost, mx+num_ghost
            amdq(m,i) = fwave(m,1,i) + fwave(m,2,i)
            apdq(m,i) = fwave(m,3,i) + fwave(m,4,i)
        end do
    end do

    return
end subroutine rpn3
