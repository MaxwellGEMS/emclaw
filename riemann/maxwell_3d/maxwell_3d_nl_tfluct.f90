! ==================================================================
subroutine tfluct3(ixyz,maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,ql,qr,auxl,auxr,amdq2)
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
!
!             num_ghost   contains the number of ghost cells
!             num_eqn     contains the number of equations
!
!   On output, amdq2 contains the decomposition of the flux difference
!              f(qr(i)) - f(ql)) - deltaxyz*psii.
!

    implicit none

    integer,          intent(in)  :: ixyz, mx, num_ghost, maxnx, num_aux, num_eqn, num_waves

    double precision, intent(in)  :: auxl(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  :: auxr(num_aux,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   ql(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in)  ::   qr(num_eqn,1-num_ghost:maxnx+num_ghost)

    double precision, intent(out) :: amdq2(num_eqn,1-num_ghost:maxnx+num_ghost)

    integer          :: i, nl, psi
    double precision :: q1i, q1im, q2i, q2im, q3i, q3im, q4i, q4im, q5i, q5im, q6i, q6im
    double precision :: dq1, dq2, dq3, dq4, dq5, dq6
    double precision :: df1, df2, df3, df4, df5, df6, psi1, psi2, psi3, psi4, psi5, psi6
    double precision :: eta1i, eta1im, eta2i, eta2im, eta3i, eta3im
    double precision :: eta4i, eta4im, eta5i, eta5im, eta6i, eta6im
    double precision :: etat1i, etat1im, etat2i, etat2im, etat3i, etat3im
    double precision :: etat4i, etat4im, etat5i, etat5im, etat6i, etat6im
    double precision :: kappa1, kappa2, kappa3, kappa4, kappa5, kappa6, zo, co, dx, dy, dz
    double precision :: eo, mo
    double precision :: chi2(6), chi3(6)

    common /cparam/  dx, dy, dz, chi2, chi3, co, zo, eo, mo, nl, psi

    amdq2(:,:) = 0.0d0

    do i = 2-num_ghost, mx+num_ghost
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

        if (psi.eq.1) then
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

            psi1 = -0.5d0*(etat1i*q1i + etat1im*q1im)
            psi2 = -0.5d0*(etat2i*q2i + etat2im*q2im)
            psi3 = -0.5d0*(etat3i*q3i + etat3im*q3im)
            psi4 = -0.5d0*(etat4i*q4i + etat4im*q4im)
            psi5 = -0.5d0*(etat5i*q5i + etat5im*q5im)
            psi6 = -0.5d0*(etat6i*q6i + etat6im*q6im)
        end if

!       flux difference
        dq1 = (q1i - q1im)
        dq2 = (q2i - q2im)
        dq3 = (q3i - q3im)
        dq4 = (q4i - q4im)
        dq5 = (q5i - q5im)
        dq6 = (q6i - q6im)

        kappa1 = 0.5d0*(eta1i + eta1im)
        kappa2 = 0.5d0*(eta2i + eta2im)
        kappa3 = 0.5d0*(eta3i + eta3im)
        kappa4 = 0.5d0*(eta4i + eta4im)
        kappa5 = 0.5d0*(eta5i + eta5im)
        kappa6 = 0.5d0*(eta6i + eta6im)

        if (nl.eq.1) then
            kappa1 = kappa1 + 0.5d0*chi2(1)*(q1i + q1im) + (3.d0/8.d0)*chi3(1)*((q1i + q1im)**2)
            kappa2 = kappa2 + 0.5d0*chi2(2)*(q2i + q2im) + (3.d0/8.d0)*chi3(2)*((q2i + q2im)**2)
            kappa3 = kappa3 + 0.5d0*chi2(3)*(q3i + q3im) + (3.d0/8.d0)*chi3(3)*((q3i + q3im)**2)
            kappa4 = kappa4 + 0.5d0*chi2(4)*(q4i + q4im) + (3.d0/8.d0)*chi3(4)*((q4i + q4im)**2)
            kappa5 = kappa5 + 0.5d0*chi2(5)*(q5i + q5im) + (3.d0/8.d0)*chi3(5)*((q5i + q5im)**2)
            kappa6 = kappa6 + 0.5d0*chi2(6)*(q6i + q6im) + (3.d0/8.d0)*chi3(6)*((q6i + q6im)**2)
        end if

        if (ixyz == 1) then
            df2 =  dq6/eo
            df3 = -dq5/eo
            df5 = -dq3/mo
            df6 =  dq2/mo

            if (psi.eq.1) then
                df2 = df2 - dx*psi2
                df3 = df3 - dx*psi3
                df5 = df5 - dx*psi5
                df6 = df6 - dx*psi6
            end if

            amdq2(2,i) = df2/kappa2
            amdq2(3,i) = df3/kappa3
            amdq2(5,i) = df5/kappa5
            amdq2(6,i) = df6/kappa6

        else if (ixyz == 2) then
            df1 = -dq6/eo
            df3 =  dq4/eo
            df4 =  dq3/mo
            df6 = -dq1/mo

            if (psi.eq.1) then
                df1 = df1 - dy*psi1
                df3 = df3 - dy*psi3
                df4 = df4 - dy*psi4
                df6 = df6 - dy*psi6
            end if

            amdq2(1,i) = df1/kappa1
            amdq2(3,i) = df3/kappa3
            amdq2(4,i) = df4/kappa4
            amdq2(6,i) = df6/kappa6

        else if (ixyz == 3) then
            df1 =  dq5/eo
            df2 = -dq4/eo
            df4 = -dq2/mo
            df5 =  dq1/mo

            if (psi.eq.1) then
                df1 = df1 - dz*psi1
                df2 = df2 - dz*psi2
                df4 = df4 - dz*psi4
                df5 = df5 - dz*psi5
            end if

            amdq2(1,i) = df1/kappa1
            amdq2(2,i) = df2/kappa2
            amdq2(4,i) = df4/kappa4
            amdq2(5,i) = df5/kappa5
        endif
    enddo

    return
end subroutine tfluct3
