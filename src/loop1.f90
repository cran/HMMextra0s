subroutine loop1(m, T, phi, pRS, gamma, logalp, lscale, tmp)
    !     first loop (forward eqns)
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    integer i, j, m, T
    real(r8) phi(m), sumphi, lscale, lscalearr(T)
    real(r8) pRS(T,m), gamma(m,m), logalp(T,m)
    real(r8) tmp(m)
    !     the above array occurs in the subroutine call for
    !     memory allocation reasons in non gfortran compilers
    !     its contents are purely internal to this subroutine
    lscale = 0.0_r8
    do i = 1,T
        if (i .gt. 1) call multi1(m, phi, gamma, tmp)
        sumphi=0.0_r8
        do j = 1,m
            phi(j) = phi(j)*pRS(i, j)
            sumphi = sumphi + phi(j)
        enddo
        do j = 1,m
            phi(j) = phi(j)/sumphi
            logalp(i,j) = phi(j)
        enddo
        lscale = lscale + dlog(sumphi)
        lscalearr(i) = lscale
    enddo

    !     Separate this loop to invert loop nest order and enable
    !     vectorisation of inner loop
    do j = 1,m
        do i = 1,T
            logalp(i,j) = dlog(logalp(i,j)) + lscalearr(i)
        enddo
    enddo
end subroutine loop1

subroutine multi1(m, a, b, c)
    !     a is row vector
    !     b is (m*m) matrix
    !     a is replaced by the matrix product of a*b
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    integer m, j, k
    real(r8) a(m), b(m, m), c(m)
    do j = 1,m
        c(j) = 0.0_r8
        do k = 1,m
            c(j) = c(j) + a(k)*b(k, j)
        enddo
    enddo
    do j = 1,m
        a(j) = c(j)
    enddo
end subroutine multi1


