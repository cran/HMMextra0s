subroutine loop2(m, T, phi, pRS, gamma, logbet, lscale, tmp)
    !     second loop (backward eqns) 
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    integer i, j, m, T
    real(r8) phi(m), sumphi, lscale, lscalearr(T-1)
    real(r8) pRS(T,m), gamma(m,m), logbet(T,m)
    real(r8) tmp(m)
    !     the above array occurs in the subroutine call for
    !     memory allocation reasons in non gfortran compilers
    !     its contents are purely internal to this subroutine
    do i = T-1,1,-1
        do j = 1,m
             phi(j) = phi(j)*pRS(i+1, j)
        enddo
        call multi2(m, gamma, phi, tmp)
        sumphi=0.0_r8
        do j = 1,m
            logbet(i,j) = phi(j)
            sumphi = sumphi + phi(j)
        enddo
        do j = 1,m
            phi(j) = phi(j)/sumphi
        enddo
        lscalearr(i) = lscale
        lscale = lscale + dlog(sumphi)
    enddo

    !     Separate this loop to invert loop nest order and enable
    !     vectorisation of inner loop
    do j = 1,m
        do i = 1,T-1
            logbet(i,j) = dlog(logbet(i,j)) + lscalearr(i)
        enddo
        logbet(T,j) = 0.0_r8
    enddo
end subroutine loop2

subroutine multi2(m, a, b, c)
    !     a is (m*m) matrix
    !     b is column vector
    !     b is replaced by the matrix product of a*b
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    integer m, j, k
    real(r8) a(m, m), b(m), c(m)
    do j = 1,m
        c(j) = 0.0_r8
        do k = 1,m
            c(j) = c(j) + a(j, k)*b(k)
        enddo
    enddo
    do j = 1,m
        b(j) = c(j)
    enddo
end subroutine multi2


