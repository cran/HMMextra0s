subroutine prsloop(m, nn, pie, R, mu, sig, Z, pRS)
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    ! sqrt(2*pi)
    real(r8), parameter :: sqrt2pi = sqrt(8.0_r8 * atan(1.0_r8))
    ! Arguments
    integer m, nn
    real(r8) pie(m), R(nn), mu(m), sig(m), Z(nn)
    real(r8) pRS(nn,m)
    ! Local variables
    integer j, k
    real(r8) oneminpie, fac1, fac2, deltaR

    do k = 1,m
        oneminpie = 1.0_r8-pie(k)
        fac1 = pie(k)/(sqrt2pi*sig(k))
        fac2 = -1.0_r8/(2.0_r8*sig(k)*sig(k))
        do j = 1,nn
            deltaR = R(j) - mu(k)
            pRS(j,k) = (fac1*exp(fac2*deltaR*deltaR))*Z(j) + oneminpie*(1.0_r8-Z(j))
        end do
    end do

end subroutine prsloop
