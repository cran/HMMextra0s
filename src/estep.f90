subroutine estep(m, nn, logalpha, logbeta, ll, pRS, gam, v, w)
    ! Calculate v_t(j) and w_t(i,j)
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    ! Arguments
    integer m, nn
    real(r8) logalpha(nn, m), logbeta(nn, m), ll, pRS(nn, m), gam(m,m)
    real(r8) v(nn, m), w(nn-1, m, m)
    ! Local variables
    integer j, k
    real(r8) loggamminll(m,m), logPbeta(nn-1)

    loggamminll(:,:) = log(gam(:,:)) - ll

    do k = 1,m

        logPbeta(:) = log( pRS(2:nn,k) ) + logbeta(2:nn,k)
        v(:,k) = exp(logalpha(:,k) + logbeta(:,k) - ll)

        do j = 1,m
            w(:,j,k) = exp(loggamminll(j,k) + logalpha(1:nn-1,j) + &
                           logPbeta(:))
        end do

    end do

end subroutine estep
