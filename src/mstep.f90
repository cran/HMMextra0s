subroutine mstep1d(n, m, nn, v, Z, R, hatpie, hatmu, hatsig)
    ! Estimate pie_{i}, mu_{ij} and sigma_{ij}
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    ! Arguments
    integer n, m, nn
    real(r8) v(nn,m), Z(nn), R(nn,n)
    real(r8) hatpie(m), hatmu(n,m), hatsig(n,m)
    ! Local variables
    integer i, j, k
    real(r8) vZ(nn), vdotZ, vsum

    do j = 1,m

        vsum = 0.0_r8
        vdotZ = 0.0_r8
        ! We assume here that Z(i) is element of {0.0,1.0}
        do i = 1,nn
            vZ(i) = v(i,j)*Z(i)
            vsum = vsum + v(i,j)
            vdotZ = vdotZ + vZ(i)
        end do

        hatpie(j) = vdotZ/vsum

        do k = 1,n
            hatmu(k,j) = dot_product(vZ(:),R(:,k))/vdotZ
            hatsig(k,j) = sqrt(dot_product(vZ(:),(R(:,k)-hatmu(k,j))**2)/vdotZ)
        end do

    end do

end subroutine mstep1d

subroutine mstep2d(n, m, nn, v, Z, R, hatpie, hatmu, hatsig)
    ! Estimate pie_{i}, mu_{ij} and sigma_{ij}
    implicit none
    integer, parameter :: r8 = selected_real_kind(15, 307)
    ! Arguments
    integer n, m, nn
    real(r8) v(nn,m), Z(nn), R(nn,n)
    real(r8) hatpie(m), hatmu(m,n), hatsig(n,n,m)
    ! Local variables
    integer i, j, k, l
    real(r8) vZ(nn), vdotZ, vsum

    do j = 1,m

        vsum = 0.0_r8
        vdotZ = 0.0_r8
        ! We assume here that Z(i) is element of {0.0,1.0}
        do i = 1,nn
            vZ(i) = v(i,j)*Z(i)
            vsum = vsum + v(i,j)
            vdotZ = vdotZ + vZ(i)
        end do

        hatpie(j) = vdotZ/vsum

        do k = 1,n
            hatmu(j,k) = dot_product(vZ(:),R(:,k))/vdotZ
            ! Variance-covariance matrix is symmetric, so compute one triangular part only
            ! and copy covariances into the other part - this also ensures that hatsig is
            ! always exactly symmetric
            do l = 1,k
                hatsig(k,l,j) = dot_product(vZ(:),(R(:,k)-hatmu(j,k))*(R(:,l)-hatmu(j,l)))/vdotZ
                if (k.ne.l) hatsig(l,k,j) = hatsig(k,l,j)
            end do
        end do

    end do

end subroutine mstep2d
