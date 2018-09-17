      subroutine loop1(m, T, phi, pRS, gamma, logalp, lscale, tmp)
c     first loop (forward eqns)
      implicit none
      integer i, j, m, T
      double precision phi(m), sumphi, lscale
      double precision pRS(T,m), gamma(m,m), logalp(T,m)
      double precision tmp(m)
c     the above array occurs in the subroutine call for
c     memory allocation reasons in non gfortran compilers
c     its contents are purely internal to this subroutine
      lscale = 0
      i = 1
      do while(i .le. T)
          if (i .gt. 1) call multi1(m, phi, gamma, tmp)
          j = 1
          sumphi=0.0
          do while(j .le. m)
              phi(j) = phi(j)*pRS(i, j)
              sumphi = sumphi + phi(j)
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              phi(j) = phi(j)/sumphi
              j = j+1
          enddo
          lscale = lscale + dlog(sumphi)
          j = 1
          do while(j .le. m)
              logalp(i,j) = dlog(phi(j)) + lscale
              j = j+1
          enddo
          i = i+1
      enddo
      end

      subroutine multi1(m, a, b, c)
c     a is row vector
c     b is (m*m) matrix
c     a is replaced by the matrix product of a*b
      implicit none
      integer m, j, k
      double precision a(m), b(m, m), c(m)
      j = 1
      do while(j .le. m)
          k = 1
          c(j) = 0
          do while(k .le. m)
              c(j) = c(j) + a(k)*b(k, j)
              k = k+1
          enddo
          j = j+1
      enddo
      j = 1
      do while(j .le. m)
          a(j) = c(j)
          j = j+1
      enddo
      end


