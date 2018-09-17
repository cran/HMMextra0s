      subroutine loop2(m, T, phi, pRS, gamma, logbet, lscale, tmp)
c     second loop (backward eqns) 
      implicit none
      integer i, j, m, T
      double precision phi(m), sumphi, lscale
      double precision pRS(T,m), gamma(m,m), logbet(T,m)
      double precision tmp(m)
c     the above array occurs in the subroutine call for
c     memory allocation reasons in non gfortran compilers
c     its contents are purely internal to this subroutine
      i = T-1
      do while(i .ge. 1)
          j = 1
          do while(j .le. m)
              phi(j) = phi(j)*pRS(i+1, j)
              j = j+1
          enddo
          call multi2(m, gamma, phi, tmp)
          j = 1
          sumphi=0.0
          do while(j .le. m)
              logbet(i,j) = dlog(phi(j)) + lscale
              sumphi = sumphi + phi(j)
              j = j+1
          enddo
          j = 1
          do while(j .le. m)
              phi(j) = phi(j)/sumphi
              j = j+1
          enddo
          lscale = lscale + dlog(sumphi)
          i = i-1
      enddo
      end

      subroutine multi2(m, a, b, c)
c     a is (m*m) matrix
c     b is column vector
c     b is replaced by the matrix product of a*b
      implicit none
      integer m, j, k
      double precision a(m, m), b(m), c(m)
      j = 1
      do while(j .le. m)
          k = 1
          c(j) = 0
          do while(k .le. m)
              c(j) = c(j) + a(j, k)*b(k)
              k = k+1
          enddo
          j = j+1
      enddo
      j = 1
      do while(j .le. m)
          b(j) = c(j)
          j = j+1
      enddo
      end


