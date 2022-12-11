      subroutine makeRHS(nx,ny,nz,omega,m0,ci,wpath,
     +        SIGMAEX,SIGMAEY,SIGMAEZ,SIGMADiffEX,SIGMADiffEY,
     +        realExBack,imagExBack,realEyBack,imagEyBack,ExBack,EyBack,
     +        e1,e2,e3,F1,F2)

      implicit none
      integer:: nx,ny,nz,i
      character(64) :: wpath
      real*8:: omega,m0,dznrm2,
     +          SIGMAEX(nx*(ny-1)*(nz-1)),SIGMADiffEX(nx*(ny-1)*(nz-1)),
     +          SIGMAEY((nx-1)*ny*(nz-1)),SIGMADiffEY((nx-1)*ny*(nz-1)),
     +          SIGMAEZ((nx-1)*(ny-1)*nz),
     +          realExBack(nx*(ny-1)*(nz-1)),
     +          imagExBack(nx*(ny-1)*(nz-1)),
     +          realEyBack((nx-1)*ny*(nz-1)),
     +          imagEyBack((nx-1)*ny*(nz-1))
 
      complex*16:: ci,ExBack(nx*(ny-1)*(nz-1)),EyBack((nx-1)*ny*(nz-1)),
     +  e1(nx*(ny-1)*(nz-1)),e2((nx-1)*ny*(nz-1)),e3((nx-1)*(ny-1)*nz),
     +  F1(nx*(ny-1)*(nz-1)),F2((nx-1)*ny*(nz-1))
      
         open(10,file=trim(wpath)//'temp/ExBack.txt')
         open(11,file=trim(wpath)//'temp/EyBack.txt')
         do i=1,nx*(ny-1)*(nz-1)
          read(10,*) realExBack(i),imagExBack(i)
          ExBack(i)=dcmplx(realExBack(i),imagExBack(i))
         enddo
         do i=1,(nx-1)*ny*(nz-1)
          read(11,*) realEyBack(i),imagEyBack(i)
          EyBack(i)=dcmplx(realEyBack(i),imagEyBack(i))
         enddo
         close(10)
         close(11)
         print*,'ExBack - EyBack copied'

c       Compute F1
!$OMP PARALLEL
!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       F1(i)=-ci*omega*m0*SIGMADiffEX(i)*ExBack(i)
      enddo
!$OMP END DO
c       Compute F2
!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       F2(i)=-ci*omega*m0*SIGMADiffEY(i)*EyBack(i)
      enddo
!$OMP END DO

!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       e1(i) = ci*omega*m0*SIGMAEX(i)
      enddo
!$OMP END DO

!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       e2(i) = ci*omega*m0*SIGMAEY(i)
      enddo
!$OMP END DO

!$OMP DO
      do i=1,(nx-1)*(ny-1)*nz
       e3(i) = ci*omega*m0*SIGMAEZ(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL

       if (1.eq.0) then 
       open(12, file='F1.txt')
       open(13, file='F2.txt')
       open(14, file='e1.txt')
       open(15, file='e2w.txt')
       open(16, file='e3.txt')
       open(17, file='sigmaEX.txt')
       open(18, file='sigmaEY.txt')
       open(19, file='sigmaEZ.txt')
       open(20, file='sigmaDiffEX.txt')
       open(21, file='sigmaDiffEY.txt')

       do i=1,nx*(ny-1)*(nz-1)
        write(12,*) F1(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(13,*) F2(i)
       enddo
       do i=1,nx*(ny-1)*(nz-1)
        write(14,*) e1(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(15,*) e2(i)
       enddo
       do i=1,(nx-1)*(ny-1)*nz
        write(16,*) e3(i)
       enddo
       do i=1,nx*(ny-1)*(nz-1)
        write(17,*) sigmaEX(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(18,*) sigmaEY(i)
       enddo
       do i=1,(nx-1)*(ny-1)*nz
        write(19,*) sigmaEZ(i)
       enddo
       do i=1,nx*(ny-1)*(nz-1)
        write(20,*) sigmaDiffEX(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(21,*) sigmaDiffEY(i)
       enddo

       close(12)
       close(13)
       close(14)
       close(15)
       close(16)
       close(17)
       close(18)
       close(19)
       close(20)
       close(21)
       stop
       endif

      return
      end
