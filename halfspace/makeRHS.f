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
      
! Subroutine makeRHS computes the right hand side F of the main equation 
! using the conductivity vectors (SIGMAEX, SIGMAEY and SIGMAEZ) computed by makeS 
! and the primary electric field Ex, Ey values (ExBack and EyBack) computed by 
! create_EbackFiles_Halfspace Matlab script

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

      return
      end
