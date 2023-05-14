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
     
! ===================================================================
! Title: makeRHS.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! makeRHS.f computes the right hand side vectors F1 and F2 of equations
! AEx+BEy+CEz=F1 and DEx+ZEy+KEz=F2
! and coefficients e1,e2,e3 needed in subroutines AS.f,ZS.f,PS.f
! using the conductivity values computed by makeS.f 
! (in vectors SIGMAEX, SIGMAEY, SIGMAEZ and SIGMADiffEX,SIGMADiffEY) 
! and the background electric field values (in vectors ExBack and EyBack) 
! computed by create_EbackFiles_Halfspace Matlab script
!====================================================================      
! Input:
! nx: number of elements in x-direction
! ny: number of elements in x-direction
! nz: number of elements in x-direction
! omega: angular frequency
! m0: magnetic permeability 
! SIGMAEX: conductivity σ in Ex nodes, dimension: nx*(ny-1)*(nz-1)
! SIGMAEY: conductivity σ in Ey nodes, dimension: (nx-1)*ny*(nz-1) 
! SIGMAEZ: conductivity σ in Ez nodes, dimension: (nx-1)*(ny-1)*nz
! SIGMADiffEX: Difference between conductivity σ and background conductivity σ0 in Ex nodes, dimension: nx*(ny-1)*(nz-1)
! SIGMADiffEY: Difference between conductivity σ and background conductivity σ0 in Ey nodes, dimension: (nx-1)*ny*(nz-1) 
! ExBack: Background electric field intensity Ex, dimension: nx*(ny-1)*(nz-1)
! EyBack: Background electric field intensity Ey, dimension: (nx-1)*ny*(nz-1) 
! Output:
! e1: i*omega*m0*SIGMAEX, dimension: nx*(ny-1)*(nz-1)
! e2: i*omega*m0*SIGMAEY, dimension: (nx-1)*ny*(nz-1) 
! e3: i*omega*m0*SIGMAEZ, dimension: (nx-1)*(ny-1)*nz
! F1: Right-hand side of equation 1 (refers to Ex nodes), dimension: nx*(ny-1)*(nz-1)
! F2: Right-hand side of equation 1 (refers to Ey nodes), dimension: (nx-1)*ny*(nz-1) 
!==================================================================== 



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

!       Compute F1
!$OMP PARALLEL
!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       F1(i)=-ci*omega*m0*SIGMADiffEX(i)*ExBack(i)
      enddo
!$OMP END DO
!       Compute F2
!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       F2(i)=-ci*omega*m0*SIGMADiffEY(i)*EyBack(i)
      enddo
!$OMP END DO
!       Compute e1 coefficient used in subroutine AS.f
!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       e1(i) = ci*omega*m0*SIGMAEX(i)
      enddo
!$OMP END DO
!       Compute e2 coefficient used in subroutine ZS.f
!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       e2(i) = ci*omega*m0*SIGMAEY(i)
      enddo
!$OMP END DO
!       Compute e3 coefficient used in subroutine PS.f
!$OMP DO
      do i=1,(nx-1)*(ny-1)*nz
       e3(i) = ci*omega*m0*SIGMAEZ(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end
