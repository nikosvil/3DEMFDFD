	subroutine NS(q,nx,ny,nz,hx,hy,hz,ci,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +             tempEZ)

	implicit none
        integer :: nx,ny,nz,i
	real*8  :: hx,hy,hz,MultVal1(*),MultVal2(*),MultVal3(*),dznrm2,
     +  a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz)
	complex*16 :: tempEZ((nx-1)*(ny-1)*nz),ci,
     +  pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     +  pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     +  pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     +  pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     +  pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     +  pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     +  RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     +  RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     +  RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     +  sv1(nx*(ny-1)),sv2((nx-1)*ny),
     +  q((nx-1)*ny*(nz-1))
     
! ===================================================================
! Title: NS.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! NS.f performs the multiplication N*q=tempEZ
! by performing the operations between terms that make up block N 
! of matrix U and the input array q
!====================================================================      
! Input:
! q: complex array, dimension (nx-1)*ny*(nz-1)
! nx: number of elements in x-direction
! ny: number of elements in y-direction
! nz: number of elements in z-direction
! hx: discretization step in x-direction
! hy: discretization step in y-direction
! hz: discretization step in z-direction
! Output:
! tempEZ: complex array, dimension (nx-1)*(ny-1)*nz
!==================================================================== 

! B13*q=
!====================================================================
	  call B_mult(q,13,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S13*=
!====================================================================
	  call S_solve(13,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

! B20*pt7=
!====================================================================
	  call B_mult(pt7,20,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S20*=
!====================================================================
	  call S_solve(20,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

!$OMP PARALLEL DO
	  do i=1,(nx-1)*(ny-1)*nz
           tempEZ(i)=(24.0d0/(hy*hz))*pt6(i)
	  enddo
!$OMP END PARALLEL DO

	  return
	  end subroutine
