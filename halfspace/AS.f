! ===================================================================
! Title: AS.f 
! Authors: N. D. Vilanakis, E. N. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
! * A 3D frequency-domain electromagnetic solver employing a high order compact finite-difference scheme
! * N. D. Vilanakis, N. Economou, E. N. Mathioudakis, A. Vafidis
! * Computers and Geosciences
!====================================================================
! AS.f performs the multiplication A*q=tempEX
! by performing the operations between terms that make up block A 
! of matrix U and the input array q
!====================================================================      
! Input:
! q: complex array, dimension nx*(ny-1)*(nz-1)
! nx: integer, number of elements in x-direction
! ny: integer, number of elements in y-direction
! nz: integer, number of elements in z-direction
! hx: real, discretization step in x-direction
! hy: real, discretization step in y-direction
! hz: real, discretization step in z-direction
! e1: complex array, dimension nx*(ny-1)*(nz-1)
! Output:
! tempEX: complex array, dimension nx*(ny-1)*(nz-1)
!==================================================================== 
	subroutine AS(q,nx,ny,nz,hx,hy,hz,ci,e1,
     +              pt1,pt2,pt3,pt4,pt5,pt6,
     +              st1,st2,st3,st4,st5,st6,
     +              RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +              MultVal1,MultVal2,MultVal3,
     +              a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +              temp1,temp11,
     +              tempEX)
        implicit none
        integer :: nx,ny,nz,i
        real*8  :: hx,hy,hz,dznrm2,
     +  a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     +  MultVal1(5*int(dlog(dfloat(nx))/dlog(2.0d0))),
     +  MultVal2(5*int(dlog(dfloat(ny))/dlog(2.0d0))),
     +  MultVal3(5*int(dlog(dfloat(nz))/dlog(2.0d0)))
        complex*16 :: q(nx*(ny-1)*(nz-1)),e1(nx*(ny-1)*(nz-1)),
     +  tempEX(nx*(ny-1)*(nz-1)),ci,
     +  pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     +  pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     +  pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     +  pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     +  pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     +  pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     +  RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     +  RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     +  RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     +  temp1(nx*(ny-1)*(nz-1)),temp11(nx*(ny-1)*(nz-1)),
     +  sv1(nx*(ny-1)),sv2((nx-1)*ny),zdotc
     
! B5*q=pt3
!====================================================================  
        call B_mult(q,5,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S5*pt3=pt3
!====================================================================  	
        call S_solve(5,nx,ny,nz,
     +           pt1,pt2,pt3,pt4,pt5,pt6,
     +           st1,st2,st3,st4,st5,st6,
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
     
! B6*pt3=pt1
!====================================================================  
        call B_mult(pt3,6,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S8*pt1=pt1
!====================================================================  
        call S_solve(8,nx,ny,nz,
     +           pt1,pt2,pt3,pt4,pt5,pt6,
     +           st1,st2,st3,st4,st5,st6,
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
     
!$OMP PARALLEL DO
        do i=1,nx*(ny-1)*(nz-1)
         temp1(i)=(24.0d0/(hz**2))*pt1(i)
        enddo
!$OMP END PARALLEL DO

! B3*q=pt2
!==================================================================== 
        call B_mult(q,3,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S3*pt2=pt2
!====================================================================
        call S_solve(3,nx,ny,nz,
     +           pt1,pt2,pt3,pt4,pt5,pt6,
     +           st1,st2,st3,st4,st5,st6,
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
     
! B4*pt2=pt1
!====================================================================
        call B_mult(pt2,4,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S4*pt1=pt1
!====================================================================
        call S_solve(4,nx,ny,nz,
     +           pt1,pt2,pt3,pt4,pt5,pt6,
     +           st1,st2,st3,st4,st5,st6,
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
     
!$OMP PARALLEL DO
        do i=1,nx*(ny-1)*(nz-1)
         temp11(i)=(24.0d0/(hy**2))*pt1(i)
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
        do i=1,nx*(ny-1)*(nz-1)
         tempEX(i)=e1(i)*q(i)-temp1(i)-temp11(i)
        enddo
!$OMP END PARALLEL DO
 
	  return
	  end
