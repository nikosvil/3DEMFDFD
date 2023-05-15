      subroutine PS(q,nx,ny,nz,hx,hy,hz,ci,e3,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +             temp3,temp33,
     +             tempEZ)

      implicit none
      integer :: nx,ny,nz,i
      real*8  :: hx,hy,hz,MultVal1(*),MultVal2(*),MultVal3(*),dznrm2,
     +  a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz)
      complex*16 :: q((nx-1)*(ny-1)*nz),e3((nx-1)*(ny-1)*nz),ci,
     +  pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     +  pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     +  pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     +  pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     +  pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     +  pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     +  RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     +  RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     +  RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     +  tempEZ((nx-1)*(ny-1)*nz),
     +  temp3((nx-1)*(ny-1)*nz),temp33((nx-1)*(ny-1)*nz),
     +  sv1(nx*(ny-1)),sv2((nx-1)*ny)
     
! ===================================================================
! Title: PS.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! PS.f performs the multiplication P*q=tempEZ
! by performing the operations between terms that make up block P 
! of matrix U and the input array q
!====================================================================      
! Input:
! q: complex array, dimension (nx-1)*(ny-1)*nz
! nx: integer, number of elements in x-direction
! ny: integer, number of elements in y-direction
! nz: integer, number of elements in z-direction
! hx: real, discretization step in x-direction
! hy: real, discretization step in y-direction
! hz: real, discretization step in z-direction
! e3: complex array, dimension (nx-1)*(ny-1)*nz
! Output:
! tempEZ: complex array, dimension (nx-1)*(ny-1)*nz
!==================================================================== 

! B19*q=pt5
!====================================================================
      call B_mult(q,19,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S19*pt5=pt5
!====================================================================
      call S_solve(19,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

! B20*pt5=pt6
!====================================================================
      call B_mult(pt5,20,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6,pt7)
! S20*pt6=pt6
!====================================================================
      call S_solve(20,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

!$OMP PARALLEL DO
      do i=1,(nx-1)*(ny-1)*nz
        temp3(i)=(24.0d0/hy**2)*pt6(i)
      enddo
!$OMP END PARALLEL DO

! B17*q=pt3
!====================================================================
      call B_mult(q,17,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S17*pt3=pt3
!====================================================================
      call S_solve(17,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
! B18*pt3=
!====================================================================
      call B_mult(pt3,18,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6)
! S18*pt3=pt3
!====================================================================
      call S_solve(18,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,
     +             st1,st2,st3,st4,st5,st6,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)
     
!$OMP PARALLEL DO
      do i=1,(nx-1)*(ny-1)*nz
        temp33(i)=(24.0d0/hx**2)*pt6(i)
      enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
      do i=1,(nx-1)*(ny-1)*nz
        tempEZ(i)=e3(i)*q(i)-temp3(i)-temp33(i)
      enddo
!$OMP END PARALLEL DO

      return
      end
