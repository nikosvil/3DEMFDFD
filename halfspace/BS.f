        subroutine BS(q,nx,ny,nz,hx,hy,hz,ci,
     +              pt1,pt2,pt3,pt4,pt5,pt6,pt7,
     +              st1,st2,st3,st4,st5,st6,st7,
     +              RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,
     +              MultVal1,MultVal2,MultVal3,
     +              a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +              tempEX)

        implicit none
        integer :: nx,ny,nz,i
        real*8  :: hx,hy,hz,
     +  a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     +  MultVal1(*),MultVal2(*),MultVal3(*),dznrm2
        complex*16 :: tempEX(nx*(ny-1)*(nz-1)),ci,
     +  pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     +  pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     +  pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     +  pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     +  pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     +  pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     +  pt7((nx-1)*(ny-1)*(nz-1)),st7((nx-1)*(ny-1)*(nz-1)),
     +  RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)),
     +  RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     +  RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     +  RHS7((nx-1)*(ny-1)*(2*nz-1)),
     +  sv1(nx*(ny-1)),sv2((nx-1)*ny),
     +  q((nx-1)*ny*(nz-1))

! ===================================================================
! Title: BS.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! BS.f performs the multiplication B*q=tempEX
! by performing the operations between terms that make up block B 
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
! tempEX: complex array, dimension nx*(ny-1)*(nz-1)
!==================================================================== 

        call B_mult(q,1,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6,pt7)
        call S_solve(1,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,pt7,
     +             st1,st2,st3,st4,st5,st6,st7,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

        call B_mult(pt2,2,nx,ny,nz,pt1,pt2,pt3,pt4,pt5,pt6,pt7)

        call S_solve(4,nx,ny,nz,
     +             pt1,pt2,pt3,pt4,pt5,pt6,pt7,
     +             st1,st2,st3,st4,st5,st6,st7,
     +             RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,
     +             MultVal1,MultVal2,MultVal3,
     +             a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)

!$OMP PARALLEL DO
      do i=1,nx*(ny-1)*(nz-1)
       tempEX(i)=(24.0d0/(hx*hy))*pt1(i)
      enddo
!$OMP END PARALLEL DO

      return
      end
