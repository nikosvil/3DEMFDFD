      subroutine U_Mult(b,n,nx,ny,nz,hx,hy,hz,ci,
     +      e1,e2,e3,tempEX,tempEY,tempEZ,
     +      pt1,pt2,pt3,pt4,pt5,pt6,
     +      st1,st2,st3,st4,st5,st6,
     +      RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +      MultVal1,MultVal2,MultVal3,
     +      a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +      temp1,temp11,temp2,temp22,temp3,temp33,
     +      tempEX_1,tempEX_2,tempEX_3,
     +      tempEY_1,tempEY_2,tempEY_3,
     +      tempEZ_1,tempEZ_2,tempEZ_3,
     +      temp)

        implicit none
        integer:: nx,ny,nz,i,n,k
        real*8:: hx,hy,hz,dznrm2,t1,t2,
     +  a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     +  MultVal1(5*int(dlog(dfloat(nx))/dlog(2.0d0))),
     +  MultVal2(5*int(dlog(dfloat(ny))/dlog(2.0d0))),
     +  MultVal3(5*int(dlog(dfloat(nz))/dlog(2.0d0)))
        complex*16:: b(n),ci,
     +  e1(nx*(ny-1)*(nz-1)),e2((nx-1)*ny*(nz-1)),e3((nx-1)*(ny-1)*nz),
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
     +  temp2((nx-1)*ny*(nz-1)),temp22((nx-1)*ny*(nz-1)),
     +  temp3((nx-1)*(ny-1)*nz),temp33((nx-1)*ny*(nz-1)),
     +  sv1(nx*(ny-1)),sv2((nx-1)*ny),
     +  tempEX(nx*(ny-1)*(nz-1)),  tempEY((nx-1)*ny*(nz-1)),
     +  tempEX_1(nx*(ny-1)*(nz-1)),tempEY_1((nx-1)*ny*(nz-1)),
     +  tempEX_2(nx*(ny-1)*(nz-1)),tempEY_2((nx-1)*ny*(nz-1)),
     +  tempEX_3(nx*(ny-1)*(nz-1)),tempEY_3((nx-1)*ny*(nz-1)),
     +  tempEZ((nx-1)*(ny-1)*nz),  tempEZ_1((nx-1)*(ny-1)*nz),
     +  tempEZ_2((nx-1)*(ny-1)*nz),tempEZ_3((nx-1)*(ny-1)*nz),
     +  temp(n)
! ===================================================================
! Title: U_Mult.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! U_Mult.f successively calls proper subroutines to perform the
! multiplication between the blocks of matrix U and the 
! respective vectors (as described in manuscript) during the
! implementation of the BiCGSTAB method
!====================================================================      
! Input:
! b: complex array, dimension n
! nx: number of elements in x-direction
! ny: number of elements in y-direction
! nz: number of elements in z-direction
! hx: discretization step in x-direction
! hy: discretization step in y-direction
! hz: discretization step in z-direction
! e1: complex array, dimension nx*(ny-1)*(nz-1), passed to subroutine AS.f
! e2: complex array, dimension (nx-1)*ny*(nz-1), passed to subroutine ZS.f
! e3: complex array, dimension (nx-1)*(ny-1)*nz, passed to subroutine PS.f
! Output:
! temp: complex array, dimension n
!====================================================================

c Line 1

        call AS(b(1:nx*(ny-1)*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,e1,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       temp1,temp11,
     +       tempEX_1)

        call BS(b(nx*(ny-1)*(nz-1)+1:nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEX_2)

        call CS(b(nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+1:n),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEX_3)

c Line 2
        call DS(b(1:nx*(ny-1)*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEY_1)

        call ZS(b(nx*(ny-1)*(nz-1)+1:nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,e2,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       temp2,temp22,
     +       tempEY_2)

        call KS(b(nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+1:n),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEY_3)

c Line 3
        call MS(b(1:nx*(ny-1)*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEZ_1)

        call NS(b(nx*(ny-1)*(nz-1)+1:nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)),
     +       nx,ny,nz,hx,hy,hz,ci,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       tempEZ_2)

        call PS(b(nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+1:n),
     +       nx,ny,nz,hx,hy,hz,ci,e3,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       temp3,temp33,
     +       tempEZ_3)

        call cpu_time(t1)
	
!$OMP PARALLEL DO
        do i=1,nx*(ny-1)*(nz-1)
         temp(i)=tempEX_1(i)+tempEX_2(i)+tempEX_3(i)
        enddo
!$OMP END PARALLEL DO
        k=1
        do i=nx*(ny-1)*(nz-1)+1,nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)
         temp(i)=tempEY_1(k)+tempEY_2(k)+tempEY_3(k)
         k=k+1
        enddo
        k=1
        do i=nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+1,n
         temp(i)=tempEZ_1(k)+tempEZ_2(k)+tempEZ_3(k)
         k=k+1
        enddo
      call cpu_time(t2)

	return
	end subroutine
