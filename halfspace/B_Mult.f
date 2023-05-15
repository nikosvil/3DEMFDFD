! ===================================================================
! Title: B_Mult.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! B_Mult.f performs the multiplication between each Bi matrix from the B class
! and the respective input array t
! Each time the subroutine is called, integer input tag specifies
! which operation is being performed using the input t array with the 
! proper dimensions and returning the respective tnew array
!====================================================================      
! Input:
! t: complex array, non-fixed dimension
! tag: integer index which indicates the respective B matrix
! nx: number of elements in x-direction
! ny: number of elements in y-direction
! nz: number of elements in z-direction
! Output:
! *One of the following arrays depending on the tag index*
! tnew1: complex array, dimension nx*(ny-1)*(nz-1)
! tnew2: complex array, dimension nx*ny*(nz-1)
! tnew3: complex array, dimension nx*(ny-1)*nz
! tnew4: complex array, dimension (nx-1)*ny*(nz-1)
! tnew5: complex array, dimension (nx-1)*ny*nz
! tnew6: complex array, dimension(nx-1)*(ny-1)*nz
!==================================================================== 
      subroutine B_Mult(t,tag,nx,ny,nz,
     +                  tnew1,tnew2,tnew3,tnew4,tnew5,tnew6)
c B_mult subroutine performs the multiplication of Matrix B_tag with
c vector t producing vector tnew
      implicit none
      integer :: nx,ny,nz,i,j,k,tag
      complex*16 :: t(*),
     +  tnew1(nx*(ny-1)*(nz-1)),
     +  tnew2(nx*ny*(nz-1)),
     +  tnew3(nx*(ny-1)*nz),
     +  tnew4((nx-1)*ny*(nz-1)),
     +  tnew5((nx-1)*ny*nz),
     +  tnew6((nx-1)*(ny-1)*nz),

!==================================================================== 
! B1*t=tnew2 tnew2: nx*ny*(nz-1), t: (nx-1)*ny*(nz-1)
!==================================================================== 
      if (tag.eq.1) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz-1
       do j=1,ny
        tnew2(1+(i-1)*ny*nx+(j-1)*nx)=-10.0d0*t((i-1)*ny*(nx-1)+(j-1)*
     +  (nx-1)+1)+36.0d0*t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+2)-
     +  6.0d0*t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+3)+
     +  t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+4)
        do k=2,nx-1
         tnew2(k+(i-1)*ny*(nx)+(j-1)*(nx))=24.0d0*(t(k+(i-1)*ny*(nx-1)+
     +   (j-1)*(nx-1))-t((k-1)+(i-1)*ny*(nx-1)+(j-1)*(nx-1)))
        enddo
        tnew2(nx+(i-1)*ny*nx+(j-1)*nx)=
     +  10.0d0*t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+(nx-1))-
     +  36.0d0*t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+(nx-2))+
     +  6.0d0*t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+(nx-3))-
     +  t((i-1)*ny*(nx-1)+(j-1)*(nx-1)+(nx-4))
       enddo
      enddo
!$OMP END PARALLEL DO
!==================================================================== 
! B2*t=tnew1 tnew1: nx*(ny-1)*(nz-1) t: nx*ny*(nz-1)
!==================================================================== 
      elseif (tag.eq.2) then

!$OMP PARALLEL DO COLLAPSE(3)
      do i=1,nz-1
       do j=1,ny-1
        do k=1,nx
         tnew1(k+(j-1)*nx+(i-1)*(ny-1)*nx)=-t(k+(j-1)*nx+(i-1)*(ny)*nx)
     +   + t(k+(j-1)*nx+(i-1)*(ny)*nx+nx)
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

c**********************************************
c B3*t=tnew2 tnew2:nx*ny*(nz-1) t:nx*(ny-1)*(nz-1)

      elseif (tag.eq.3) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz-1
       do j=1,ny
        if (j.eq.1) then
         do k=1,nx
         tnew2((i-1)*ny*nx+k)=-10.0d0*t((i-1)*(ny-1)*nx+k)
     +   +36.0d0*t((i-1)*(ny-1)*nx+nx+k)-6.0d0*
     +   t((i-1)*(ny-1)*nx+2*nx+k)+t((i-1)*(ny-1)*nx+3*nx+k)
         enddo
        elseif (j.eq.ny) then
         do k=1,nx
          tnew2((ny-1)*nx+(i-1)*ny*nx+k)=10.0d0*t((ny-2)*nx+(i-1)*
     +    (ny-1)*nx+k)-36.0d0*t((ny-3)*nx+(i-1)*(ny-1)*nx+k)+6.0d0*
     +    t((ny-4)*nx+(i-1)*(ny-1)*nx+k)-t((ny-5)*nx+
     +    (i-1)*(ny-1)*nx+k)
         enddo
        else
         do k=1,nx
          tnew2((j-1)*nx+(i-1)*ny*nx+k)=-24.0d0*(t((j-2)*nx+
     +    (i-1)*(ny-1)*nx+k)-t((j-1)*nx+(i-1)*(ny-1)*nx+k))
         enddo
        endif
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B4*t=tnew1 tnew1: nx*(ny-1)*(nz-1) t: nx*ny*(nz-1)
!==================================================================== 
      elseif (tag.eq.4) then
!$OMP PARALLEL DO COLLAPSE(3)
       do i=1,nz-1
        do j=1,ny-1
         do k=1,nx
          tnew1(k+(j-1)*nx+(i-1)*(ny-1)*nx) = -t(k+(j-1)*nx+(i-1)*
     +    (ny)*nx)+t(k+(j-1)*nx+(i-1)*(ny)*nx+nx)
          enddo
        enddo
       enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B5*t=tnew3 tnew3: nx*(ny-1)*nz t: nx*(ny-1)*(nz-1)
!==================================================================== 
      elseif (tag.eq.5) then
       do i = 1,nz
        if (i.eq.1) then
!$OMP PARALLEL DO 
         do k=1,nx*(ny-1)
          tnew3(k)=-10.0d0*t(k)+36.0d0*t((ny-1)*nx+k)-
     +    6.0d0*t(2*(ny-1)*nx+k)+t(3*(ny-1)*nx+k)
         enddo
!$OMP END PARALLEL DO
        elseif (i.eq.nz)then
!$OMP PARALLEL DO 
          do k=1,nx*(ny-1)
           tnew3((nz-1)*(ny-1)*nx+k)=10.0d0*t((nz-2)*(ny-1)*nx+k)-
     +     36.0d0*t((nz-3)*(ny-1)*nx+k)+6.0*t((nz-4)*(ny-1)*nx+k)-
     +     t((nz-5)*(ny-1)*nx+k)
          enddo
!$OMP END PARALLEL DO
        else
!$OMP PARALLEL DO 
          do k=1,nx*(ny-1)
            tnew3((i-1)*(ny-1)*nx+k)=-24.0d0*(t((i-2)*(ny-1)*nx+k)-
     +      t((i-1)*(ny-1)*nx+k))
          enddo
!$OMP END PARALLEL DO
        endif
      enddo

!==================================================================== 
! B6*t=tnew1 tnew1: nx*(ny-1)*(nz-1) t: nx*(ny-1)*nz
!==================================================================== 
      elseif (tag.eq.6) then
!$OMP PARALLEL DO COLLAPSE(3)
       do i=1,nz-1
        do j=1,ny-1
         do k=1,nx
          tnew1(k+(j-1)*nx+(i-1)*(ny-1)*nx)=-t(k+(j-1)*nx+(i-1)*(ny-1)
     +    *nx)+t(k+(j-1)*nx+(i-1)*(ny-1)*nx+nx*(ny-1))
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B7*t=tnew3 tnew3: nx*(ny-1)*nz t: (nx-1)*(ny-1)*nz
!==================================================================== 
      elseif (tag.eq.7) then
!$OMP PARALLEL DO 
       do i=1,nz*(ny-1)
        tnew3(1+(i-1)*(nx))=-10.0d0*t(1+(i-1)*(nx-1))+
     +  36.0d0*t(2+(i-1)*(nx-1))-6.0d0*t(3+(i-1)*(nx-1))+
     +  t(4+(i-1)*(nx-1))
        do k=2,nx-1
         tnew3(k+(i-1)*nx)=24.0d0*(t(k+(i-1)*(nx-1))-t(k-1+
     +   (i-1)*(nx-1)))
        enddo
        tnew3(nx+(i-1)*(nx))=10.0d0*t(nx-1+(i-1)*(nx-1))-36.0d0*
     +  t(nx-2+(i-1)*(nx-1))+6.0d0*t(nx-3+(i-1)*(nx-1))-t(nx-4+(i-1)*
     +  (nx-1))
       enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B8*t=tnew1 tnew1: nx*(ny-1)*(nz-1) t: nx*(ny-1)*nz
!==================================================================== 
      elseif (tag.eq.8) then
!$OMP PARALLEL DO COLLAPSE(3)
       do i=1,nz-1
        do j=1,ny-1
         do k=1,nx
          tnew1(k+(j-1)*nx+(i-1)*(ny-1)*nx)=-t(k+(j-1)*nx+(i-1)*(ny-1)*
     +    nx)+t(k+(j-1)*nx+(i-1)*(ny-1)*nx+nx*(ny-1))
          enddo
         enddo
        enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B9*t=tnew2 tnew2: nx*ny*(nz-1) t: nx*(ny-1)*(nz-1)
!==================================================================== 
      elseif (tag.eq.9) then
!$OMP PARALLEL DO COLLAPSE(2)
       do i=1,nz-1
        do j=1,ny
         if (j.eq.1) then
          do k=1,nx
           tnew2((i-1)*ny*nx+k)=-10.0d0*t((i-1)*(ny-1)*nx+k)+36.0d0*
     +     t((i-1)*(ny-1)*nx+nx+k)-6.0d0*t((i-1)*(ny-1)*nx+2*nx+k)+
     +     t((i-1)*(ny-1)*nx+3*nx+k)
          enddo
         elseif (j.eq.ny) then
          do k=1,nx
           tnew2((ny-1)*nx+(i-1)*ny*nx+k)=10.0d0*t((ny-2)*nx+
     +     (i-1)*(ny-1)*nx+k)-36.0d0*t((ny-3)*nx+(i-1)*(ny-1)*nx+k)+
     +     6.0d0*t((ny-4)*nx+(i-1)*(ny-1)*nx+k)-
     +     t((ny-5)*nx+(i-1)*(ny-1)*nx+k)
          enddo
         else
          do k=1,nx
           tnew2((j-1)*nx+(i-1)*ny*nx+k)=-24.0d0*(t((j-2)*nx+
     +     (i-1)*(ny-1)*nx+k)-t((j-1)*nx+(i-1)*(ny-1)*nx+k))
          enddo
         endif
        enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B10*t=tnew tnew4: (nx-1)*ny*(nz-1) t: nx*ny*(nz-1)
!==================================================================== 
      elseif (tag.eq.10) then
!$OMP PARALLEL DO COLLAPSE(3)
       do i=1,nz-1
        do j=1,ny
         do k=1,nx-1
          tnew4(k+(j-1)*(nx-1)+(i-1)*ny*(nx-1))=-t(k+(j-1)*nx+(i-1)
     +    *nx*ny) + t(k+1+(j-1)*nx+(i-1)*ny*nx)
         enddo
        enddo
       enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B11*t=tnew tnew2: nx*ny*(nz-1) t: (nx-1)*ny*(nz-1)
!==================================================================== 
	elseif (tag.eq.11) then
       do i=1,(nz-1)*ny
        tnew2(1+(i-1)*(nx))=-10.0d0*t(1+(i-1)*(nx-1))+36.0d0*
     +  t(2+(i-1)*(nx-1))-6.0d0*t(3+(i-1)*(nx-1))+t(4+(i-1)*(nx-1))
        do k=2,nx-1
         tnew2(k+(i-1)*nx) = 24*(t(k+(i-1)*(nx-1))-t(k-1+(i-1)*(nx-1)))
        enddo
        tnew2(nx+(i-1)*(nx)) = 10.0d0*t(nx-1+(i-1)*(nx-1))-36.0d0*
     +  t(nx-2+(i-1)*(nx-1))+6.0d0*t(nx-3+(i-1)*(nx-1))-
     +  t(nx-4+(i-1)*(nx-1))
       enddo

!==================================================================== 
! B12*t=tnew tnew4: (nx-1)*ny*(nz-1) t: nx*ny*(nz-1)
!==================================================================== 
      elseif (tag.eq.12) then
!$OMP PARALLEL DO COLLAPSE(3)
      do i=1,nz-1
       do j=1,ny
        do k=1,nx-1
         tnew4(k+(j-1)*(nx-1)+(i-1)*ny*(nx-1))=
     +   -t(k+(j-1)*nx+(i-1)*nx*ny)+t(k+1+(j-1)*nx+(i-1)*ny*nx)
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B13*t=tnew tnew5:(nx-1)*ny*nz t:(nx-1)*ny*(nz-1)
!==================================================================== 
      elseif (tag.eq.13) then
!$OMP PARALLEL
!$OMP DO
      do k=1,(nx-1)*ny
       tnew5(k)=-10.0*t(k)+36.0d0*t(k+ny*(nx-1))
     +  -6.0d0*t(k+2*ny*(nx-1))+t(k+3*ny*(nx-1))
      enddo
!$OMP END DO
!$OMP DO COLLAPSE(2)
      do i=2,nz-1
       do k=1,(nx-1)*ny
        tnew5(k+(i-1)*(nx-1)*ny)=-24.0d0*(t(k+(i-2)*(nx-1)*ny)-
     +   t(k+(i-1)*(nx-1)*ny))
       enddo
      enddo
!$OMP END DO
!$OMP DO
      do k=1,(nx-1)*ny
       tnew5(k+(nz-1)*(nx-1)*ny) = 10.0d0*t(k+(nz-2)*(nx-1)*ny)-36.0d0*
     +  t(k+(nz-3)*ny*(nx-1))+6.0d0*t(k+(nz-4)*ny*(nx-1))-
     +  t(k+(nz-5)*ny*(nx-1))
      enddo
!$OMP END DO
!$OMP END PARALLEL

!==================================================================== 
! B14*t=tnew tnew4: (nx-1)*ny*(nz-1) t: (nx-1)*ny*nz
!==================================================================== 
      elseif (tag.eq.14) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz-1
       do k=1,ny*(nx-1)
        tnew4(k+(i-1)*(nx-1)*ny)=
     +     -t(k+(i-1)*(nx-1)*ny)+t(k+i*(nx-1)*ny)
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B15*t=tnew tnew5: (nx-1)*ny*nz t: (nx-1)*(ny-1)*nz
!==================================================================== 
      elseif (tag.eq.15) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz
       do j=1,ny
        if (j.eq.1) then
         do k=1,nx-1
          tnew5(k+(i-1)*ny*(nx-1))=-10.0d0*t(k+(i-1)*(nx-1)*(ny-1))+
     +     36.0d0*t(k+(i-1)*(nx-1)*(ny-1)+(nx-1))-
     +     6.0d0*t(k+(i-1)*(nx-1)*(ny-1)+2*(nx-1))+
     +     t((k+(i-1)*(nx-1)*(ny-1)+3*(nx-1)))
         enddo
        elseif (j.eq.ny) then
         do k=1,nx-1
          tnew5(k+(ny-1)*(nx-1)+(i-1)*(nx-1)*ny)=
     +     10.0d0*t(k+(ny-2)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +     36.0d0*t(k+(ny-3)*(nx-1)+(i-1)*(nx-1)*(ny-1))+
     +     6.0d0*t(k+(ny-4)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +     t((k+(ny-5)*(nx-1)+(i-1)*(nx-1)*(ny-1)))
         enddo
        else
         do k=1,nx-1
          tnew5(k+(j-1)*(nx-1)+(i-1)*(nx-1)*ny)=-24.0d0*
     +     (t(k+(j-2)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +      t(k+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1)))
         enddo
        endif
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B16*t=tnew tnew:(nx-1)*ny*(nz-1) t: (nx-1)*ny*(nz-1)
!==================================================================== 
	elseif (tag.eq.16) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz-1
       do k=1,(nx-1)*ny
        tnew4(k+(i-1)*(nx-1)*ny)=-t(k+(i-1)*(nx-1)*ny)+t(k+i*(nx-1)*ny)
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B17*t=tnew tnew3: nx*(ny-1)*nz t: (nx-1)*(ny-1)*nz
!==================================================================== 
	elseif (tag.eq.17) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz
       do j=1,ny-1
        tnew3(1+(j-1)*nx+(i-1)*(ny-1)*nx)=
     +  -10.0d0*t(1+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))+
     +  36.0d0*t(2+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +  6.0d0*t(3+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))+
     +  t(4+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))
        do k=2,nx-1
         tnew3(k+(j-1)*nx+(i-1)*(ny-1)*nx)=24.0d0*(t(k+(j-1)*(nx-1)+
     +   (i-1)*(nx-1)*(ny-1))-t(k-1+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1)))
        enddo
        tnew3(nx+(j-1)*nx+(i-1)*(ny-1)*nx)=10.0d0*t(nx-1+(j-1)*(nx-1)+
     +  (i-1)*(nx-1)*(ny-1))-36.0d0*t(nx-2+(j-1)*(nx-1)+(i-1)*(nx-1)*
     +  (ny-1))+6.0d0*t(nx-3+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))-t(nx-4+
     +  (j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1))
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B18*t=tnew tnew6: (nx-1)*(ny-1)*nz t: nx*(ny-1)*nz
!==================================================================== 
	elseif (tag.eq.18) then
!$OMP PARALLEL DO COLLAPSE(3)
      do i=1,nz
       do j=1,ny-1
        do k=1,nx-1
         tnew6(k+(j-1)*(nx-1)+(i-1)*(ny-1)*(nx-1)) = -t(k+(j-1)*nx
     +   + (i-1)*nx*(ny-1)) + t(k+1+(j-1)*nx+(i-1)*(ny-1)*nx)
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B19*t=tnew tnew5: (nx-1)*ny*nz t: (nx-1)*(ny-1)*nz
!==================================================================== 
	elseif (tag.eq.19) then
!$OMP PARALLEL DO COLLAPSE(2)
      do i=1,nz
       do j=1,ny
        if (j.eq.1) then
         do k=1,nx-1
          tnew5(k+(i-1)*ny*(nx-1))=-10.0d0*t(k+(i-1)*(nx-1)*(ny-1))+
     +    36.0d0*t(k+(i-1)*(nx-1)*(ny-1)+(nx-1))-
     +    6.0d0*t(k+(i-1)*(nx-1)*(ny-1)+2*(nx-1))+
     +    t((k+(i-1)*(nx-1)*(ny-1) + 3*(nx-1)))
         enddo
        elseif (j.eq.ny) then
         do k=1,nx-1
          tnew5(k+(ny-1)*(nx-1)+(i-1)*(nx-1)*ny)=
     +    10.0d0*t(k+(ny-2)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +    36.0d0*t(k+(ny-3)*(nx-1)+(i-1)*(nx-1)*(ny-1))+
     +    6.0d0*t(k+(ny-4)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +    t((k+(ny-5)*(nx-1)+(i-1)*(nx-1)*(ny-1)))
         enddo
        else
         do k=1,nx-1
          tnew5(k+(j-1)*(nx-1)+(i-1)*(nx-1)*ny)=
     +    -24.0d0*(t(k+(j-2)*(nx-1)+(i-1)*(nx-1)*(ny-1))-
     +    t(k+(j-1)*(nx-1)+(i-1)*(nx-1)*(ny-1)))
         enddo
        endif
       enddo
      enddo
!$OMP END PARALLEL DO

!==================================================================== 
! B20*t=tnew tnew6: (nx-1)*(ny-1)*nz t: (nx-1)*ny*nz
!==================================================================== 
	elseif (tag.eq.20) then
!$OMP PARALLEL DO COLLAPSE(3)
      do i=1,nz
       do j=1,ny-1
        do k=1,nx-1
         tnew6(k+(j-1)*(nx-1)+(i-1)*(ny-1)*(nx-1))=-t(k+(j-1)*(nx-1)+
     +   (i-1)*(nx-1)*ny)+t(k+j*(nx-1)+(i-1)*(nx-1)*ny)
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO
      endif
	  
	  
      return
      end subroutine
