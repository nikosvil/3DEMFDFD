! ===================================================================
! Title: S_solve.f 
! Authors: N. D. Vilanakis, E. N. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
! * A 3D frequency-domain electromagnetic solver employing a high order compact finite-difference scheme
! * N. D. Vilanakis, N. Economou, E. N. Mathioudakis, A. Vafidis
! * Computers and Geosciences
!====================================================================
! S_solve.f computes the solution tnew of the linear system Snum*tnew=tnew
! using classic Cyclic Reduction or the Fourier method
! depending on the structure of the Snum matrix
! Each time the subroutine is being called, integer input Snum specifies
! which operation is being performed using the input y array with the 
! proper dimensions and returning the respective tnew array
!====================================================================      
! Input:
! Snum: index which refers to the respective S matrix
! nx: number of elements in x-direction
! ny: number of elements in y-direction
! nz: number of elements in z-direction
! *and one of the following arrays depending on the Snum index*
! tnew1: complex array, dimension nx*(ny-1)*(nz-1)
! tnew2: complex array, dimension nx*ny*(nz-1)
! tnew3: complex array, dimension nx*(ny-1)*nz
! tnew4: complex array, dimension (nx-1)*ny*(nz-1)
! tnew5: complex array, dimension (nx-1)*ny*nz
! tnew6: complex array, dimension(nx-1)*(ny-1)*nz

! Local auxiliary variables:
! RHS1: complex array, dimension (2*nx-1)*ny*(nz-1)
! RHS2: complex array, dimension nx*(2*ny-1)*(nz-1)
! RHS3: complex array, dimension (2*nx-1)*(ny-1)*nz
! RHS4: complex array, dimension (nx-1)*(2*ny-1)*nz
! RHS5: complex array, dimension (nx-1)*ny*(2*nz-1)
! RHS6: complex array, dimension nx*(ny-1)*(2*nz-1)
! y1: complex array, dimension nx*(ny-1)*(nz-1)
! y2: complex array, dimension nx*ny*(nz-1))
! y3: complex array, dimension nx*(ny-1)*nz
! y4: complex array, dimension (nx-1)*ny*(nz-1)
! y5: complex array, dimension (nx-1)*ny*nz
! y6: complex array, dimension (nx-1)*(ny-1)*nz
! MultVal1: real array, dimension 5*int(dlog(dfloat(nx))/dlog(2.0d0))
! MultVal2: real array, dimension 5*int(dlog(dfloat(ny))/dlog(2.0d0))
! MultVal3: real array, dimension 5*int(dlog(dfloat(nz))/dlog(2.0d0))
! a1: real array, dimension nx
! a2: real array, dimension ny
! a3: real array, dimension nz
! b1: real array, dimension nx
! b2: real array, dimension ny
! b3: real array, dimension nz
! c1: real array, dimension nx
! c2: real array, dimension ny
! c3: real array, dimension nz
! sv1: real array, dimension nx*(ny-1)
! sv2: real array, dimension (nx-1)*ny
! sv: real scalar
! Output:
! *One of the following arrays depending on the Snum index*
! tnew1: complex array, dimension nx*(ny-1)*(nz-1)
! tnew2: complex array, dimension nx*ny*(nz-1)
! tnew3: complex array, dimension nx*(ny-1)*nz
! tnew4: complex array, dimension (nx-1)*ny*(nz-1)
! tnew5: complex array, dimension (nx-1)*ny*nz
! tnew6: complex array, dimension(nx-1)*(ny-1)*nz
!==================================================================== 
	subroutine S_solve(Snum,nx,ny,nz,
     +                   tnew1,tnew2,tnew3,tnew4,tnew5,tnew6,
     +                   y1,y2,y3,y4,y5,y6,
     +                   RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +                   MultVal1,MultVal2,MultVal3,
     +                   a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2)


      implicit none
      real*8, parameter :: pi=4.0d0*datan(1.0d0)
      integer:: nx,ny,nz,m,i,j,k,ii,iz,Snum,lgx,lgy,lgz
      real*8:: s,t,pdn,d,q,dsqrt,dcos,dlog,dsin,dznrm2,jv1,jv2,
     + a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     + MultVal1(5*int(dlog(dfloat(nx))/dlog(2.0d0))),
     + MultVal2(5*int(dlog(dfloat(ny))/dlog(2.0d0))),
     + MultVal3(5*int(dlog(dfloat(nz))/dlog(2.0d0)))
       complex*16::
     + tnew1(nx*(ny-1)*(nz-1)),y1(nx*(ny-1)*(nz-1)),
     + tnew2(nx*ny*(nz-1)),y2(nx*ny*(nz-1)),
     + tnew3(nx*(ny-1)*nz),y3(nx*(ny-1)*nz),
     + tnew4((nx-1)*ny*(nz-1)),y4((nx-1)*ny*(nz-1)),
     + tnew5((nx-1)*ny*nz),y5((nx-1)*ny*nz),
     + tnew6((nx-1)*(ny-1)*nz),y6((nx-1)*(ny-1)*nz), 
     + RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     + RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     + RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     + sv1(nx*(ny-1)),sv2((nx-1)*ny),sv

!==================================================================== 
! Solves S1*tnew2=tnew2 using classic CR
!==================================================================== 
      if (Snum.eq.1) then
!$omp parallel do
       do i=1,size(RHS1)
        RHS1(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,ny*(nz-1)
         do j=1,nx
           RHS1(j+(i-1)*(2*nx-1))=tnew2(j+(i-1)*nx)
         enddo
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(tnew2)
         tnew2(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(MultVal1)
         MultVal1(i)=0.0d0
       enddo
!$omp end parallel do

	k=int(dlog(dfloat(nx))/dlog(2.0d0))

       do ii=1,ny*(nz-1)
        do i=1,nx
         a1(i)=22.0d0
         b1(i)=1.0d0
         c1(i)=1.0d0
        enddo
        a1(1)=24.0d0
        a1(nx)=24.0d0
        b1(1)=24.0d0
        b1(nx)=0.0d0
        c1(1)=0.0d0
        c1(nx)=24.0d0
        do j=1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(nx/(2**j))
         do i=1,m
          s=c1(2*i)/a1(2*i-1)
          if (i.eq.m) then
           t=0.0d0
           MultVal1(5*(j-1)+3)=a1(2*i)-s*b1(2*i-1)
          else
           t=b1(2*i)/a1(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal1(5*(j-1)+i)=a1(2*i)-s*b1(2*i-1)- t*c1(2*i+1)
           endif
          endif

          if (i.eq.1) then
           MultVal1(5*(j-1) + 4) = -t*b1(2*i+1)
           a1(i)=MultVal1(5*(j-1) + 1)
           b1(i)=MultVal1(5*(j-1) + 4)
           c1(i)=0.0d0
          elseif (i.eq.m) then
           MultVal1(5*(j-1) + 5) = -s*c1(2*i-1)
           a1(i)=MultVal1(5*(j-1) + 3)
           b1(i)=0.0d0
           c1(i)=MultVal1(5*(j-1) + 5)
          else
           a1(i)=MultVal1(5*(j-1) + 2)
           b1(i)=MultVal1(5*(j-1) + 4)
           c1(i)=MultVal1(5*(j-1) + 4)
          endif
c MultVal Complete
          if (i.eq.m) then
           RHS1(int(jv1*nx+(nx/2**j)+(ii-1)*(2*nx-1))) =
     +     RHS1(int(jv1*nx+(ii-1)*(2*nx-1))) - 
     +     s*RHS1(int(jv1*nx-1+(ii-1)*(2*nx-1)))
          else
           RHS1(int(jv1*nx+i+(ii-1)*(2*nx-1))) = 
     +     RHS1(int(jv2*nx+2*i+(ii-1)*(2*nx-1))) -
     +     s*RHS1(int(jv2*nx+2*i-1+(ii-1)*(2*nx-1))) -
     +     t*RHS1(int(jv2*nx+2*i+1+(ii-1)*(2*nx-1)))
          endif
         enddo

         MultVal1(5*(k-1) + 1) = MultVal1(5*(k-1) + 3)
         MultVal1(5*(k-1) + 3) = 0.0d0
        enddo

        tnew2(nx+(ii-1)*nx)=RHS1(2*nx-1+(ii-1)*
     +  (2*nx-1))/MultVal1(5*(k-1)+1)
       enddo

!$omp parallel do shared(tnew2,RHS1) private(iz,j,i,jv2,ii)
       do iz=1,ny*(nz-1)
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do  i=2**(j-1),nx-2**(j-1),2**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew2(i+(iz-1)*nx)=(RHS1(1+(iz-1)*(2*nx-1))-
     +     24.0d0*tnew2(2+(iz-1)*nx))/24.0d0
          elseif (i.eq.2**(j-1)) then
           tnew2(i+(iz-1)*nx)=(RHS1(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew2(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+1)
          elseif ((i.eq.nx-2**(j-1)).and.(j.eq.2)) then
           tnew2(i+(iz-1)*nx)=(RHS1(3*nx/2-1+(iz-1)*(2*nx-1))-
     +     MultVal1(4)*tnew2(i-2+(iz-1)*nx)-MultVal1(4)*
     +     tnew2(i+2+(iz-1)*nx))/MultVal1(2)
          elseif (j.eq.1) then
           tnew2(i+(iz-1)*nx)=(RHS1(i+(iz-1)*(2*nx-1))-
     +     tnew2(i-1+(iz-1)*nx)-tnew2(i+1+(iz-1)*nx))/22.0d0
          else
           tnew2(i+(iz-1)*nx)=(RHS1(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew2(i-2**(j-1)+(iz-1)*nx)-
     +     MultVal1(5*(j-2)+4)*tnew2(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+2)
          endif
         enddo
        enddo
       enddo
!$omp end parallel do

!==================================================================== 
! Solves S2*tnew1=tnew1 using Fourier method
!==================================================================== 
      elseif (Snum.eq.2) then
 
       call zcopy(size(tnew1),tnew1,1,y1,1)
!$omp parallel do
       do i=1,size(tnew1)
        tnew1(i)=0.0d0
       enddo
!$omp end parallel do
       pdn=pi/dfloat(nz)
       q=dsqrt(2.0d0/nz)

       do i=1,nz-1
!$omp parallel do
        do k=1,nx*(ny-1)
          sv1(k)=dcmplx(0.0d0,0.0d0)
        enddo
!$omp end parallel do
!$omp critical
        do j=1,nz-1
         sv1(1:nx*(ny-1))=
     +   sv1(1:nx*(ny-1))+y1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))*
     +   q*dsin(i*j*pdn)
        enddo
        sv1(1:nx*(ny-1))=
     +  sv1(1:nx*(ny-1))*(1.0d0/(22.0d0+2.0d0*dcos(i*pdn)))
!$omp end critical
!$omp parallel do shared(tnew1)
        do j=1,nz-1
         tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))=
     +   tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))+
     +   sv1(1:nx*(ny-1))*q*dsin(i*j*pdn)
        enddo
!$omp end parallel do
       enddo

!==================================================================== 
! Solves S3*tnew2=tnew2 using classic CR
!==================================================================== 
	elseif (Snum.eq.3) then
!$omp parallel do
       do i=1,size(RHS2)
         RHS2(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,nz-1
          RHS2((i-1)*(2*ny-1)*nx+1:(i-1)*(2*ny-1)*nx+nx*ny)=
     +    tnew2((i-1)*nx*ny+1:i*nx*ny)
       enddo
!$omp end parallel do

!$omp parallel do
	 do i=1,size(tnew2)
         tnew2(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(MultVal2)
         MultVal2(i)=0.0d0
       enddo
!$omp end parallel do

       k=int(dlog(dfloat(ny))/dlog(2.0d0))


       do iz=1,nz-1
        do i = 1,ny
         a2(i) = 22.0d0
         b2(i) = 1.0d0
         c2(i) = 1.0d0
        enddo
        a2(1) = 24.0d0
        a2(ny) = 24.0d0
        b2(ny) = 0.0d0
        c2(1) = 0.0d0
        c2(ny) = 24.0d0
        b2(1) = 24.0d0

        do j = 1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(ny/2**j)
         do i = 1,m
          s = c2(2*i)/a2(2*i-1)
          if (i.eq.m) then
           t = 0.0d0
           MultVal2((j-1)*5 + 3)=a2(2*i)-s*b2(2*i-1)
          else
           t = b2(2*i)/a2(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal2((j-1)*5 + i) = a2(2*i)-s*b2(2*i-1)-t*c2(2*i+1)
           endif
          endif
            
          if (i.eq.1) then
           MultVal2((j-1)*5 + 4) = -t*b2(2*i+1)
           a2(i)=MultVal2((j-1)*5 + 1)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=0.0d0
          elseif (i.eq.m) then
           MultVal2((j-1)*5 + 5) = -s*c2(2*i-1)
           a2(i)=MultVal2((j-1)*5 + 3)
           b2(i)=0.0d0
           c2(i)=MultVal2((j-1)*5 + 5)
          else
           a2(i)=MultVal2((j-1)*5 + 2)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=MultVal2((j-1)*5 + 4)
          endif
c  Multval complete
          if (i.eq.m) then
           RHS2(int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)-1)*nx+1):
     +          int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)+0)*nx)) = 
     +     RHS2(int((jv1*ny+(iz-1)*(2*ny-1)-1)*nx+1):
     +          int((jv1*ny+(iz-1)*(2*ny-1)-1)*nx+nx)) - s*
     +     RHS2(int((jv1*ny+(iz-1)*(2*ny-1)-2)*nx+1):
     +          int((jv1*ny+(iz-1)*(2*ny-1)-2)*nx+nx))
          else
           RHS2(int((jv1*ny+i+(iz-1)*(2*ny-1)-1)*nx+1):
     +          int((jv1*ny+i+(iz-1)*(2*ny-1))*nx)) = 
     +     RHS2(int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*nx+1):
     +          int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*nx+nx))-s*
     +     RHS2(int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*nx+1):
     +          int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*nx+nx))-t*
     +     RHS2(int((jv2*ny+2*i+(iz-1)*(2*ny-1))*nx+1):
     +          int((jv2*ny+2*i+(iz-1)*(2*ny-1))*nx+nx))            
          endif
         enddo
         MultVal2((k-1)*5 + 1) = MultVal2((k-1)*5 + 3)
         MultVal2((k-1)*5 + 3) = 0.0d0
        enddo
        tnew2(((ny-1)+(iz-1)*ny)*nx+1:((ny-1)+(iz-1)*ny)*nx+nx)=
     +  RHS2(nx*(2*ny-2+(iz-1)*(2*ny-1))+1:nx*(2*ny-1+(iz-1)*
     +  (2*ny-1)))*1.0d0/MultVal2(5*(k-1)+1)
       enddo

!$omp parallel do shared(tnew2,RHS2) private(iz,j,i,jv2,ii)
       do iz=1,nz-1
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/2.0d0**(j-2)
         ii=0
         do i=2.0d0**(j-1),ny-2.0d0**(j-1),2.0d0**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew2((nx)*(iz-1)*ny+1:(nx)*(iz-1)*ny+(nx))=
     +     (RHS2(1+(nx)*(iz-1)*(2*ny-1):nx+(nx)*(iz-1)*(2*ny-1))-24.0d0*
     +     tnew2(1+(nx)+(nx)*(iz-1)*ny:2*(nx)+(nx)*(iz-1)*ny))*
     +     (1.0d0/24.0d0)
          elseif (i.eq.2**(j-1)) then
           tnew2((i-1)*(nx)+(iz-1)*ny*(nx)+1:
     +           (i-1)*(nx)+(iz-1)*ny*(nx)+(nx))=
     +            (1.0d0/MultVal2(5*(j-2)+1))*
     +     (RHS2(int((nx)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+1):
     +           int((nx)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+(nx)))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew2(((i-1)+2**(j-1)+(iz-1)*ny)*(nx)+1:
     +            ((i-1)+2**(j-1)+(iz-1)*ny)*(nx)+(nx)))    
          elseif ((i.eq.(ny-2**(j-1))).and.(j.eq.2)) then
           tnew2(((i-1)+(iz-1)*ny)*(nx)+1:((i-1)+(iz-1)*ny)*(nx)+(nx))=
     +     (RHS2((3*ny/2-2+(iz-1)*(2*ny-1))*(nx)+1:
     +           (3*ny/2-2+(iz-1)*(2*ny-1))*(nx)+(nx))-
     +     MultVal2(4)*tnew2((i-3+(iz-1)*ny)*(nx)+1:
     +                       (i-3+(iz-1)*ny)*(nx)+(nx))-
     +     MultVal2(4)*tnew2((i+1+(iz-1)*ny)*(nx)+1:
     +                       (i+1+(iz-1)*ny)*(nx)+(nx)))*
     +     (1.0d0/MultVal2(2))     
          elseif (j.eq.1) then
           tnew2(((i-1)+(iz-1)*ny)*(nx)+1:(i-1+(iz-1)*ny)*(nx)+(nx))=
     +     (RHS2((i-1+(iz-1)*(2*ny-1))*(nx)+1:
     +           (i-1+(iz-1)*(2*ny-1))*(nx)+(nx))-
     +      tnew2((i-2+(iz-1)*ny)*(nx)+1:
     +            (i-2+(iz-1)*ny)*(nx)+(nx))-
     +      tnew2((i+(iz-1)*ny)*(nx)+1:
     +            (i+(iz-1)*ny)*(nx)+(nx)))*(1.0d0/22.0d0)           
          else
           tnew2((i-1+(iz-1)*ny)*(nx)+1:(i-1+(iz-1)*ny)*(nx)+(nx))=
     +     (RHS2(int((ny*jv2+(iz-1)*(2*ny-1))*(nx)+(nx)*(2*ii-2)+1):
     +           int((ny*jv2+(iz-1)*(2*ny-1))*(nx)+(nx)*(2*ii-2)+nx))-
     +      MultVal2(5*(j-2)+4)*
     +     tnew2((i-1-2**(j-1)+(iz-1)*ny)*(nx)+1:
     +           (i-1-2**(j-1)+(iz-1)*ny)*(nx)+(nx))-
     +      MultVal2(5*(j-2)+4)*
     +     tnew2((i-1+2**(j-1)+(iz-1)*ny)*(nx)+1:
     +           (i-1+2**(j-1)+(iz-1)*ny)*(nx)+(nx)))*
     +      (1.0d0/MultVal2(5*(j-2)+2))  
          endif
         enddo
        enddo
       enddo
!$omp end parallel do

!==================================================================== 
! Solves S4*tnew1=tnew1 using Fourier method
!==================================================================== 
	elseif (Snum.eq.4) then
	
	call zcopy(size(tnew1),tnew1,1,y1,1)
!$omp parallel do
	do i=1,size(tnew1)
	 tnew1(i)=dcmplx(0.0d0,0.0d0)
	enddo
!$omp end parallel do

	q=dsqrt(2.0d0/dfloat(ny))
    	pdn=pi/dfloat(ny)
!$omp parallel do private(k,i,iz,sv,j)
      do iz=0,nz-2
       do k=1,nx
        do i=1,ny-1
         sv=dcmplx(0.0d0,0.0d0)
         do j=0,ny-2
          sv=sv+y1(j*nx+iz*nx*(ny-1)+k)*q*dsin(i*(j+1)*pdn)
         enddo
         sv=sv*(1.0d0/(22.0d0+2.0d0*dcos(i*pdn)))
         do j=1,ny-1
          tnew1(k+(j-1)*nx+iz*nx*(ny-1))=
     +    tnew1(k+(j-1)*nx+iz*nx*(ny-1))+sv*q*dsin(i*j*pdn)
         enddo
        enddo
       enddo
      enddo
!$omp end parallel do

!==================================================================== 
! Solves S5*tnew3=tnew3 using classic CR 
!==================================================================== 
	elseif (Snum.eq.5) then
!$omp parallel do
       do i=1,size(RHS6)
         RHS6(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(tnew3)
        RHS6(i)=tnew3(i)
       enddo
!$omp end parallel do

!$omp parallel do
      do i=1,size(tnew3)
       tnew3(i)=dcmplx(0.0d0,0.0d0)
      enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(MultVal3)
         MultVal3(i)=0.0d0
       enddo
!$omp end parallel do

	k=int(dlog(dfloat(nz))/dlog(2.0d0))
 
!$omp parallel do
       do i=1,nz
        a3(i)=22.0d0
        b3(i)=1.0d0
        c3(i)=1.0d0
       enddo
!$omp end parallel do
    	 a3(1) = 24.0d0
    	 a3(nz) = 24.0d0
    	 b3(nz) = 0.0d0
    	 c3(1) = 0.0d0
    	 c3(nz) = 24.0d0
    	 b3(1) = 24.0d0
       do j=1,k
        jv1=(2.0d0**j-1)/(2.0d0**(j-1))
        jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
        m = int(nz/2**j)
        do i = 1,m
         s = c3(2*i)/a3(2*i-1)
         if (i.eq.m) then
          t=0.0d0
          s=c3(2*i)/a3(2*i-1)
          MultVal3((j-1)*5 + 3) = a3(2*i)-s*b3(2*i-1)
         else	
          t=b3(2*i)/a3(2*i+1)
          if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
           MultVal3((j-1)*5 + i)=a3(2*i)-s*b3(2*i-1)-t*c3(2*i+1)
          endif
         endif
         if (i.eq.1) then
          MultVal3((j-1)*5 + 4) = -t*b3(2*i+1)
          a3(i)=MultVal3((j-1)*5 + 1)
          b3(i)=MultVal3((j-1)*5 + 4)
          c3(i)=0.0d0
         elseif (i.eq.m) then
          MultVal3((j-1)*5 + 5) = -s*c3(2*i-1)
          a3(i)=MultVal3((j-1)*5 + 3)
          b3(i)=0.0d0
          c3(i)=MultVal3((j-1)*5 + 5)
         else
          a3(i)=MultVal3((j-1)*5 + 2)
          b3(i)=MultVal3((j-1)*5 + 4)
          c3(i)=MultVal3((j-1)*5 + 4)
         endif
c MultVal complete
         if (i.eq.m) then
          RHS6(int((jv1*nz+(nz/2**j)-1)*nx*(ny-1)+1):
     +         int((jv1*nz+(nz/2**j))*nx*(ny-1)))=
     +    RHS6(int((jv1*nz-1)*nx*(ny-1)+1):
     +         int((jv1*nz-1)*nx*(ny-1)+nx*(ny-1)))-
     +    s*RHS6(int((jv1*nz-2)*nx*(ny-1)+1):
     +           int((jv1*nz-2)*nx*(ny-1)+nx*(ny-1)))
         else
          RHS6(int((jv1*nz+i-1)*nx*(ny-1)+1):int((jv1*nz+i)*nx*(ny-1)))=
     +    RHS6(int((jv2*nz-1+2*i)*nx*(ny-1)+1):
     +         int((jv2*nz-1+2*i)*nx*(ny-1)+nx*(ny-1)))-
     +    s*RHS6(int((jv2*nz+2*i-2)*nx*(ny-1)+1):
     +           int((jv2*nz+2*i-2)*nx*(ny-1)+nx*(ny-1)))-
     +    t*RHS6(int((jv2*nz+2*i)*nx*(ny-1)+1):
     +           int((jv2*nz+2*i)*nx*(ny-1)+nx*(ny-1)))        
         endif
        enddo
        MultVal3((k-1)*5 + 1) = MultVal3((k-1)*5 + 3)
        MultVal3((k-1)*5 + 3) = 0.0d0
    	 enddo
       tnew3((nz-1)*nx*(ny-1)+1:nz*nx*(ny-1))=RHS6(2*(nz-1)*nx*(ny-1)+1:
     +      (2*nz-1)*nx*(ny-1))*1.0d0/MultVal3(5*(k-1)+1)
  
   	 do j=k,1,-1
        jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
        ii=0
        do i=2**(j-1),nz-2**(j-1),2**j
         ii=ii+1
         if ((i.eq.1).and.(j.eq.1)) then	
          tnew3(1:nx*(ny-1))=(RHS6(1:nx*(ny-1))-
     +        24.0d0*tnew3(nx*(ny-1)+1:2*nx*(ny-1)))/24.0d0
         elseif (i.eq.(2**(j-1))) then
          tnew3((i-1)*nx*(ny-1)+1:i*nx*(ny-1))=
     +    (RHS6(int((nz*jv2+2*ii-2)*nx*(ny-1)+1):
     +          int((nz*jv2+2*ii-1)*nx*(ny-1)))-
     +    MultVal3(5*(j-2)+4)*tnew3((i-1+2**(j-1))*nx*(ny-1)+1:
     +    (i+2**(j-1))*nx*(ny-1)))/MultVal3(5*(j-2)+1)
         elseif ((i.eq.(nz-2**(j-1))).and.(j.eq.2)) then
          tnew3((i-1)*nx*(ny-1)+1:i*nx*(ny-1))=
     +    (RHS6((3*nz/2-2)*nx*(ny-1)+1:(3*nz/2-1)*nx*(ny-1))-
     +    MultVal3(4)*
     +    tnew3((i-3)*nx*(ny-1)+1:(i-2)*nx*(ny-1))-
     +    MultVal3(4)*
     +    tnew3((i+1)*nx*(ny-1)+1:(i+2)*nx*(ny-1)))/MultVal3(2)
         elseif (j.eq.1) then
          tnew3((i-1)*nx*(ny-1)+1:i*nx*(ny-1))=
     +    (RHS6((i-1)*nx*(ny-1)+1:i*nx*(ny-1))-
     +    tnew3((i-2)*nx*(ny-1)+1:(i-1)*nx*(ny-1))-
     +    tnew3(i*nx*(ny-1)+1:(i+1)*nx*(ny-1)))/22.0d0
         else  
          tnew3((i-1)*nx*(ny-1)+1:i*nx*(ny-1))=
     +    (RHS6(int((nz*jv2+2*ii-2)*nx*(ny-1)+1):
     +          int((nz*jv2+2*ii-1)*nx*(ny-1)))-
     +    MultVal3(5*(j-2)+4)*
     +    tnew3((i-2**(j-1)-1)*nx*(ny-1)+1:(i-2**(j-1))*nx*(ny-1))-
     +    MultVal3(5*(j-2)+4)*
     +    tnew3((i+2**(j-1)-1)*nx*(ny-1)+1:(i+2**(j-1))*nx*(ny-1)))/
     +    MultVal3(5*(j-2)+2)
         endif
       enddo
      enddo

!==================================================================== 
! Solves S6*tnew1=tnew1 using Fourier method
!==================================================================== 
	elseif (Snum.eq.6) then
       
	 call zcopy(size(tnew1),tnew1,1,y1,1)
!$omp parallel do 
      do i=1,size(tnew1)
       tnew1(i)=dcmplx(0.0d0,0.0d0)
      enddo
!$omp end parallel do

    	 pdn=pi/dfloat(nz)
    	 q=dsqrt(2.0d0/dfloat(nz))

        do i=1,nz-1
!$omp parallel do 
	   do iz=1,nx*(ny-1)
            sv1(iz)=dcmplx(0.0d0,0.0d0)
           enddo
!$omp end parallel do
!$omp critical
       do j=1,nz-1
         sv1(1:nx*(ny-1))=sv1(1:nx*(ny-1))+ 
     +   y1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))*q*dsin(i*j*pdn)
       enddo
        sv1(1:nx*(ny-1))=
     +  sv1(1:nx*(ny-1))*(1.0d0/(22.0d0+2.0d0*dcos(i*pdn)))
!$omp end critical

!$omp parallel do shared(tnew1)
        do j=1,nz-1
         tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))=
     +   tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))+
     +   sv1(1:nx*(ny-1))*q*dsin(i*j*pdn)
        enddo
!$omp end parallel do
       enddo
       
!==================================================================== 
! Solves S7*tnew3=tnew3 using classic CR
!==================================================================== 
	elseif (Snum.eq.7) then
!$omp parallel do
    	 do i=1,size(RHS3)
	  RHS3(i)=dcmplx(0.0d0,0.0d0)
 	 enddo       
!$omp end parallel do
!$omp parallel do
      do ii=1,nz*(ny-1)
        do i=1,nx
         RHS3(i+(ii-1)*(2*nx-1))=tnew3(i+(ii-1)*nx)
        enddo
       enddo
!$omp end parallel do
!$omp parallel do
	 do i=1,size(tnew3)
        tnew3(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do
!$omp parallel do
       do i=1,size(MultVal1)
         MultVal1(i)=0.0d0
       enddo
!$omp end parallel do
	 k=int(dlog(dfloat(nx))/dlog(2.0d0))

    	 do ii=1,(ny-1)*nz
        do i = 1,nx
         a1(i) = 22.0d0
         b1(i) = 1.0d0
         c1(i) = 1.0d0
        enddo
        a1(1) = 24.0d0
        a1(nx) = 24.0d0
        b1(nx) = 0.0d0
        c1(1) = 0.0d0
        c1(nx) = 24.0d0
        b1(1) = 24.0d0
        do j = 1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(nx/2**j)
         do i=1,m
          s=c1(2*i)/a1(2*i-1)
          if (i.eq.m) then
           t=0.0d0
           MultVal1((j-1)*5 + 3)=a1(2*i)-s*b1(2*i-1)
          else
           t=b1(2*i)/a1(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal1((j-1)*5 + i)=a1(2*i)-s*b1(2*i-1)-t*c1(2*i+1)
           endif
          endif
            
          if (i.eq.1) then
           MultVal1((j-1)*5 + 4) = -t*b1(2*i+1)
           a1(i)=MultVal1((j-1)*5 + 1)
           b1(i)=MultVal1((j-1)*5 + 4)
           c1(i)=0.0d0
          elseif (i.eq.m) then
           MultVal1((j-1)*5 + 5) = -s*c1(2*i-1)
           a1(i)=MultVal1((j-1)*5 + 3)
           b1(i)=0.0d0
           c1(i)=MultVal1((j-1)*5 + 5)
          else
           a1(i)=MultVal1((j-1)*5 + 2)
           b1(i)=MultVal1((j-1)*5 + 4)
           c1(i)=MultVal1((j-1)*5 + 4)
          endif
c MultVal complete
          if (i.eq.m) then
           RHS3(int(jv1*nx+(nx/2**j))+(ii-1)*(2*nx-1))=
     +     RHS3(int(jv1*nx+(ii-1)*(2*nx-1)))-
     +     s*RHS3(int(jv1*nx-1)+(ii-1)*(2*nx-1))
          else 
           RHS3(int(jv1*nx+i)+(ii-1)*(2*nx-1))=
     +     RHS3(int(jv2*nx)+2*i+(ii-1)*(2*nx-1))-
     +     s*RHS3(int(jv2*nx)+2*i-1+(ii-1)*(2*nx-1))-
     +     t*RHS3(int(jv2*nx)+2*i+1+(ii-1)*(2*nx-1))             
          endif
         enddo
         MultVal1((k-1)*5 + 1) = MultVal1((k-1)*5 + 3)
         MultVal1((k-1)*5 + 3) = 0.0d0
        enddo

        tnew3(nx+(ii-1)*nx)=RHS3(2*nx-1+(ii-1)*(2*nx-1))/
     +  MultVal1(5*(k-1)+1)
       enddo


!$omp parallel do shared(tnew3,RHS3) private(iz,j,i,jv2,ii)
       do iz=1,(ny-1)*nz
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do i=2**(j-1),nx-2**(j-1),2**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew3(i+(iz-1)*nx)=
     +     (RHS3(1+(iz-1)*(2*nx-1))-24.0d0*tnew3(2+(iz-1)*nx))/24.0d0
          elseif (i.eq.(2**(j-1))) then
	     tnew3(i+(iz-1)*nx)=(RHS3(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew3(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+1)
          elseif ((i.eq.nx-2**(j-1)).and.(j.eq.2)) then
           tnew3(i+(iz-1)*nx)=(RHS3(3*nx/2-1+(iz-1)*(2*nx-1))-
     +     MultVal1(4)*tnew3(i-2+(iz-1)*nx)-
     +     MultVal1(4)*tnew3(i+2+(iz-1)*nx))/MultVal1(2)
          elseif (j.eq.1) then
           tnew3(i+(iz-1)*nx)=(RHS3(i+(iz-1)*(2*nx-1))-
     +     tnew3(i-1+(iz-1)*nx)-tnew3(i+1+(iz-1)*nx))/22.0d0
          else
           tnew3(i+(iz-1)*nx)=(RHS3(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew3(i-2**(j-1)+(iz-1)*nx)-
     +     MultVal1(5*(j-2)+4)*tnew3(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+2)
          endif
         enddo
        enddo
       enddo
!$omp end parallel do 

!==================================================================== 
! Solves S8*tnew1=tnew1 using Fourier method 
!==================================================================== 
	elseif (Snum.eq.8) then

	 call zcopy(size(tnew1),tnew1,1,y1,1)
!$omp parallel do    
	 do i=1,size(tnew1)
	  tnew1(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 
       pdn=pi/dfloat(nz)
       q=dsqrt(2.0d0/dfloat(nz))
       do iz=1,nz-1
!$omp parallel do    
        do i=1,nx*(ny-1)
         sv1(i)=dcmplx(0.0d0,0.0d0)
	enddo
!$omp end parallel do
!$omp critical
        do j=1,nz-1
         sv1(1:nx*(ny-1))=sv1(1:nx*(ny-1))+
     +   y1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))*q*dsin(iz*j*pdn)
        enddo
        sv1(1:nx*(ny-1))=
     +  sv1(1:nx*(ny-1))*(1.0d0/(22.0d0+2.0d0*dcos(iz*pdn)))
!$omp end critical
!$omp parallel do shared(tnew1)
        do j=1,nz-1
         tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))=
     +   tnew1((j-1)*nx*(ny-1)+1:j*nx*(ny-1))+
     +   sv1(1:nx*(ny-1))*q*dsin(iz*j*pdn)
        enddo
!$omp end parallel do 
       enddo

!==================================================================== 
! Solves S9*tnew2=tnew2 using classic CR
!==================================================================== 
	elseif (Snum.eq.9) then
!$omp parallel do 
       do i=1,size(RHS2)
        RHS2(i)=dcmplx(0.0d0,0.0d0)
	   enddo
!$omp end parallel do
!$omp parallel do
   	 do i=1,nz-1
        RHS2((i-1)*(2*ny-1)*nx+1:(i-1)*(2*ny-1)*nx+nx*ny)=
     +   tnew2((i-1)*nx*ny+1:i*nx*ny)
       enddo
!$omp end parallel do
!$omp parallel do  
       do i=1,size(tnew2)
	  tnew2(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do  
    	 do i=1,size(MultVal2)
        MultVal2(i)=0.0d0
	 enddo 
!$omp end parallel do 
	 k=int(dlog(dfloat(ny))/dlog(2.0d0))
c       k=lgy

       do iz=1,nz-1
        do i = 1,ny
         a2(i) = 22.0d0
         b2(i) = 1.0d0
         c2(i) = 1.0d0
        enddo
        a2(1) = 24.0d0
        a2(ny) = 24.0d0
        b2(1) = 24.0d0
        b2(ny) = 0.0d0
        c2(1) = 0.0d0
        c2(ny) = 24.0d0

        do j = 1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(ny/2**j)
         do i=1,m
          s=c2(2*i)/a2(2*i-1)
          if (i.eq.m) then
           t=0.0d0 
           MultVal2((j-1)*5 + 3) = a2(2*i)-s*b2(2*i-1)
          else
           t=b2(2*i)/a2(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal2((j-1)*5 + i) = a2(2*i)-s*b2(2*i-1)-t*c2(2*i+1)
           endif
          endif
            
          if (i.eq.1) then
          MultVal2((j-1)*5 + 4) = -t*b2(2*i+1)
          a2(i)=MultVal2((j-1)*5 + 1)
          b2(i)=MultVal2((j-1)*5 + 4)
          c2(i)=0.0d0
         elseif (i.eq.m) then
          MultVal2((j-1)*5 + 5) = -s*c2(2*i-1)
          a2(i)=MultVal2((j-1)*5 + 3)
          b2(i)=0.0d0
          c2(i)=MultVal2((j-1)*5 + 5)
         else
          a2(i)=MultVal2((j-1)*5 + 2)
          b2(i)=MultVal2((j-1)*5 + 4)
          c2(i)=MultVal2((j-1)*5 + 4)
         endif
c MultVal complete
         if (i.eq.m) then
          RHS2(int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)-1)*nx+1):
     +    int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)+0)*nx))=
     +    RHS2(int((jv1*ny+(iz-1)*(2*ny-1)-1)*nx+1):
     +    int((jv1*ny+(iz-1)*(2*ny-1)-1)*nx+nx))-
     +    s*RHS2(int((jv1*ny+(iz-1)*(2*ny-1)-2)*nx+1):
     +    int((jv1*ny+(iz-1)*(2*ny-1)-2)*nx+nx))
         else
          RHS2(int((jv1*ny+i+(iz-1)*(2*ny-1)-1)*nx+1):
     +    int((jv1*ny+i+(iz-1)*(2*ny-1))*nx))=
     +    RHS2(int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*nx+1):
     +    int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*nx+nx))-
     +    s*RHS2(int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*nx+1):
     +    int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*nx+nx))-
     +    t*RHS2(int((jv2*ny+2*i+(iz-1)*(2*ny-1))*nx+1):
     +    int((jv2*ny+2*i+(iz-1)*(2*ny-1))*nx+nx))             
          endif
         enddo
         MultVal2((k-1)*5 + 1) = MultVal2((k-1)*5 + 3)
         MultVal2((k-1)*5 + 3) = 0.0d0
        enddo
        tnew2(((ny-1)+(iz-1)*ny)*nx+1:((ny-1)+(iz-1)*ny)*nx+nx)=
     +  RHS2(nx*(2*ny-2+(iz-1)*(2*ny-1))+1:
     +       nx*(2*ny-1+(iz-1)*(2*ny-1)))*1.0d0/MultVal2(5*(k-1)+1)
       enddo
!$omp parallel do shared(tnew2,RHS2) private(iz,i,j,jv2,ii)
       do iz=1,nz-1
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do i=2**(j-1),ny-2**(j-1),2**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew2(nx*(iz-1)*ny+1:nx*(iz-1)*ny+nx)=
     +     (RHS2(1+nx*(iz-1)*(2*ny-1):nx+nx*(iz-1)*(2*ny-1))-
     +     24.0d0*tnew2(1+nx+nx*(iz-1)*ny:2*nx+nx*(iz-1)*ny))*1.0d0/
     +     24.0d0
          elseif (i.eq.(2**(j-1))) then
           tnew2((i-1)*nx+(iz-1)*ny*nx+1:(i-1)*nx+(iz-1)*ny*nx+nx)=
     +     1.0d0/MultVal2(5*(j-2)+1)*
     +     (RHS2(int(nx*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+1):
     +           int(nx*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+nx))-
     +     MultVal2(5*(j-2)+4)*tnew2(((i-1)+2**(j-1)+(iz-1)*ny)*nx+1:
     +     ((i-1)+2**(j-1)+(iz-1)*ny)*nx+nx))
          elseif ((i.eq.(ny-2**(j-1))).and.(j.eq.2)) then
           tnew2(((i-1)+(iz-1)*ny)*nx+1:((i-1)+(iz-1)*ny)*nx+nx)=
     +     (RHS2((3*ny/2-1+(iz-1)*(2*ny-1))*nx+1:
     +     (3*ny/2-1+(iz-1)*(2*ny-1))*nx+nx)-MultVal2(4)*
     +     tnew2((i-3+(iz-1)*ny)*nx+1:(i-3+(iz-1)*ny)*nx+nx)-
     +     MultVal2(4)*tnew2((i+1+(iz-1)*ny)*nx+1:
     +     (i+1+(iz-1)*ny)*nx+nx))*1.0d0/MultVal2(2)
          elseif (j.eq.1) then
           tnew2(((i-1)+(iz-1)*ny)*nx+1:(i-1+(iz-1)*ny)*nx+nx)=
     +     (RHS2((i-1+(iz-1)*(2*ny-1))*nx+1:
     +     (i-1+(iz-1)*(2*ny-1))*nx+nx)-tnew2((i-2+(iz-1)*ny)*nx+1:
     +     (i-2+(iz-1)*ny)*nx+nx)-tnew2((i+(iz-1)*ny)*nx+1:
     +     (i+(iz-1)*ny)*nx+nx))*1.0d0/22.0d0
          else
           tnew2((i-1+(iz-1)*ny)*nx+1:(i-1+(iz-1)*ny)*nx+nx)=
     +     (RHS2(int((ny*jv2+(iz-1)*(2*ny-1))*nx+nx*(2*ii-2)+1):
     +           int((ny*jv2+(iz-1)*(2*ny-1))*nx+nx*(2*ii-2)+nx))-
     +     MultVal2(5*(j-2)+4)*tnew2((i-1-2**(j-1)+(iz-1)*ny)*nx+1:
     +     (i-1-2**(j-1)+(iz-1)*ny)*nx+nx)-
     +     MultVal2(5*(j-2)+4)*tnew2((i-1+2**(j-1)+(iz-1)*ny)*nx+1:
     +     (i-1+2**(j-1)+(iz-1)*ny)*nx+nx))*1.0d0/MultVal2(5*(j-2)+2)
          endif
         enddo
        enddo
    	 enddo
!$omp end parallel do

!==================================================================== 
! Solves S10*tnew4=tnew4 using Fourier method
!==================================================================== 
	elseif (Snum.eq.10) then
        
	 call zcopy(size(tnew4),tnew4,1,y4,1)
!$omp parallel do 
	 do i=1,size(tnew4)
	  tnew4(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 

    	 pdn=pi/dfloat(nx)
    	 q=dsqrt(2.0d0/dfloat(nx))

!$omp parallel do shared(tnew4) private (k,i,sv,j)
       do k=1,ny*(nz-1)
        do i=1,nx-1
         sv=dcmplx(0.0d0,0.0d0)
         do j=1,nx-1
          sv=sv+y4((k-1)*(nx-1)+j)*q*dsin(i*j*pdn)
         enddo
         sv=sv*(1.0d0/(22.0d0+2.0d0*dcos(i*pdn)))
         do j=1,nx-1
          tnew4((k-1)*(nx-1)+j)=tnew4((k-1)*(nx-1)+j)+sv*q*dsin(i*j*pdn)
         enddo
        enddo
       enddo
!$omp end parallel do 

!==================================================================== 
! Solves S11*tnew2=tnew2 using classic CR
!==================================================================== 
	elseif (Snum.eq.11) then
!$omp parallel do 
	 do i=1,size(RHS1)
	  RHS1(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do 
    	 do ii=1,ny*(nz-1)
        do i=1,nx
         RHS1(i+(ii-1)*(2*nx-1))=tnew2(i+(ii-1)*nx)
        enddo
    	 enddo
!$omp end parallel do 
!$omp parallel do 
       do i=1,size(tnew2)
	  tnew2(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do 
	 do i=1,size(MultVal1)
	  MultVal1(i)=0.0d0
	 enddo
!$omp end parallel do 
	 k=int(dlog(dfloat(nx))/dlog(2.0d0))

c RHS 
    	 do iz=1,(nz-1)*ny
        do i=1,nx
         a1(i) = 22.0d0
         b1(i) = 1.0d0
         c1(i) = 1.0d0
        enddo
        a1(1)=24.0d0
        a1(nx)=24.0d0
        b1(nx)=0.0d0
        c1(1) = 0.0d0
        c1(nx)=24.0d0
        b1(1)=24.0d0
        do j=1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(nx/2**j)
c	   print*,'m',m
         do i=1,m
          s=c1(2*i)/a1(2*i-1)
          if (i.eq.m) then
           t=0.0d0
           MultVal1((j-1)*5 + 3) = a1(2*i)-s*b1(2*i-1)
          else
           t=b1(2*i)/a1(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal1((j-1)*5 + i) = a1(2*i) - s*b1(2*i-1) - t*c1(2*i+1)
           endif
          endif
         
          if (i.eq.1) then
           MultVal1((j-1)*5 + 4) = -t*b1(2*i+1)
           a1(i)=MultVal1((j-1)*5 + 1)
           b1(i)=MultVal1((j-1)*5 + 4)
           c1(i)=0.0d0
          elseif (i.eq.m) then
           MultVal1((j-1)*5 + 5) = -s*c1(2*i-1)
           a1(i)=MultVal1((j-1)*5 + 3)
           b1(i)=0.0d0
           c1(i)=MultVal1((j-1)*5 + 5)
          else
           a1(i)=MultVal1((j-1)*5 + 2)
           b1(i)=MultVal1((j-1)*5 + 4)
           c1(i)=MultVal1((j-1)*5 + 4)
          endif
c MultVal complete
          if (i.eq.m) then
           RHS1(int(jv1*nx+(nx/2**j)+(iz-1)*(2*nx-1)))=
     +     RHS1(int(jv1*nx+(iz-1)*(2*nx-1)))-
     +     s*RHS1(int(jv1*nx-1+(iz-1)*(2*nx-1)))
          else 
           RHS1(int(jv1*nx+i+(iz-1)*(2*nx-1)))=
     +     RHS1(int(jv2*nx+2*i+(iz-1)*(2*nx-1)))
     +     -s*RHS1(int(jv2*nx+2*i-1+(iz-1)*(2*nx-1)))
     +     -t*RHS1(int(jv2*nx+2*i+1+(iz-1)*(2*nx-1)))               
          endif
         enddo
         MultVal1((k-1)*5 + 1)=MultVal1((k-1)*5 + 3)
         MultVal1((k-1)*5 + 3)=0.0d0
        enddo
        tnew2(nx+(iz-1)*nx)=
     +  RHS1(2*nx-1+(iz-1)*(2*nx-1))/MultVal1(5*(k-1)+1)
       enddo

  
!$omp parallel do shared(tnew2,RHS1) private(iz,j,i,jv2,ii)
    	 do iz=1,(nz-1)*ny
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do i=2**(j-1),nx-2**(j-1),2**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew2(i+(iz-1)*nx)=
     +     (RHS1(1+(iz-1)*(2*nx-1))-24.0d0*tnew2(2+(iz-1)*nx))/24.0d0
          elseif (i.eq.(2**(j-1))) then
           tnew2(i+(iz-1)*nx)=
     +     (RHS1(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew2(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+1)
          elseif ((i.eq.nx-2**(j-1)).and.(j.eq.2)) then
           tnew2(i+(iz-1)*nx)=
     +     (RHS1(3*nx/2-1+(iz-1)*(2*nx-1))-
     +      MultVal1(4)*tnew2(i-2+(iz-1)*nx)-
     +      MultVal1(4)*tnew2(i+2+(iz-1)*nx))/MultVal1(2)
          elseif (j.eq.1) then
           tnew2(i+(iz-1)*nx)=(RHS1(i+(iz-1)*(2*nx-1))-
     +     tnew2(i-1+(iz-1)*nx)-tnew2(i+1+(iz-1)*nx))/22.0d0
          else
           tnew2(i+(iz-1)*nx)=(RHS1(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +     MultVal1(5*(j-2)+4)*tnew2(i-2**(j-1)+(iz-1)*nx)-
     +     MultVal1(5*(j-2)+4)*tnew2(i+2**(j-1)+(iz-1)*nx))/
     +     MultVal1(5*(j-2)+2)
          endif
         enddo
        enddo
       enddo
!$omp end parallel do

!==================================================================== 
! Solves S12*tnew4=tnew4 using Fourier method 
!==================================================================== 
	elseif (Snum.eq.12) then

	 call zcopy(size(tnew4),tnew4,1,y4,1)
!$omp parallel do 
   	 do i=1,size(tnew4)
	  tnew4(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 

       pdn=pi/dfloat(nx)
       q=dsqrt(2.0d0/dfloat(nx))

!$omp parallel do shared(tnew4) private (k,i,sv,j)
       do k=1,ny*(nz-1)
        do i=1,nx-1
         sv=dcmplx(0.0d0,0.0d0)
         do j=1,nx-1
          sv=sv+y4((k-1)*(nx-1)+j)*q*dsin(i*j*pdn)
         enddo
         sv=sv*(1.0d0)/(22.0d0+2.0d0*dcos(i*pdn))
         do j=1,nx-1
          tnew4((k-1)*(nx-1)+j)=tnew4((k-1)*(nx-1)+j)+sv*q*dsin(i*j*pdn)
         enddo
        enddo
       enddo
!$omp end parallel do

!==================================================================== 
! Solves S13*tnew5=tnew5 using classic CR
!==================================================================== 
	elseif (Snum.eq.13) then
!$omp parallel do 
	 do i=1,size(RHS5)
	  RHS5(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do 
    	 do i=1,size(tnew5)
        RHS5(i)=tnew5(i)
       enddo
!$omp end parallel do
!$omp parallel do 
	 do i=1,size(tnew5)
	  tnew5(i)=dcmplx(0.0d0,0.0d0)
 	 enddo
!$omp end parallel do
!$omp parallel do 
	 do i=1,size(MultVal3)
	  MultVal3(i)=0.0d0
	 enddo
!$omp end parallel do 
    	 k=int(dlog(dfloat(nz))/dlog(2.0d0))
c	 k=lgz

!$omp parallel do 
       do iz=1,nz
        a3(iz)=22.0d0
        b3(iz)=1.0d0
        c3(iz)=1.0d0
       enddo
!$omp end parallel do 
         a3(1) = 24.0d0
    	 a3(nz) = 24.0d0
    	 b3(1) = 24.0d0
    	 b3(nz) = 0.0d0
    	 c3(1) = 0.0d0
    	 c3(nz) = 24.0d0
    	do j=1,k
        jv1=(2.0d0**j-1)/(2.0d0**(j-1))
        jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
c	2
        m=int(nz/2**j)
        do i=1,m
         s=c3(2*i)/a3(2*i-1)
         if (i.eq.m) then
          t=0.0d0 
          MultVal3((j-1)*5 + 3)=a3(2*i)-s*b3(2*i-1)
         else
          t=b3(2*i)/a3(2*i+1)
          if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
           MultVal3((j-1)*5 + i)=a3(2*i)-s*b3(2*i-1)-t*c3(2*i+1)
          endif
         endif
            
         if (i.eq.1) then
          MultVal3((j-1)*5 + 4)=-t*b3(2*i+1)
          a3(i)=MultVal3((j-1)*5 + 1)
          b3(i)=MultVal3((j-1)*5 + 4)
          c3(i)=0.0d0
         elseif (i.eq.m) then
          MultVal3((j-1)*5 + 5) = -s*c3(2*i-1)
          a3(i)=MultVal3((j-1)*5 + 3)
          b3(i)=0.0d0
          c3(i)=MultVal3((j-1)*5 + 5)
         else
          a3(i)=MultVal3((j-1)*5 + 2)
          b3(i)=MultVal3((j-1)*5 + 4)
          c3(i)=MultVal3((j-1)*5 + 4)
         endif

c Multval complete
         if (i.eq.m) then
          RHS5(int((jv1*nz+(nz/2**j)-1)*(nx-1)*ny+1):
     +         int((jv1*nz+(nz/2**j))*(nx-1)*ny))=
     +    RHS5(int((jv1*nz-1)*(nx-1)*ny+1):
     +         int((jv1*nz-1)*(nx-1)*ny+(nx-1)*ny))-
     +    s*RHS5(int((jv1*nz-2)*(nx-1)*ny+1):
     +           int((jv1*nz-2)*(nx-1)*ny+(nx-1)*ny))
         else
          RHS5(int((jv1*nz+i-1)*(nx-1)*ny+1):int((jv1*nz+i)*(nx-1)*ny))=
     +    RHS5(int((jv2*nz-1+2*i)*(nx-1)*ny+1):
     +         int((jv2*nz-1+2*i)*(nx-1)*ny+(nx-1)*ny))-
     +    s*RHS5(int((jv2*nz+2*i-2)*(nx-1)*ny+1):
     +           int((jv2*nz+2*i-2)*(nx-1)*ny+(nx-1)*ny))-
     +    t*RHS5(int((jv2*nz+2*i)*(nx-1)*ny+1):
     +           int((jv2*nz+2*i)*(nx-1)*ny+(nx-1)*ny))           
         endif
        enddo
        MultVal3((k-1)*5 + 1) = MultVal3((k-1)*5 + 3)
        MultVal3((k-1)*5 + 3) = 0.0d0
       enddo

       tnew5((nz-1)*(nx-1)*ny+1:nz*(nx-1)*ny)=
     + RHS5(2*(nz-1)*(nx-1)*ny+1:(2*nz-1)*(nx-1)*ny)*1.0d0/
     + MultVal3(5*(k-1)+1)

    	do j=k,1,-1
       jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
       ii=0
       do i=2**(j-1),nz-2**(j-1),2**j
        ii=ii+1
        if ((i.eq.1).and.(j.eq.1)) then
         tnew5(1:(nx-1)*ny)=(RHS5(1:(nx-1)*ny)-
     +   24.0d0*tnew5((nx-1)*ny+1:2*(nx-1)*ny))/24.0d0
        elseif (i.eq.(2**(j-1))) then
         tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)=
     +   (RHS5(int((nz*jv2+2*ii-2)*(nx-1)*ny+1):
     +         int((nz*jv2+2*ii-1)*(nx-1)*ny))-
     +   MultVal3(5*(j-2)+4)*
     +   tnew5((i-1+2**(j-1))*(nx-1)*ny+1:(i+2**(j-1))*(nx-1)*ny))/
     +   MultVal3(5*(j-2)+1)
        elseif ((i.eq.(nz-2**(j-1))).and.(j.eq.2)) then
         tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)=
     +   (RHS5((3*nz/2-2)*(nx-1)*ny+1:(3*nz/2-1)*(nx-1)*ny)-
     +   MultVal3(4)*tnew5((i-3)*(nx-1)*ny+1:(i-2)*(nx-1)*ny)-
     +   MultVal3(4)*tnew5((i+1)*(nx-1)*ny+1:(i+2)*(nx-1)*ny))/
     +   MultVal3(2)
        elseif (j.eq.1) then
         tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)=
     +   (RHS5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)-
     +   tnew5((i-2)*(nx-1)*ny+1:(i-1)*(nx-1)*ny)-
     +   tnew5(i*(nx-1)*ny+1:(i+1)*(nx-1)*ny))/22.0d0
        else
         tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)=
     +   (RHS5(int((nz*jv2+2*ii-2)*(nx-1)*ny+1):
     +         int((nz*jv2+2*ii-1)*(nx-1)*ny))-
     +   MultVal3(5*(j-2)+4)*
     +   tnew5((i-2**(j-1)-1)*(nx-1)*ny+1:(i-2**(j-1))*(nx-1)*ny)-
     +   MultVal3(5*(j-2)+4)*
     +   tnew5((i+2**(j-1)-1)*(nx-1)*ny+1:(i+2**(j-1))*(nx-1)*ny))/
     +   MultVal3(5*(j-2)+2)
        endif
       enddo
      enddo

!==================================================================== 
! Solves S14*tnew4=tnew4 using Fourier method
!==================================================================== 
	elseif (Snum.eq.14) then
        
	 call zcopy(size(tnew4),tnew4,1,y4,1)
c	 print*,dznrm2(size(y4),y4,1)
!$omp parallel do 
	 do i=1,size(tnew4)
	  tnew4(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 

    	 pdn=pi/dfloat(nz)
    	 q=dsqrt(2.0d0/dfloat(nz))

         do iz=1,nz-1
!$omp parallel do 
          do i=1,ny*(nx-1)
	   sv2(i)=dcmplx(0.0d0,0.0d0)
	  enddo
!$omp end parallel do 
!$omp critical
          do j=1,nz-1
           sv2(1:(nx-1)*ny)=sv2(1:(nx-1)*ny)+
     +     y4((j-1)*ny*(nx-1)+1:j*ny*(nx-1))*q*dsin(iz*j*pdn)
          enddo
        sv2(1:(nx-1)*ny)=
     +  sv2(1:(nx-1)*ny)*(1.0d0/(22.0d0+2.0d0*dcos(iz*pdn)))
!$omp end critical

!$omp parallel do shared(tnew4)
        do j=1,nz-1
         tnew4((j-1)*(ny)*(nx-1)+1:j*ny*(nx-1))=
     +   tnew4((j-1)*ny*(nx-1)+1:j*ny*(nx-1))+
     +   sv2(1:(nx-1)*ny)*q*dsin(iz*j*pdn)
        enddo
!$omp end parallel do 
    	 enddo

!==================================================================== 
! Solves S15*tnew5=tnew5 using classic CR
!==================================================================== 
	elseif (Snum.eq.15) then
!$omp parallel do 
	 do i=1,size(RHS4)
	  RHS4(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do 
    	 do i=1,nz
        RHS4((i-1)*(2*ny-1)*(nx-1)+1:(i-1)*(2*ny-1)*(nx-1)+(nx-1)*ny)=
     +     tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)
    	 enddo
!$omp end parallel do
!$omp parallel do 
 	 do i=1,size(tnew5)
	  tnew5(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 
!$omp parallel do 
       do i=1,size(MultVal2)
	  MultVal2(i)=0.0d0
	 enddo
!$omp end parallel do 
	 k = int(dlog(dfloat(ny))/dlog(2.0d0))

    	 do iz=1,nz
        do i=1,ny
         a2(i)=22.0d0
         b2(i)=1.0d0
         c2(i)=1.0d0
        enddo
        a2(1)=24.0d0
        a2(ny)=24.0d0
        b2(1)=24.0d0
        b2(ny)=0.0d0
        c2(1)=0.0d0
        c2(ny)=24.0d0

        do j=1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(ny/2**j)
         do i=1,m
          s=c2(2*i)/a2(2*i-1)
          if (i.eq.m) then
           t=0.0d0 
           MultVal2((j-1)*5+3)=a2(2*i)-s*b2(2*i-1)
          else
           t=b2(2*i)/a2(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
             MultVal2((j-1)*5 + i)=a2(2*i)-s*b2(2*i-1)-t*c2(2*i+1)
           endif
          endif
            
          if (i.eq.1) then
           MultVal2((j-1)*5 + 4) = -t*b2(2*i+1)
           a2(i)=MultVal2((j-1)*5 + 1)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=0.0d0
          elseif (i.eq.m) then
           MultVal2((j-1)*5 + 5) = -s*c2(2*i-1)
           a2(i)=MultVal2((j-1)*5 + 3)
           b2(i)=0.0d0
           c2(i)=MultVal2((j-1)*5 + 5)
          else
           a2(i)=MultVal2((j-1)*5 + 2)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=MultVal2((j-1)*5 + 4)
          endif
c Multval complete
          if (i.eq.m) then
           RHS4(int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +          int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)+0)*(nx-1)))=
     +     RHS4(int((jv1*ny+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +          int((jv1*ny+(iz-1)*(2*ny-1)-1)*(nx-1)+(nx-1)))-
     +     s*RHS4(int((jv1*ny+(iz-1)*(2*ny-1)-2)*(nx-1)+1):
     +            int((jv1*ny+(iz-1)*(2*ny-1)-2)*(nx-1)+(nx-1)))
          else
           RHS4(int((jv1*ny+i+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +          int((jv1*ny+i+(iz-1)*(2*ny-1))*(nx-1)))=
     +     RHS4(int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*(nx-1)+1):
     +          int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))-
     +     s*RHS4(int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*(nx-1)+1):
     +            int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))-
     +     t*RHS4(int((jv2*ny+2*i+(iz-1)*(2*ny-1))*(nx-1)+1):
     +            int((jv2*ny+2*i+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))   
          endif
         enddo
        MultVal2((k-1)*5 + 1) = MultVal2((k-1)*5 + 3)
        MultVal2((k-1)*5 + 3) = 0.0d0
        enddo

c tnew evaluation
        tnew5(((ny-1)+(iz-1)*ny)*(nx-1)+1:
     +        ((ny-1)+(iz-1)*ny)*(nx-1)+(nx-1))=
     +  RHS4((nx-1)*(2*ny-2+(iz-1)*(2*ny-1))+1:
     +       (nx-1)*(2*ny-1+(iz-1)*(2*ny-1)))*1.0d0/
     +  MultVal2(5*(k-1)+1)
    	 enddo

!$omp parallel do shared(tnew5,RHS4) private(iz,j,i,jv2,ii)
    	 do iz=1,nz
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do i=2**(j-1),ny-2**(j-1),2**j
          ii=ii+1
          if ((i.eq.1).and.(j.eq.1)) then
           tnew5((nx-1)*(iz-1)*ny+1:(nx-1)*(iz-1)*ny+(nx-1))=
     +     (RHS4(1+(nx-1)*(iz-1)*(2*ny-1):(nx-1)+(nx-1)*(iz-1)*(2*ny-1))
     +     -24.0d0*tnew5(1+(nx-1)+(nx-1)*(iz-1)*ny:
     +           2*(nx-1)+(nx-1)*(iz-1)*ny))*1.0d0/24.0d0
          elseif (i.eq.(2**(j-1))) then
           tnew5((i-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:
     +     (i-1)*(nx-1)+(iz-1)*ny*(nx-1)+(nx-1))=1.0d0/
     +     MultVal2(5*(j-2)+1)*
     +     (RHS4(int((nx-1)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+1):
     +           int((nx-1)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+(nx-1)))-
     +     MultVal2(5*(j-2)+4)*
     +     tnew5(((i-1)+2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +           ((i-1)+2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1)))
          elseif ((i.eq.(ny-2**(j-1))).and.(j.eq.2)) then
           tnew5(((i-1)+(iz-1)*ny)*(nx-1)+1:
     +           ((i-1)+(iz-1)*ny)*(nx-1)+(nx-1))=
     +      (RHS4((3*ny/2-2+(iz-1)*(2*ny-1))*(nx-1)+1:
     +            (3*ny/2-2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1))-
     +      MultVal2(4)*tnew5((i-3+(iz-1)*ny)*(nx-1)+1:
     +                      (i-3+(iz-1)*ny)*(nx-1)+(nx-1))-
     +      MultVal2(4)*tnew5((i+1+(iz-1)*ny)*(nx-1)+1:
     +                        (i+1+(iz-1)*ny)*(nx-1)+(nx-1)))*1.0d0/
     +      MultVal2(2)
          elseif (j.eq.1) then
           tnew5(((i-1)+(iz-1)*ny)*(nx-1)+1:
     +           (i-1+(iz-1)*ny)*(nx-1)+(nx-1))=
     +     (RHS4((i-1+(iz-1)*(2*ny-1))*(nx-1)+1:
     +          (i-1+(iz-1)*(2*ny-1))*(nx-1)+(nx-1))-
     +     tnew5((i-2+(iz-1)*ny)*(nx-1)+1:
     +           (i-2+(iz-1)*ny)*(nx-1)+(nx-1))-
     +     tnew5((i+(iz-1)*ny)*(nx-1)+1:(i+(iz-1)*ny)*(nx-1)+(nx-1)))*
     +           1.0d0/22.0d0
          else
           tnew5((i-1+(iz-1)*ny)*(nx-1)+1:
     +           (i-1+(iz-1)*ny)*(nx-1)+(nx-1))=(RHS4(
     +     int((ny*jv2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)*(2*ii-2)+1):
     +     int((ny*jv2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)*(2*ii-2)+(nx-1)))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew5((i-1-2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +          (i-1-2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew5((i-1+2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +           (i-1+2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1)))*1.0d0/
     +      MultVal2(5*(j-2)+2)
          endif
         enddo
        enddo
       enddo
!$omp end parallel do
 
!==================================================================== 
! Solves S16*tnew4=tnew4 using Fourier method
!==================================================================== 
	elseif (Snum.eq.16) then

	 call zcopy(size(tnew4),tnew4,1,y4,1)
!$omp parallel do 
	 do i=1,size(tnew4)
	  tnew4(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 
       pdn=pi/dfloat(nz)
       q=dsqrt(2.0d0/dfloat(nz))

       do i=1,nz-1
!$omp parallel do 
        do ii=1,ny*(nx-1)
	  sv2(ii)=dcmplx(0.0d0,0.0d0)
	enddo
!$omp end parallel do 
!$omp critical 
        do j=1,nz-1
         sv2(1:(nx-1)*ny)=sv2(1:(nx-1)*ny)+
     +   y4((j-1)*ny*(nx-1)+1:j*ny*(nx-1))*q*dsin(i*j*pdn)
        enddo
        sv2(1:(nx-1)*ny)=
     +  sv2(1:(nx-1)*ny)*(1.0d0/(22.0d0+2.0d0*dcos(i*pdn)))
!$omp end critical 
!$omp parallel do shared(tnew4) 
        do j=1,nz-1
         tnew4((j-1)*ny*(nx-1)+1:j*ny*(nx-1))=
     +   tnew4((j-1)*ny*(nx-1)+1:j*ny*(nx-1))+
     +   sv2(1:(nx-1)*ny)*q*dsin(i*j*pdn)
        enddo
!$omp end parallel do 
       enddo
	 
!==================================================================== 
! Solves S17*tnew3=tnew3 using classic CR
!==================================================================== 
	elseif (Snum.eq.17) then

!$omp parallel do 
	 do i=1,size(RHS3)
	  RHS3(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do
!$omp parallel do 
	 do ii=1,(ny-1)*nz
    	  do i=1,nx
         RHS3(i+(ii-1)*(2*nx-1))=tnew3(i+(ii-1)*nx)
        enddo
	 enddo
!$omp end parallel do
!$omp parallel do 
	 do i=1,size(tnew3)
	  tnew3(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 
!$omp parallel do 
	 do i=1,size(MultVal1)
	  MultVal1(i)=0.0d0
	 enddo
!$omp end parallel do 

      k=int(dlog(dfloat(nx))/dlog(2.0d0))


	 do iz=1,(ny-1)*nz
    	  do i=1,nx
         a1(i)=22.0d0
         b1(i)=1.0d0
         c1(i)=1.0d0
    	  enddo
    	  a1(1)=24.0d0
    	  a1(nx)=24.0d0
        b1(1)=24.0d0
        b1(nx)=0.0d0
        c1(1)=0.0d0
        c1(nx)=24.0d0
 
        do j=1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
          m=int(nx/2**j)
 
          do i=1,m
           s=c1(2*i)/a1(2*i-1)
           if (i.eq.m) then
            t=0.0d0
            MultVal1((j-1)*5 + 3)=a1(2*i)-s*b1(2*i-1)
           else
            t=b1(2*i)/a1(2*i+1)
            if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
             MultVal1((j-1)*5+i)=a1(2*i)-s*b1(2*i-1)-t*c1(2*i+1)
            endif
           endif
            
           if (i.eq.1) then
            MultVal1((j-1)*5 + 4) = -t*b1(2*i+1)
            a1(i) = MultVal1((j-1)*5 + 1)
            b1(i) = MultVal1((j-1)*5 + 4)
            c1(i) = 0.0d0
           elseif (i.eq.m) then
            MultVal1((j-1)*5 + 5) = -s*c1(2*i-1)
            a1(i) = MultVal1((j-1)*5 + 3)
            b1(i) = 0.0d0
            c1(i) = MultVal1((j-1)*5 + 5)
           else
            a1(i) = MultVal1((j-1)*5 + 2)
            b1(i) = MultVal1((j-1)*5 + 4)
            c1(i) = MultVal1((j-1)*5 + 4)
           endif
c MultVal complete
           if (i.eq.m) then
            RHS3(int(jv1*nx+(nx/2**j)+(iz-1)*(2*nx-1)))=
     +      RHS3(int(jv1*nx+(iz-1)*(2*nx-1)))-
     +      s*RHS3(int(jv1*nx-1+(iz-1)*(2*nx-1)))
           else 
            RHS3(int(jv1*nx+i+(iz-1)*(2*nx-1)))=
     +      RHS3(int(jv2*nx+2*i+(iz-1)*(2*nx-1)))-
     +      s*RHS3(int(jv2*nx+2*i-1+(iz-1)*(2*nx-1)))-
     +      t*RHS3(int(jv2*nx+2*i+1+(iz-1)*(2*nx-1)))              
           endif
          enddo
 
          MultVal1((k-1)*5 + 1) = MultVal1((k-1)*5 + 3)
          MultVal1((k-1)*5 + 3) = 0.0d0
         enddo
         tnew3(nx+(iz-1)*nx)=RHS3(2*nx-1+(iz-1)*(2*nx-1))/
     +   MultVal1(5*(k-1)+1)
	  enddo


!$omp parallel do shared(tnew3,RHS3) private(iz,j,i,jv2,ii)
	  do iz=1,(ny-1)*nz
	   do j=k,1,-1
          jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
          ii=0
          do i=2**(j-1),nx-2**(j-1),2**j
           ii=ii+1
           if ((i.eq.1).and.(j.eq.1)) then
            tnew3(i+(iz-1)*nx)=(RHS3(1+(iz-1)*(2*nx-1))-24.0d0*
     +      tnew3(2+(iz-1)*nx))/24.0d0
           elseif (i.eq.2**(j-1)) then
           tnew3(i+(iz-1)*nx)=(RHS3(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +      MultVal1(5*(j-2)+4)*tnew3(i+2**(j-1)+(iz-1)*nx))/
     +      MultVal1(5*(j-2)+1)
           elseif ((i.eq.nx-2**(j-1)).and.(j.eq.2)) then
            tnew3(i+(iz-1)*nx)=(RHS3(3*nx/2-1+(iz-1)*(2*nx-1))-
     +      MultVal1(4)*tnew3(i-2+(iz-1)*nx)-
     +      MultVal1(4)*tnew3(i+2+(iz-1)*nx))/MultVal1(2)
           elseif (j.eq.1) then
            tnew3(i+(iz-1)*nx)=(RHS3(i+(iz-1)*(2*nx-1))-
     +      tnew3(i-1+(iz-1)*nx)-tnew3(i+1+(iz-1)*nx))/22.0d0
           else
           tnew3(i+(iz-1)*nx)=(RHS3(int(nx*jv2+2*ii-1+(iz-1)*(2*nx-1)))-
     +      MultVal1(5*(j-2)+4)*tnew3(i-2**(j-1)+(iz-1)*nx)-
     +      MultVal1(5*(j-2)+4)*tnew3(i+2**(j-1)+(iz-1)*nx))/
     +      MultVal1(5*(j-2)+2)
           endif
    	    enddo
	   enddo
	  enddo
!$omp end parallel do

!==================================================================== 
! Solves S18*tnew6=tnew6 using Fourier method
!==================================================================== 
	elseif (Snum.eq.18) then

	call zcopy(size(tnew6),tnew6,1,y6,1)
!$omp parallel do 
      do i=1,size(tnew6)
    	 tnew6(i)=dcmplx(0.0d0,0.0d0)
	enddo
!$omp end parallel do 

    	pdn=pi/dfloat(nx)
    	q=dsqrt(2.0d0/dfloat(nx))

!$omp parallel do shared(tnew6) private(k,i,sv,j) 
       do k=1,(ny-1)*nz
        do i=1,nx-1
         sv=dcmplx(0.0d0,0.0d0)
         do j=1,nx-1
          sv=sv+y6((k-1)*(nx-1)+j)*q*dsin(i*j*pdn)
         enddo
         sv=sv*1.0d0/(22.0d0+2.0d0*dcos(i*pdn))
         do j=1,nx-1
          tnew6((k-1)*(nx-1)+j)=tnew6((k-1)*(nx-1)+j)+sv*q*dsin(i*j*pdn)
         enddo
        enddo
       enddo
!$omp end parallel do 
!==================================================================== 
! Solves S19*tnew5=tnew5 using classic CR
!==================================================================== 
	elseif (Snum.eq.19) then

!$omp parallel do 
       do i=1,size(RHS4)
	  RHS4(i)=dcmplx(0.0d0,0.0d0)
       enddo
!$omp end parallel do 

!$omp parallel do 
       do i=1,nz
        RHS4((i-1)*(2*ny-1)*(nx-1)+1:(i-1)*(2*ny-1)*(nx-1)+(nx-1)*ny)=
     +          tnew5((i-1)*(nx-1)*ny+1:i*(nx-1)*ny)
       enddo
!$omp end parallel do 
c	 print*,'HERERHS4_1',RHS4,size(RHS4)
!$omp parallel do 
	 do i=1,size(tnew5)
	  tnew5(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 

!$omp parallel do 
	 do i=1,size(MultVal2)
	  MultVal2(i)=0.0d0
	 enddo
!$omp end parallel do 

 	 k=int(dlog(dfloat(ny))/dlog(2.0d0))
c        k=lgy

       do iz=1,nz
        do i=1,ny
         a2(i)=22.0d0
         b2(i)=1.0d0
         c2(i)=1.0d0
        enddo
        a2(1)=24.0d0
        a2(ny)=24.0d0
        b2(1)=24.0d0
        b2(ny)=0.0d0
        c2(1)=0.0d0
        c2(ny)=24.0d0

        do j=1,k
         jv1=(2.0d0**j-1)/(2.0d0**(j-1))
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         m=int(ny/2**j)
         do i=1,m
          s=c2(2*i)/a2(2*i-1)
          if (i.eq.m) then
           t=0.0d0 
           MultVal2((j-1)*5 + 3)=a2(2*i)-s*b2(2*i-1)
          else
           t=b2(2*i)/a2(2*i+1)
           if ((i.eq.1).or.(i.eq.2).and.(j.ne.k)) then
            MultVal2((j-1)*5 + i)=a2(2*i)-s*b2(2*i-1)-t*c2(2*i+1)
           endif
          endif

          if (i.eq.1) then
           MultVal2((j-1)*5 + 4)=-t*b2(2*i+1)
           a2(i)=MultVal2((j-1)*5 + 1)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=0.0d0
          elseif (i.eq.m) then
           MultVal2((j-1)*5 + 5) = -s*c2(2*i-1)
           a2(i)=MultVal2((j-1)*5 + 3)
           b2(i)=0.0d0
           c2(i)=MultVal2((j-1)*5 + 5)
          else
           a2(i)=MultVal2((j-1)*5 + 2)
           b2(i)=MultVal2((j-1)*5 + 4)
           c2(i)=MultVal2((j-1)*5 + 4)
          endif

          if (i.eq.m) then
           RHS4(int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +          int((jv1*ny+(ny/2**j)+(iz-1)*(2*ny-1)+0)*(nx-1)))= 
     +     RHS4(int((jv1*ny+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +          int((jv1*ny+(iz-1)*(2*ny-1)-1)*(nx-1)+(nx-1)))-
     +     s*RHS4(int((jv1*ny+(iz-1)*(2*ny-1)-2)*(nx-1)+1):
     +          int((jv1*ny+(iz-1)*(2*ny-1)-2)*(nx-1)+(nx-1)))
          else
           RHS4(int((jv1*ny+i+(iz-1)*(2*ny-1)-1)*(nx-1)+1):
     +           int((jv1*ny+i+(iz-1)*(2*ny-1))*(nx-1)))=
     +     RHS4(int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*(nx-1)+1):
     +           int((jv2*ny-1+2*i+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))-
     +     s*RHS4(int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*(nx-1)+1):
     +           int((jv2*ny+2*i-2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))-
     +     t*RHS4(int((jv2*ny+2*i+(iz-1)*(2*ny-1))*(nx-1)+1):
     +           int((jv2*ny+2*i+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)))       
          endif
         enddo
         MultVal2((k-1)*5 + 1) = MultVal2((k-1)*5 + 3)
         MultVal2((k-1)*5 + 3) = 0.0d0

        enddo
        tnew5(((ny-1)+(iz-1)*ny)*(nx-1)+1:
     +       ((ny-1)+(iz-1)*ny)*(nx-1)+(nx-1))=
     +  RHS4((nx-1)*(2*ny-2+(iz-1)*(2*ny-1))+1:
     +       (nx-1)*(2*ny-1+(iz-1)*(2*ny-1)))*1.0d0/MultVal2(5*(k-1)+1)
    	 enddo
!$omp parallel do shared(tnew5,RHS4) private(iz,j,i,jv2,ii)
    	 do iz=1,nz
        do j=k,1,-1
         jv2=(2.0d0**(j-1)-1)/(2.0d0**(j-2))
         ii=0
         do i=2**(j-1),ny-2**(j-1),2**j
          ii=ii+1
           if ((i.eq.1).and.(j.eq.1)) then
            tnew5((nx-1)*(iz-1)*ny+1:
     +            (nx-1)*(iz-1)*ny+(nx-1))=
     +      (RHS4(1+(nx-1)*(iz-1)*(2*ny-1):
     +              (nx-1)+(nx-1)*(iz-1)*(2*ny-1))-24.0d0*
     +      tnew5(1+(nx-1)+(nx-1)*(iz-1)*ny:
     +            2*(nx-1)+(nx-1)*(iz-1)*ny))*1.0d0/24.0d0
           elseif (i.eq.(2**(j-1))) then
            tnew5((i-1)*(nx-1)+(iz-1)*ny*(nx-1)+1:
     +            (i-1)*(nx-1)+(iz-1)*ny*(nx-1)+(nx-1))=1.0d0/
     +      MultVal2(5*(j-2)+1)*
     +      (RHS4(int((nx-1)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+1):
     +            int((nx-1)*(ny*jv2+(ii-1)+(iz-1)*(2*ny-1))+(nx-1)))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew5(((i-1)+2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +            ((i-1)+2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1)))   
           elseif ((i.eq.(ny-2**(j-1))).and.(j.eq.2)) then
            tnew5(((i-1)+(iz-1)*ny)*(nx-1)+1:
     +            ((i-1)+(iz-1)*ny)*(nx-1)+(nx-1))=
     +      (RHS4((3*ny/2-2+(iz-1)*(2*ny-1))*(nx-1)+1:
     +            (3*ny/2-2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1))-
     +       MultVal2(4)*tnew5((i-3+(iz-1)*ny)*(nx-1)+1:
     +                       (i-3+(iz-1)*ny)*(nx-1)+(nx-1))-
     +       MultVal2(4)*tnew5((i+1+(iz-1)*ny)*(nx-1)+1:
     +                         (i+1+(iz-1)*ny)*(nx-1)+(nx-1)))*1.0d0/
     +       MultVal2(2)
           elseif (j.eq.1) then
            tnew5(((i-1)+(iz-1)*ny)*(nx-1)+1:
     +             (i-1+(iz-1)*ny)*(nx-1)+(nx-1))=
     +      (RHS4((i-1+(iz-1)*(2*ny-1))*(nx-1)+1:
     +           ( i-1+(iz-1)*(2*ny-1))*(nx-1)+(nx-1))-
     +      tnew5((i-2+(iz-1)*ny)*(nx-1)+1:
     +            (i-2+(iz-1)*ny)*(nx-1)+(nx-1))-
     +      tnew5((i+(iz-1)*ny)*(nx-1)+1:
     +            (i+(iz-1)*ny)*(nx-1)+(nx-1)))*1.0d0/22.0d0        
           else
            tnew5((i-1+(iz-1)*ny)*(nx-1)+1:
     +           (i-1+(iz-1)*ny)*(nx-1)+(nx-1))=
     +     (RHS4(int((ny*jv2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)*(2*ii-2)+1):
     +     int((ny*jv2+(iz-1)*(2*ny-1))*(nx-1)+(nx-1)*(2*ii-2)+(nx-1)))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew5((i-1-2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +           (i-1-2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1))-
     +      MultVal2(5*(j-2)+4)*
     +      tnew5((i-1+2**(j-1)+(iz-1)*ny)*(nx-1)+1:
     +            (i-1+2**(j-1)+(iz-1)*ny)*(nx-1)+(nx-1)))*1.0d0/
     +      MultVal2(5*(j-2)+2)
           endif
          enddo
         enddo
        enddo
!$omp end parallel do

!==================================================================== 
! Solves S20*tnew6=tnew6 using Fourier method
!==================================================================== 
	elseif (Snum.eq.20) then
         call zcopy(size(tnew6),tnew6,1,y6,1)
!$omp parallel do 
	 do i=1,size(tnew6)
 	  tnew6(i)=dcmplx(0.0d0,0.0d0)
	 enddo
!$omp end parallel do 
    	 q=dsqrt(2.0d0/dfloat(ny))
    	 pdn=pi/dfloat(ny)
!$omp parallel do shared(tnew6) private (k,i,iz,sv,j)
       do iz=0,nz-1
        do k=1,nx-1
         do i=1,ny-1
          sv=dcmplx(0.0d0,0.0d0)
          do j=0,ny-2
           sv=sv+y6(j*(nx-1)+iz*(nx-1)*(ny-1)+k)*q*dsin(i*(j+1)*pdn)
          enddo
          sv=sv*1.0d0/(22.0d0+2.0d0*dcos(i*pdn))
          do j=1,ny-1
           tnew6(k+(j-1)*(nx-1)+iz*(nx-1)*(ny-1))=
     +     tnew6(k+(j-1)*(nx-1)+iz*(nx-1)*(ny-1))+sv*q*dsin(i*j*pdn)
          enddo
         enddo
        enddo
       enddo
!$omp end parallel do


	 
	endif
	return
	end subroutine

