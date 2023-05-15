! ===================================================================
! Title: bicgstabU.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! bicgstabU.f implements the BiCGSTAB iterative method for the
! solution of the linear system UE=b
!====================================================================      
! Input:
! n: integer, total number of unknown elements in the grid
! nx: integer,number of elements in x-direction
! ny: integer,number of elements in y-direction
! nz: integer,number of elements in z-direction
! hx: real, discretization step in x-direction
! hy: real, discretization step in y-direction
! hz: real, discretization step in z-direction
! e1: complex array, dimension: nx*(ny-1)*(nz-1)
! e2: complex array, dimension: (nx-1)*ny*(nz-1) 
! e3: complex array, dimension: (nx-1)*(ny-1)*nz
! F1: complex array, dimension: nx*(ny-1)*(nz-1)
! F2: complex array, dimension: (nx-1)*ny*(nz-1)  
! maxstepU: integer, maximum number of iterations to be performed
! tolU: real, tolerance of the method
! Local auxiliary variables:
! v_i,r_i1,rh,p_i,t,s,x_i1,b0: complex arrays, dimension n
! Auxiliary variables used by other subroutines:
! pt1,pt2,pt3,pt4,pt5,pt6,st1,st2,st3,st4,st5,st6
! RHS1,RHS2,RHS3,RHS4,RHS5,RHS6
! MultVal1,MultVal2,MultVal3
! a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2
! temp1,temp11,temp2,temp22,temp3,temp33
! tempEX,tempEY,tempEZ
! tempEX_1,tempEY_1,tempEZ_1
! tempEX_2,tempEY_2,tempEZ_2
! tempEX_3,tempEY_3,tempEZ_3
! Output:
! x: complex array, dimension: n
!==================================================================== 
      subroutine bicgstabU(n,nx,ny,nz,hx,hy,hz,ci,
     +           e1,e2,e3,F1,F2,
     +           temp1,temp11,temp2,temp22,temp3,temp33,
     +           pt1,pt2,pt3,pt4,pt5,pt6,
     +           st1,st2,st3,st4,st5,st6,
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +           maxstepU,tolU,
     +           tempEX,tempEY,tempEZ,
     +           tempEX_1,tempEY_1,tempEZ_1,
     +           tempEX_2,tempEY_2,tempEZ_2,
     +           tempEX_3,tempEY_3,tempEZ_3,
     +           v_i,r_i1,rh,p_i,t,s,x_i1,b0,
     +           x)

      implicit none
      integer::n,nx,ny,nz,maxstepU,i,istep,k
      real*8 ::hx,hy,hz,dznrm2,errU,
     + MultVal1(5*int(dlog(dfloat(nx))/dlog(2.0d0))),
     + MultVal2(5*int(dlog(dfloat(ny))/dlog(2.0d0))),
     + MultVal3(5*int(dlog(dfloat(nz))/dlog(2.0d0))),
     + a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     + t1,t2,tolU,nr
      complex*16 ::
     + e1(nx*(ny-1)*(nz-1)),e2((nx-1)*ny*(nz-1)),e3((nx-1)*(ny-1)*nz),
     + F1(nx*(ny-1)*(nz-1)),F2((nx-1)*ny*(nz-1)),
     + EX(nx*(ny-1)*(nz-1)),tempEX(nx*(ny-1)*(nz-1)),ci, 
     + EY((nx-1)*ny*(nz-1)),tempEY((nx-1)*ny*(nz-1)),
     + EZ((nx-1)*(ny-1)*nz),tempEZ((nx-1)*(ny-1)*nz),
     + pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     + pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     + pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     + pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     + pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     + pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     + temp1(nx*(ny-1)*(nz-1)),temp11(nx*(ny-1)*(nz-1)),
     + temp2((nx-1)*ny*(nz-1)),temp22((nx-1)*ny*(nz-1)),
     + temp3((nx-1)*(ny-1)*nz),temp33((nx-1)*(ny-1)*nz),
     + tempEX_1(nx*(ny-1)*(nz-1)),tempEX_2(nx*(ny-1)*(nz-1)),
     + tempEY_1((nx-1)*ny*(nz-1)),tempEY_2((nx-1)*ny*(nz-1)),
     + tempEZ_1((nx-1)*(ny-1)*nz),tempEZ_2((nx-1)*(ny-1)*nz),
     + tempEX_3(nx*(ny-1)*(nz-1)),tempEY_3((nx-1)*ny*(nz-1)),
     + tempEZ_3((nx-1)*(ny-1)*nz),
     + RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     + RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     + RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     + w_i,a_i,b_i1,po_i1,po_i2,zdotc,
     + sv1(nx*(ny-1)),sv2((nx-1)*ny),
     + v_i(n),r_i1(n),rh(n),p_i(n),b(n),t(n),s(n),x_i1(n),x(n),
     + b0(n)

       istep=1
       errU=1.0d0
       
       k=1
       do i=1,nx*(ny-1)*(nz-1)
        b(i)=F1(k)
        k=k+1
       enddo
       k=1
       do i=nx*(ny-1)*(nz-1)+1,nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)
        b(i)=F2(k)
        k=k+1
       enddo
       do i=nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+1,n
        b(i)=dcmplx(0.0d0,0.0d0)
       enddo

       nr=dznrm2(n,b,1)
       call U_Mult(b,n,nx,ny,nz,hx,hy,hz,ci,
     +       e1,e2,e3,tempEX,tempEY,tempEZ,
     +       pt1,pt2,pt3,pt4,pt5,pt6,
     +       st1,st2,st3,st4,st5,st6,
     +       RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +       MultVal1,MultVal2,MultVal3,
     +       a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +       temp1,temp11,temp2,temp22,temp3,temp33,
     +       tempEX_1,tempEX_2,tempEX_3,
     +       tempEY_1,tempEY_2,tempEY_3,
     +       tempEZ_1,tempEZ_2,tempEZ_3,
     +       p_i)
c      print*,'Norm p_i',dznrm2(size(p_i),p_i,1),size(p_i)

!$OMP PARALLEL DO        
        do i=1,n
         rh(i)=b(i)-p_i(i)
        enddo
!$OMP END PARALLEL DO        

        call zcopy(n,rh,1,r_i1,1)
        call zcopy(n, b,1,x_i1,1)

        do while ((errU.gt.tolU).and.(istep.lt.maxstepU))
         po_i1=zdotc(n,rh,1,r_i1,1)
         if (cdabs(po_i1).eq.(0.0d0)) then
          print*,'BiCGSTAB U fails'
          exit
         endif
         if (istep.eq.1) then
          call zcopy(n,rh,1,p_i,1)
         else
          b_i1=(po_i1/po_i2)*(a_i/w_i)
!$OMP PARALLEL DO        
         do i=1,n
          p_i(i)=r_i1(i)+b_i1*(p_i(i)-w_i*v_i(i))
         enddo
!$OMP END PARALLEL DO        
         endif
        
         call U_Mult(p_i,n,nx,ny,nz,hx,hy,hz,ci,
     +        e1,e2,e3,tempEX,tempEY,tempEZ,
     +        pt1,pt2,pt3,pt4,pt5,pt6,
     +        st1,st2,st3,st4,st5,st6,
     +        RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +        MultVal1,MultVal2,MultVal3,
     +        a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +        temp1,temp11,temp2,temp22,temp3,temp33,
     +        tempEX_1,tempEX_2,tempEX_3,
     +        tempEY_1,tempEY_2,tempEY_3,
     +        tempEZ_1,tempEZ_2,tempEZ_3,
     +        v_i)

       a_i=po_i1/zdotc(n,rh,1,v_i,1)
!$OMP PARALLEL DO
       do i=1,n
        s(i)=r_i1(i)-a_i*v_i(i)
       enddo
!$OMP END PARALLEL DO
       if (dznrm2(size(s),s,1).eq.(0.0)) then
!$OMP PARALLEL DO
        do i=1,n
         x(i)=x_i1(i)+a_i*p_i(i)
        enddo
!$OMP END PARALLEL DO
        print*,'BiCGSTAB U exits for norm(s)'
        exit
       endif

      call U_Mult(s,n,nx,ny,nz,hx,hy,hz,ci,
     +        e1,e2,e3,tempEX,tempEY,tempEZ,
     +        pt1,pt2,pt3,pt4,pt5,pt6,
     +        st1,st2,st3,st4,st5,st6,
     +        RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,
     +        MultVal1,MultVal2,MultVal3,
     +        a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +        temp1,temp11,temp2,temp22,temp3,temp33,
     +        tempEX_1,tempEX_2,tempEX_3,
     +        tempEY_1,tempEY_2,tempEY_3,
     +        tempEZ_1,tempEZ_2,tempEZ_3,
     +        t)

       w_i=zdotc(n,t,1,s,1)/zdotc(n,t,1,t,1)
       if (cdabs(w_i).eq.(0.0d0)) then
        print*,'BiCGSTAB U exits for wi'
        exit
       endif

!$OMP PARALLEL
!$OMP DO        
       do i=1,n
        x(i)=x_i1(i)+a_i*p_i(i)+w_i*s(i)
       enddo
!$OMP DO        
       do i=1,n
        r_i1(i)=s(i)-w_i*t(i)
       enddo
!$OMP END PARALLEL        

       call zcopy(n,x,1,x_i1,1)
       po_i2=po_i1
       errU=dznrm2(size(r_i1),r_i1,1)/nr
       istep=istep+1
      enddo

      return
      end subroutine
