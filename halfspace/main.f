      program main
      implicit none
      integer,parameter:: nx=128,ny=32,nz=128,
     +                    maxstepU=1000, 
     +                    lgx=int(dlog(nx)/dlog(2.0d0)),
     +                    lgy=int(dlog(nx)/dlog(2.0d0)),
     +                    lgz=int(dlog(nx)/dlog(2.0d0)),
     + n=nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+(nx-1)*(ny-1)*nz
      real*8,parameter :: pi=4.0d0*datan(1.0d0),
     +                    m0=4.0d0*pi*1.0d-7, 
     +                    tolU=1.0d-16,
     +                    sigmaAIR=1.0d-8
      character(64):: date,Scen,wpath,iargc
      integer:: i,j,ik,ix,iy,iz,num_threads,MaxThreads,
     +          omp_get_max_threads,omp_get_num_threads,
     +          Lx,Ly,Lz,NT,f
      real*8:: SIGMAEX(nx*(ny-1)*(nz-1)),SIGMADiffEX(nx*(ny-1)*(nz-1)),
     + SIGMAEY((nx-1)*ny*(nz-1)),SIGMADiffEY((nx-1)*ny*(nz-1)),
     + SIGMAEZ((nx-1)*nz*(ny-1)),bct1,bct2,
     + zz,zzz,zcounter,dznrm2,dnrm2,tstart,tfinish,rhstime,rhsend,
     + realExBack(nx*(ny-1)*(nz-1)),imagExBack(nx*(ny-1)*(nz-1)),
     + realEyBack((nx-1)*ny*(nz-1)),imagEyBack((nx-1)*ny*(nz-1)),
     + realETX(nx*(ny-1)*(nz-1)),imagETX(nx*(ny-1)*(nz-1)),
     + realETY((nx-1)*ny*(nz-1)),imagETY((nx-1)*ny*(nz-1)),
     + MultVal1(5*lgx),MultVal2(5*lgy),MultVal3(5*lgz),
     + a1(nx),b1(nx),c1(nx),a2(ny),b2(ny),c2(ny),a3(nz),b3(nz),c3(nz),
     + PosTx,PosTy,PosTz,sigma0,sigma,gl,hx,hy,hz,omega
      complex*16:: ci, 
     + ExBack(nx*(ny-1)*(nz-1)),EX(nx*(ny-1)*(nz-1)),
     + ETX(nx*(ny-1)*(nz-1)),EtotalX(nx*(ny-1)*(nz-1)),
     + F1(nx*(ny-1)*(nz-1)),DiffEX(nx*(ny-1)*(nz-1)),
     + e1(nx*(ny-1)*(nz-1)),tempEX(nx*(ny-1)*(nz-1)),
     + EyBack((nx-1)*ny*(nz-1)),EY((nx-1)*ny*(nz-1)),
     + F2((nx-1)*ny*(nz-1)),e2((nx-1)*ny*(nz-1)),
     + ETY((nx-1)*ny*(nz-1)),EtotalY((nx-1)*ny*(nz-1)),
     + tempEY((nx-1)*ny*(nz-1)),DiffEY((nx-1)*ny*(nz-1)),
     + EZ((nx-1)*nz*(ny-1)),tempEZ((nx-1)*(ny-1)*nz),
     + e3((nx-1)*nz*(ny-1)),HtotalZ(nx*ny*(nz-1)),
     + pt1(nx*(ny-1)*(nz-1)),st1(nx*(ny-1)*(nz-1)),
     + pt2(nx*ny*(nz-1)),st2(nx*ny*(nz-1)),
     + pt3(nx*(ny-1)*nz),st3(nx*(ny-1)*nz),
     + pt4((nx-1)*ny*(nz-1)),st4((nx-1)*ny*(nz-1)),
     + pt5((nx-1)*ny*nz),st5((nx-1)*ny*nz),
     + pt6((nx-1)*(ny-1)*nz),st6((nx-1)*(ny-1)*nz),
     + pt7((nx-1)*(ny-1)*(nz-1)),st7((nx-1)*(ny-1)*(nz-1)),
     + RHS1((2*nx-1)*ny*(nz-1)),RHS2(nx*(2*ny-1)*(nz-1)), 
     + RHS3((2*nx-1)*(ny-1)*nz),RHS4((nx-1)*(2*ny-1)*nz),
     + RHS5((nx-1)*ny*(2*nz-1)),RHS6(nx*(ny-1)*(2*nz-1)),
     + RHS7((nx-1)*(ny-1)*(2*nz-1)),
     + temp1(nx*(ny-1)*(nz-1)),temp11(nx*(ny-1)*(nz-1)),
     + temp2((nx-1)*ny*(nz-1)),temp22((nx-1)*ny*(nz-1)),
     + temp3((nx-1)*(ny-1)*nz),temp33((nx-1)*(ny-1)*nz),
     + tempEX_1(nx*(ny-1)*(nz-1)),tempEX_2(nx*(ny-1)*(nz-1)),
     + tempEY_1((nx-1)*ny*(nz-1)),tempEY_2((nx-1)*ny*(nz-1)),
     + tempEZ_1((nx-1)*(ny-1)*nz),tempEZ_2((nx-1)*(ny-1)*nz),
     + tempEX_3(nx*(ny-1)*(nz-1)),tempEY_3((nx-1)*ny*(nz-1)),
     + tempEZ_3((nx-1)*(ny-1)*nz),zdotc,tempHZ(nx*ny*(nz-1)),
     + sv1(nx*(ny-1)),sv2((nx-1)*ny),
     + bcg1(n),bcg2(n),bcg3(n),bcg4(n),bcg5(n),bcg6(n),bcg7(n),
     + bcg8(n),E(n)

      call date_and_time(date)
      NT=omp_get_max_threads()
      call omp_set_num_threads(NT)
      
      wpath=iargc()
      call getarg(1,wpath)
      call params(nx,ny,nz,Lx,Ly,Lz,hx,hy,hz,f,omega,pi,
     +            sigma0,sigma,sigmaAIR,
     +            PosTx,PosTy,PosTz,gl,Scen,tolU,maxstepU,
     +            wpath)        
      call cpu_time(tstart)
      ci=dcmplx(0.0d0,1.0d0)
      call makeS(nx,ny,nz,hz,gl,
     +           SIGMAEX,SIGMAEY,SIGMAEZ,sigmaAIR,sigma,sigma0,
     +           SIGMADiffEX,SIGMADiffEY)

      call makeRHS(nx,ny,nz,omega,m0,ci,wpath,
     +        SIGMAEX,SIGMAEY,SIGMAEZ,SIGMADiffEX,SIGMADiffEY,
     +        realExBack,imagExBack,realEyBack,imagEyBack,ExBack,EyBack,
     +        e1,e2,e3,F1,F2)

      call cpu_time(rhsend)
      print*,'ExBack norm',
     +         dznrm2(size(ExBack),ExBack,1)/(nx*(ny-1)*(nz-1))
      print*,'EyBack norm',
     +         dznrm2(size(EyBack),EyBack,1)/((nx-1)*ny*(nz-1))
      print*, 'NORM F1',dznrm2(nx*(ny-1)*(nz-1),F1,1)/(nx*(ny-1)*(nz-1))
      print*, 'NORM F2',dznrm2(ny*(nx-1)*(nz-1),F2,1)/(ny*(nx-1)*(nz-1))
      print*, 'NORM e1',dznrm2(nx*(ny-1)*(nz-1),e1,1)/(nx*(ny-1)*(nz-1))
      print*, 'NORM e2',dznrm2(ny*(nx-1)*(nz-1),e2,1)/(ny*(nx-1)*(nz-1))
      print*, 'NORM e3',dznrm2(nz*(nx-1)*(ny-1),e3,1)/(nz*(nx-1)*(ny-1))
      print*, 'NORM SIGMAEX',dnrm2(nx*(ny-1)*(nz-1),SIGMAEX,1)/
     +                       (nx*(ny-1)*(nz-1))
      print*, 'NORM SIGMAEY',dnrm2(ny*(nx-1)*(nz-1),SIGMAEY,1)/
     +                       (ny*(nx-1)*(nz-1))
      print*, 'NORM SIGMAEZ',dnrm2(nz*(nx-1)*(ny-1),SIGMAEZ,1)/
     +                       (nz*(nx-1)*(ny-1))
      
      call cpu_time(bct1)

      call bicgstabU(n,nx,ny,nz,hx,hy,hz,ci,
     +     e1,e2,e3,F1,F2,
     +     temp1,temp11,temp2,temp22,temp3,temp33,
     +     pt1,pt2,pt3,pt4,pt5,pt6,pt7,
     +     st1,st2,st3,st4,st5,st6,st7,  
     +     RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,
     +     MultVal1,MultVal2,MultVal3,
     +     a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +     maxstepU,tolU,
     +     tempEX,tempEY,tempEZ,
     +     tempEX_1,tempEY_1,tempEZ_1,
     +     tempEX_2,tempEY_2,tempEZ_2,
     +     tempEX_3,tempEY_3,tempEZ_3,
     +     bcg1,bcg2,bcg3,bcg4,bcg5,bcg6,bcg7,bcg8,
     +     E)
      call cpu_time(bct2)


!$OMP PARALLEL
!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       EX(i)=E(i)
      enddo
!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       EY(i)=E(nx*(ny-1)*(nz-1)+i)
      enddo
!$OMP DO
      do i=1,(nx-1)*(ny-1)*nz
       EZ(i)=E(nx*(ny-1)*(nz-1)+(nx-1)*ny*(nz-1)+i)
      enddo
!$OMP END PARALLEL


!$OMP PARALLEL
!$OMP DO
      do i=1,nx*(ny-1)*(nz-1)
       EtotalX(i)=EX(i)+ExBack(i)
      enddo
!$OMP DO
      do i=1,(nx-1)*ny*(nz-1)
       EtotalY(i)=EY(i)+EyBack(i)
      enddo
!$OMP END PARALLEL

      print*,'EX norm',dznrm2(size(EX),EX,1)/(nx*(ny-1)*(nz-1))
      print*,'EY norm',dznrm2(size(EY),EY,1)/((nx-1)*ny*(nz-1))
      print*,'EZ norm',dznrm2(size(EZ),EZ,1)/((nx-1)*(ny-1)*nz)
      print*,'EtotalX norm',
     +         dznrm2(size(EtotalX),EtotalX,1)/(nx*(ny-1)*(nz-1))
      print*,'EtotalY norm',
     +         dznrm2(size(EtotalY),EtotalY,1)/((nx-1)*ny*(nz-1))
c     print*,'DiffEX norm',
c    +         dznrm2(size(DiffEX),DiffEX,1)/(nx*(ny-1)*(nz-1))
c     print*,'DiffEY norm',
c    +         dznrm2(size(DiffEY),DiffEY,1)/((nx-1)*ny*(nz-1))
c     print*,'ETX norm',dznrm2(size(ETX),ETX,1)/(nx*(ny-1)*(nz-1))
c     print*,'ETY norm',dznrm2(size(ETY),ETY,1)/((nx-1)*ny*(nz-1))
       
      call makeH(EtotalX,EtotalY,nx,ny,nz,hx,hy,
     +           m0,omega,ci,tempHZ,
     +           pt1,pt2,pt3,pt4,pt5,pt6,pt7,
     +           st1,st2,st3,st4,st5,st6,st7, 
     +           RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,
     +           MultVal1,MultVal2,MultVal3,
     +           a1,a2,a3,b1,b2,b3,c1,c2,c3,sv1,sv2,
     +           HtotalZ)
      print*,'makeH COMPLETED'
      print*,'HtotalZ norm',
     +         dznrm2(size(HtotalZ),HtotalZ,1)/(nx*ny*(nz-1))

      call cpu_time(tfinish)

       print*,'Time =',tfinish-tstart,'seconds'
       open(12, file=trim(wpath)//'temp/EX.txt')
       open(13, file=trim(wpath)//'temp/EY.txt')
       open(14, file=trim(wpath)//'temp/EZ.txt')
       open(15, file=trim(wpath)//'temp/HZ.txt')
       open(16, file=trim(wpath)//'temp/EtotalX.txt')
       open(17, file=trim(wpath)//'temp/EtotalY.txt')
       open(18, file=trim(wpath)//'temp/info.txt')
       print*,'Done writing in', trim(wpath)

       do i=1,nx*(ny-1)*(nz-1)
        write(12,*) EX(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(13,*) EY(i)
       enddo
       do i=1,nz*(ny-1)*(nx-1)
        write(14,*) EZ(i)
       enddo
       do i=1,nx*ny*(nz-1)
        write(15,*) HtotalZ(i)
       enddo
       do i=1,nx*(ny-1)*(nz-1)
        write(16,*) EtotalX(i)
       enddo
       do i=1,(nx-1)*ny*(nz-1)
        write(17,*) EtotalY(i)
       enddo
       write(18,*) nx,ny,nz
       write(18,*) Lx,Ly,Lz
       write(18,*) PosTx,PosTy,PosTz
       write(18,*) gl
       write(18,*) sigma
       write(18,*) sigma0
       write(18,*) f 
       write(18,*) sigmaAIR
       write(18,*) tolU 
       write(18,*) maxstepU 

       close(12) 
       close(13) 
       close(14) 
       close(15) 
       close(16) 
       close(17) 
       close(18) 
       end


