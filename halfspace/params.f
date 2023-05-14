      subroutine params(nx,ny,nz,Lx,Ly,Lz,hx,hy,hz,f,omega,pi,
     +                 sigma0,sigma,sigmaAIR,
     +                 PosTx,PosTy,PosTz,gl,tolU,maxstepU,
     +                 wpath)

       integer:: f,nx,ny,nz,Lx,Ly,Lz,maxstepU
       real*8:: PosTx,PosTy,PosTz,gl,sigma0,sigma,sigmaAIR,tolU,omega,
     +          hx,hy,hz,pi
       character(64):: wpath

! ===================================================================
! Title: params.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! params.f parses the models physical parameters from config.dat and
! returns them to main.f
!====================================================================      
! Input:
! nx: number of elements in x-direction
! ny: number of elements in y-direction
! nz: number of elements in z-direction
! Output:
! Lx: length of edge in x-direction
! Ly: length of edge in y-direction
! Lz: length of edge in z-direction
! hx: discretization step in x-direction
! hy: discretization step  in y-direction
! hz: discretization step  in z-direction
! sigma: conductivity 
! sigma0: background conductivity
! sigmaAIR: conductivity of the air
! f: frequency
! omega: angular frequency
!==================================================================== 

        open(1,file=trim(wpath)//"config.dat")
        read(1,*) Lx,Ly,Lz
        read(1,*) PosTx,PosTy,PosTz
        read(1,*) gl
        read(1,*) sigma0
        read(1,*) sigma
        read(1,*) f
        print*,'Problem parameters parsed'
        close(1)

        hx=Lx/nx
        hy=Ly/ny
        hz=Lz/nz
        omega=2*pi*f

        print*,'nx',nx,'ny',ny,'nz',nz
        print*,'PosTx',PosTx,'PosTy',PosTy,'PosTz',PosTz
        print*,'Lx',Lx,'Ly',Ly,'Lz',Lz
        print*,'Ground Line',gl,'sigmaAIR',sigmaAIR
        print*,'f',f,'sigma',sigma,'sigmaO',sigma0
        print*,'tolU',tolU, 'maxU',maxstepU
        end subroutine
