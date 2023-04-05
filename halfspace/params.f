      subroutine params(nx,ny,nz,Lx,Ly,Lz,hx,hy,hz,f,omega,pi,
     +                 sigma0,sigma,sigmaAIR,
     +                 PosTx,PosTy,PosTz,gl,Scen,tolU,maxstepU,
     +                 wpath)

       integer:: f,nx,ny,nz,Lx,Ly,Lz,maxstepU
       real*8:: PosTx,PosTy,PosTz,gl,sigma0,sigma,sigmaAIR,tolU,omega,
     +          hx,hy,hz,pi
       character(64):: wpath

! Subroutine params.f parses the problem physical parameters from config.dat

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
