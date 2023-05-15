      subroutine makeS(nx,ny,nz,hz,gl,
     +           SIGMAEX,SIGMAEY,SIGMAEZ,sigmaAIR,sigma,sigma0,
     +           SIGMADiffEX,SIGMADiffEY)

      implicit none
      integer:: nx,ny,nz,ik,ix,iy,iz
      real*8 :: hz,sigmaAIR,sigma0,sigma,gl,zcounter,
     +          SIGMAEX(nx*(ny-1)*(nz-1)),SIGMAEY((nx-1)*ny*(nz-1)),
     +          SIGMAEZ((nx-1)*(ny-1)*nz),SIGMADiffEX(nx*(ny-1)*(nz-1)),
     +          SIGMADiffEY((nx-1)*ny*(nz-1)),dnrm2


! ===================================================================
! Title: makeS.f 
! Authors: N. Vilanakis, E. Mathioudakis
! Details: Applied Mathematics and Computers Lab, Technical University of Crete
!====================================================================
! makeS.f computes the necessary conductivity vector coefficients 
! σx (SIGMAEX), σy (SIGMAEY), σz (SIGMAEZ)
! depending on the values of σ (sigma) and σ0  (sigma0) and the ground level (gl)
! specified by the user in config.dat file
!====================================================================      
! Input:
! nx: integer, number of elements in x-direction
! ny: integer, number of elements in y-direction
! nz: integer, number of elements in z-direction
! gl: real, ground level in the cuboid
! sigma: real, conductivity 
! sigma0: real, background conductivity
! sigmaAIR: real, conductivity of the air
! Output:
! SIGMAEX: real array, conductivity σ in Ex nodes, dimension: nx*(ny-1)*(nz-1)
! SIGMAEY: real array, conductivity σ in Ey nodes, dimension: (nx-1)*ny*(nz-1) 
! SIGMAEZ: real array, conductivity σ in Ez nodes, dimension: (nx-1)*(ny-1)*nz
! SIGMADiffEX: real array, Difference between conductivity σ and background conductivity σ0 in Ex nodes, dimension: nx*(ny-1)*(nz-1)
! SIGMADiffEY: real array, Difference between conductivity σ and background conductivity σ0 in Ey nodes, dimension: (nx-1)*ny*(nz-1) 
!==================================================================== 
c SIGMA-EX
!$OMP PARALLEL DO PRIVATE(ik,zcounter) SHARED(SIGMAEX,SIGMADiffEX)
      do iz=1,nz-1
       zcounter=iz*hz
       do iy=1,ny-1
        do ix=1,nx
         ik=(iz-1)*(ny-1)*nx+(iy-1)*nx+ix
c     Below Ground
         if (zcounter.le.gl) then
          SIGMADiffEX(ik)=sigma-sigma0
          SIGMAEX(ik)=sigma
         elseif (zcounter.gt.gl) then
c     Above Ground
          SIGMADiffEX(ik)=sigmaAIR
          SIGMAEX(ik)=sigmaAIR
         endif
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

c SIGMA-EY
!$OMP PARALLEL DO PRIVATE(ik,zcounter) SHARED(SIGMAEY,SIGMADiffEY) 
      do iz=1,nz-1
       zcounter=iz*hz
       do iy=1,ny
        do ix=1,nx-1
         ik=(iz-1)*ny*(nx-1)+(iy-1)*(nx-1)+ix
c     Below Ground
         if (zcounter.le.gl) then
          SIGMADiffEY(ik)=sigma-sigma0
          SIGMAEY(ik)=sigma
c     Above Ground
         else
          SIGMADiffEY(ik)=sigmaAIR
          SIGMAEY(ik)=sigmaAIR
         endif
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO 

c SIGMA-EZ
!$OMP PARALLEL DO PRIVATE(ik,zcounter) SHARED(SIGMAEZ)
      do iz=1,nz
       zcounter=iz*hz
       do iy=1,ny-1
        do ix=1,nx-1
         ik=(iz-1)*(ny-1)*(nx-1)+(iy-1)*(nx-1)+ix
         if (zcounter.le.gl) then
          SIGMAEZ(ik)=sigma
         else
          SIGMAEZ(ik)=sigmaAIR
         endif
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
