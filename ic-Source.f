      subroutine ic(t,Nx,Ny,Nc,beta,gamma,ro,cells)
      
      implicit none

      double precision t
      integer Nx, Ny,Nc

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision beta(Nc),gamma(Nx,Ny),ro(Nc)
      double precision gamma0(10),beta0(10),ro0(10)
      double precision dke(Nx,Ny),dsigma(Nx,Ny)
      double precision cells(Nc,2)
      integer nfix, i, j
      integer grid(Nx,Ny)
      integer factor

!     %%%%% Initial State for ke=9.5 sigma=0.55
!         nfix=1
!         gamma0(1)=0.006
!         beta0(1)=0.3167
!         ro0(1)=0.8953

!     %%%%% Initial State for ke=7.0 sigma=0.55
         nfix=1
         gamma0(1)=0.359
         beta0(1)=13.91
         ro0(1)=0.285
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do i=1,nfix
      write(6,*)  'gamma=',gamma0(i),'  beta=',beta0(i),'  ro=',ro0(i)
      enddo

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)

      do i=1,Nx
       do j=1,Ny
           if(grid(i,j) .gt. 0.5)then
            gamma(i,j)=gamma0(1)
            else
            gamma(i,j)=0.0
            endif
       enddo
      enddo

      do i=1,Nc
        if(cells(i,1) .lt. 0 .and. cells(i,2) .lt. 0)then
         ro(i)=0.0
         beta(i)=0.0
        endif
         ro(i)=ro0(1)
         beta(i)=beta0(1)
      enddo

!%%%%%%%%%%%%%%%%%Adding perturbation at the center
!      do j=1,Ny
!         do i=1,Nx
!!            if(((i-(Nx/4))*(i-(Nx/4))+(j-(Ny/4))*(j-(Ny/4)))
!!     .       .le. 25)then
!
!            if(((i-(2*Nx/5))*(i-(2*Nx/5))+(j-(2*Ny/5))*(j-(2*Ny/5)))
!     .       .le. 25)then
!
!!            if(((i-(Nx/2))*(i-(Nx/2))+(j-(Ny/2))*(j-(Ny/2)))
!!     .       .le. 25)then
!
!
!               if(grid(i,j) .gt. 0.5)then
!                  factor=1.0
!               else
!                  factor=0.0
!               endif
!               gamma(i,j)=factor*(gamma0(1)+3)
!            endif
!         enddo
!      enddo



         gamma01=gamma0(1)
         beta01=beta0(1)
         ro01=ro0(1)

      return

      end

