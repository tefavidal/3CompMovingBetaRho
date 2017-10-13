      program main
      
      implicit none


      INTEGER, PARAMETER :: Nx=200
      INTEGER, PARAMETER :: Ny=200

!           ZERO PILLAR 35 PERCENT
!      INTEGER, PARAMETER :: Nc=14000

!           ONE PILLAR 40 PERCENT
!      INTEGER, PARAMETER :: Nc=15878
!           ONE PILLAR 35 PERCENT
!      INTEGER, PARAMETER :: Nc=13893

!           ONE PILLAR 30 PERCENT
!      INTEGER, PARAMETER :: Nc=11908

!         FOUR PILLARS   40 PERCENT
!      INTEGER, PARAMETER :: Nc=15512

!         FOUR PILLARS   35 PERCENT
      INTEGER, PARAMETER :: Nc=13573

!         FOUR PILLARS   30 PERCENT
!      INTEGER, PARAMETER :: Nc=11634

!     Nc decided as (Nx*Ny-Pillars)*percentageofcells

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision beta(Nc),gamma(Nx,Ny),ro(Nc)
      double precision t
      double precision cells(Nc,2)
      integer pgbeg, i, j, counter, factor
      integer grid(Nx,Ny)
      character(len=4) ct1
      character(len=54) ct2
      character(len=58) ct3

      double precision gamma0(10),ro0(10),beta0(10)

      t=0.d0
      counter=0;
      call anfang(t,Nx,Ny,Nc,beta,gamma,ro,cells)
!      ct2='/data.lfpn/evidal/OutputData2D/data'
      ct2='/data.lfpn/evidal/3CompMovingBetaRho/OutputData2D/data'
      write(ct1,'(I4)') counter
      ct3 = ct2 // ct1
          open(10,file =ct3
     .       ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,Nc,gamma,ro,beta,cells)
          close(10)
      

 5    continue


      call ODE(t,Nx,Ny,Nc,beta,gamma,ro,cells)
      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      write(ct1,'(I4)') counter
      ct3 = ct2 // ct1

      write(6,*) ct3

          open(10,file =ct3
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,Nc,gamma,ro,beta,cells)
          close(10)

!      if(counter .eq. 65)then
!          write(6,*) 'Cutting the spiral'
!      if(mod(counter,70) .eq. 0)then
!         write(6,*) 'Adding Perturbation'
!         call FromCellToGrid(Nx,Ny,Nc,cells,grid)
!         do j=1,Ny
!            do i=1,Nx
!!          do j=(Ny/4-5),(Ny/4+5)
!!              do i=(Nx/4-5),(Nx/4+5)
!!               if(((i-(Nx/4))*(i-(Nx/4))+(j-(Ny/4))*(j-(Ny/4)))
!!     .          .le. 25)then
!!               if(((i-(Nx/2))*(i-(Nx/2))+(j-(Ny/2))*(j-(Ny/2)))
!!     .          .gt. (90*90))then
!
!               if(j .lt. 10)then
!
!!          do j=(Ny/2-5),(Ny/2+5)
!!              do i=(Nx/2-5),(Nx/2+5)
!!         if(((i-(Nx/2))*(i-(Nx/2))+(j-(Ny/2))*(j-(Ny/2))) .le. 25)then
!
!
!                  if(grid(i,j) .gt. 0.5)then
!                     factor=1.0
!                  else
!                     factor=0.0
!                  endif
!                  gamma(i,j)=factor*(gamma01+4)
!               endif
!            enddo
!         enddo
!      endif



      if (t+dt .lt. tend) then
         goto 5

      endif

      close(10)

!

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
!      open(42,file ='OutputData2D/Final-State'
!     . ,status = 'unknown',form = 'formatted')
!
!      open(43,file ='OutputData2D/Final-Positions'
!     . ,status = 'unknown',form = 'formatted')
!
!      call outFinal(t,Nx,Ny,Nc,gamma,ro,beta,cells)
!      close(42)
!      close(43)

      end
      


