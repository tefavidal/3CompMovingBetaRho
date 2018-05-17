      subroutine rs(t,Nx,Ny,Nc,beta,gamma,ro,betaprime,gammaprime,
     . roprime,cells,vdx)

      implicit none
      double precision t, YesCell, aux, obstac
      integer Nx, Ny, i, j,ID, Nc
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision beta(Nc),gamma(Nx,Ny),ro(Nc)
      double precision betagrid(Nx,Ny), rhogrid(Nx,Ny)
      double precision betaprime(Nc),gammaprime(Nx,Ny),roprime(Nc)
      double precision f1,f2,Phi,Y,noise,r1,r2, keB, keU
      double precision gLaplace(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny),vdy
      double precision gammalist(Nc)
      double precision cells(Nc,2)
      integer grid(Nx,Ny)
  ! ----- variables for portable seed setting -----
!      INTEGER :: i_seed
!      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
!      INTEGER, DIMENSION(1:8) :: dt_seed
!  ! ----- end of variables for seed setting -----
!
! ! ----- Set up random seed portably -----
!      CALL RANDOM_SEED(size=i_seed)
!      ALLOCATE(a_seed(1:i_seed))
!      CALL RANDOM_SEED(get=a_seed)
!      CALL DATE_AND_TIME(values=dt_seed)
!      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
!     .*dt_seed(6)
!!       write(6,*) 'seed=',a_seed(i_seed)
!      CALL RANDOM_SEED(put=a_seed)
!      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----


      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      call GridifyingBetaRho(Nx,Ny,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)
      call GammaToList(Nx,Ny,Nc,cells,gamma,gammalist)


          keB=5.0
          keU=2.0
      do j=1,Ny
       do i=1,Nx
!       Extra variables calculation
        if(grid(i,j) .gt. 0.5)then
            YesCell=1.0
        else
            YesCell=0.0
        endif

        if(grid(i,j) .lt. -0.5)then
            obstac=1.0
        else
            obstac=0.0
        endif

        vdy=0.0
!       Right hand side
!        call random_number(r1)
!        call random_number(r2)
!        noise=sqrt(-2*log(r1))*cos(2*Pi*r2);
        gammaprime(i,j)=YesCell*s2*betagrid(i,j)/depsilon
     .     -(YesCell*keB*grid(i,j)+keU)/dke0*gamma(i,j)/depsilon
     .                  +depsilon*gLaplace(i,j)
!     .                  +abs(noise)*0.1
!     .          -  (vdx(i,j)*xgradeC(i,j))

       enddo
      enddo

      do i=1,Nc
         aux=gammalist(i)
         f1=(1.d0+dk*aux)/(1.d0+aux)
         f2=(dL1+dk*dL2*dc*aux)/(1.d0+dc*aux)
         Y=ro(i)*aux/(1.d0+aux)
         Phi=(dlambda1+Y*Y)/(dlambda2+Y*Y)
         betaprime(i)=(s1*Phi-beta(i))/depsilonp
         roprime(i)=(-f1*ro(i)+f2*(1.d0-ro(i)))
      enddo

      return

      end
!      ***********************************************************
!     *************************************************

