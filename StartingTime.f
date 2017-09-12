      subroutine StartingTime(t,Nx,Ny,TS)
      
      implicit none
      integer Nx, Ny
      double precision t, TimeGap
      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision TS(Nx,Ny)
      integer ClumpX, ClumpY, i, j
      real aux, aux2
  ! ----- variables for portable seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

 ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
     .*dt_seed(6)
       write(6,*) 'seed=',a_seed(i_seed)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----


        ClumpX=ceiling(0.1/(dx/dk1*(dke0*Diffgamma)**0.5))
        ClumpY=ceiling(0.1/(dy/dk1*(dke0*Diffgamma)**0.5))
        TimeGap=100.0

      write(6,*) 'Clump x=',ClumpX
      write(6,*) 'Clump y=',ClumpY
      write(6,*) 'Time Gap=',TimeGap

      do i=1,Nx
        if (mod(i-1,ClumpX) .eq. 0) then
            do j=1,Ny
                if (mod(j-1,ClumpY) .eq. 0) then
                    call random_number(aux)
                    call random_number(aux2)
                endif

!      %%%%%%%%%% Exponential distribution
                if(dke0 .eq. 2.5)then
                        TS(i,j)=-25*log(aux)+TimeGap
                else
                         TS(i,j)=-100*log(aux)+TimeGap
               endif
            enddo

        else
            do j=1,Ny
                TS(i,j)=TS(i-1,j)
            enddo
        endif
      enddo

      return
      end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FromCellToGrid(Nx,Ny,Nc,cells,grid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx,Ny, Nc,i, j,k
      double precision cells(Nc,2)
      double precision d
      integer grid(Nx,Ny)

!      do j=1,Ny
!        do i=1,Nx
!            grid(i,j)=0
!        enddo
!      enddo
        d=1.0
       call Pillars(Nx,Ny,d,grid)

      do k=1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        if(i .gt.Nx)then
            cells(k,1)=cells(k,1)-Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(i .lt.1)then
            cells(k,1)=cells(k,1)+Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(j .gt. Ny)then
            cells(k,2)=cells(k,2)-Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        if(j .lt. 1)then
            cells(k,2)=cells(k,2)+Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        grid(i,j)=grid(i,j)+1
      enddo

      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine GridifyingBetaRho(Nx,Ny,Nc,cells,grid,
     . beta,ro,betagrid,rhogrid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx,Ny, Nc,i, j,k
      double precision beta(Nc),ro(Nc), cells(Nc,2)
      double precision betagrid(Nx,Ny) , rhogrid(Nx,Ny)
      double precision d
      integer grid(Nx,Ny)

      do j=1,Ny
        do i=1,Nx
            betagrid(i,j)=0.0
            rhogrid(i,j)=0.0
        enddo
      enddo

        d=1.0
       call Pillars(Nx,Ny,d,grid)

      do k=1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        if(i .gt.Nx)then
            cells(k,1)=cells(k,1)-Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(i .lt.1)then
            cells(k,1)=cells(k,1)+Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(j .gt. Ny)then
            cells(k,2)=cells(k,2)-Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        if(j .lt. 1)then
            cells(k,2)=cells(k,2)+Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        betagrid(i,j)=beta(k)
        grid(i,j)=grid(i,j)+1
        rhogrid(i,j)=ro(k)
      enddo

      return
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine GammaToList(Nx,Ny,Nc,cells,gamma,gammalist)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx,Ny, Nc,i, j,k
      double precision cells(Nc,2)
      double precision d
      double precision gamma(Nx,Ny), gammalist(Nc)




      do k=1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        if(i .gt.Nx)then
            cells(k,1)=cells(k,1)-Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(i .lt.1)then
            cells(k,1)=cells(k,1)+Nx*dx
            i=ceiling(cells(k,1)/dx)
        endif

        if(j .gt. Ny)then
            cells(k,2)=cells(k,2)-Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        if(j .lt. 1)then
            cells(k,2)=cells(k,2)+Ny*dy
            j=ceiling(cells(k,2)/dy)
        endif

        gammalist(k)=gamma(i,j)
      enddo

      return
      end


!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine initialDiscreteDistribution(Nx,Ny,Nc,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      integer Nx,Ny, Nc,i, j,k
      double precision cells(Nc,2), percentage, d
      real aux
      integer grid(Nx,Ny), auxInt

  ! ----- variables for portable seed setting -----
      INTEGER :: i_seed
      INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
      INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----

 ! ----- Set up random seed portably -----
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)
     .*dt_seed(6)
       write(6,*) 'seed=',a_seed(i_seed)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
  ! ----- Done setting up random seed ----

        d=1.0;
       call Pillars(Nx,Ny,d,grid)
            k=1;
            percentage=real(Nc)/(real(Nx)*real(Ny))
      do i=1,Nx
        do j=1,Ny
            call random_number(aux)
            if(aux .lt. percentage .and. grid(i,j) .ge. 0.0)then
                if (k .gt. Nc)then
                    exit
                endif
                cells(k,1)=(i-0.5)*dx
                cells(k,2)=(j-0.5)*dy
                k=k+1
            endif
        enddo
      enddo

      auxInt=k-1

      do i=k,Nc
        cells(i,1)=cells(1,1)
        cells(i,2)=cells(1,2)
      enddo

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)


      do while (k .le. Nc)
        call random_number(aux)
        i=ceiling(aux*Nx)
        call random_number(aux)
        j=ceiling(aux*Ny)

        if(grid(i,j) .eq. 0)then
            cells(k,1)=(i-0.2)*dx
            cells(k,2)=(j-0.2)*dy
            k=k+1
            grid(i,j)=1
        endif


      enddo




      return
      end

!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
