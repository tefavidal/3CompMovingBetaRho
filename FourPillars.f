      subroutine Pillars(Nx,Ny,d,grid)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      integer Nx,Ny,i,j

      double precision d,r,r0,r1,r2,r3
      integer grid(Nx,Ny)
      integer x0, y0, x1,y1

!      d=1.0
      r=1.0/2.0*dk1/(dke0*Diffgamma)**0.5
      r=r*r

      x0=ceiling(Nx/4.0)
      y0=ceiling(Ny/4.0)

      x1=ceiling(3.0*Nx/4.0)
      y1=ceiling(3.0*Ny/4.0)

      do j= 1,Ny
        do i= 1, Nx
            r0=(i-x0)*(i-x0)*dx*dx + (j-y0)*(j-y0)*dy*dy
            r1=(i-x1)*(i-x1)*dx*dx + (j-y1)*(j-y1)*dy*dy
            r2=(i-x0)*(i-x0)*dx*dx + (j-y1)*(j-y1)*dy*dy
            r3=(i-x1)*(i-x1)*dx*dx + (j-y0)*(j-y0)*dy*dy
       if (r0 .le. r .or. r1 .le. r .or.r2 .le. r .or. r3 .le. r) then
                grid(i,j)=-1
            else
                grid(i,j)=0
                endif
        enddo
      enddo

      return
      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


      subroutine initialHigherAroundPillarDistribution(Nx,Ny,Nc,cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      integer Nx,Ny, Nc,i, j,k
      double precision cells(Nc,2), percentage, d, r, r0,r1,r2,r3
      real aux
      integer grid(Nx,Ny), auxInt, x0,y0,x1,y1

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

      r=1.5/2.0*dk1/(dke0*Diffgamma)**0.5
      r=r*r

      x0=ceiling(Nx/4.0)
      y0=ceiling(Ny/4.0)

      x1=ceiling(3.0*Nx/4.0)
      y1=ceiling(3.0*Ny/4.0)



      do j=1,Ny
        do i=1,Nx
            call random_number(aux)
            r0=(i-x0)*(i-x0)*dx*dx + (j-y0)*(j-y0)*dy*dy
            r1=(i-x1)*(i-x1)*dx*dx + (j-y1)*(j-y1)*dy*dy
            r2=(i-x0)*(i-x0)*dx*dx + (j-y1)*(j-y1)*dy*dy
            r3=(i-x1)*(i-x1)*dx*dx + (j-y0)*(j-y0)*dy*dy
       if (r0 .le. r .or. r1 .le. r .or.r2 .le. r .or. r3 .le. r) then
                percentage=2.0*real(Nc)/(real(Nx)*real(Ny))
            else
                percentage=real(Nc)/(real(Nx)*real(Ny))
            endif
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

