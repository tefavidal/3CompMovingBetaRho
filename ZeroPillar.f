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
      double precision d
      integer grid(Nx,Ny)



      do j= 1,Ny
        do i= 1, Nx
                grid(i,j)=0
        enddo
      enddo

      return
      end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      subroutine initialCenterIsland(Nx,Ny,Nc,cells,d)

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
      double precision cells(Nc,2), percentage, d, r, r0
      real aux
      integer grid(Nx,Ny), auxInt, x0,y0


      r=d/2.0*dk1/(dke0*Diffgamma)**0.5
      r=r*r

      x0=ceiling(Nx/2.0)
      y0=ceiling(Ny/2.0)


      k=1

      do j=1,Ny
        do i=1,Nx
            r0=(i-x0)*(i-x0)*dx*dx + (j-y0)*(j-y0)*dy*dy
            if (r0 .le. r)then
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
        cells(i,1)=-1.0
        cells(i,2)=-1.0
      enddo

      call FromCellToGrid(Nx,Ny,Nc,cells,grid)




      return
      end
