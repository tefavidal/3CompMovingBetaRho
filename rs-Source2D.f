      subroutine rs(t,Nx,Ny,Nc,beta,gamma,ro,betaprime,gammaprime,
     . roprime,cells,vdx)

      implicit none
      double precision t, factor, aux
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
      double precision f1,f2,Phi,Y
      double precision gLaplace(Nx,Ny)
      double precision xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision vdx(Nx,Ny),vdy
      double precision gammalist(Nc)
      double precision cells(Nc,2)
      integer grid(Nx,Ny), factor2



      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      call GridifyingBetaRho(Nx,Ny,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)
      call GammaToList(Nx,Ny,Nc,cells,gamma,gammalist)

      do j=1,Ny
       do i=1,Nx
!       Extra variables calculation
        if(grid(i,j) .gt. 0.5)then
            factor=1.0
        else
            factor=0.0
        endif
        vdy=0.0
!       Right hand side
        gammaprime(i,j)=factor/depsilon*s2*betagrid(i,j)
     .             -1.0*gamma(i,j)/depsilon
     .                  +depsilon*gLaplace(i,j)
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
      subroutine functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      implicit none
      integer Nx, Ny,i,j

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision gamma(Nx,Ny)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision gammaim2,gammaim1,gammaip1,gammajm1,gammajp1
      double precision gLapX,gLapY,thetai,thetaim1,psii,psiim1, d
      integer grid(Nx,Ny)

       d=1.0
      call Pillars(Nx,Ny,d,grid)

      do j=1,Ny
       do i=1,Nx
!       No-Flux boundary condition
!       if(i .eq. 1) then
       if(i .eq. 1 .or. grid(i-1,j) .le. -0.5) then
        gammaim2=gamma(i+1,j)
        gammaim1=gamma(i+1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0
!        gammaim1=0

       elseif(i .eq. 2) then
        gammaim2=-gamma(i,j)+2*gamma(i-1,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

!        gammaim2=0

!       elseif(i .eq. Nx) then
       elseif(i .eq. Nx .or. grid(i+1,j) .le. -0.5) then
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i-1,j)
       else
        gammaim2=gamma(i-2,j)
        gammaim1=gamma(i-1,j)
        gammaip1=gamma(i+1,j)

       endif

!      Periodic Boundary
!       if(i .eq. 1) then
!        gammaim2=gamma(Nx-1,j)
!        gammaim1=gamma(Nx,j)
!        gammaip1=gamma(2,j)
!
!       elseif(i .eq. 2) then
!        gammaim2=gamma(Nx,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(i+1,j)
!
!
!       elseif(i .eq. Nx) then
!        gammaim2=gamma(i-2,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(1,j)
!       else
!        gammaim2=gamma(i-2,j)
!        gammaim1=gamma(i-1,j)
!        gammaip1=gamma(i+1,j)
!
!       endif


!     Non Flux Boundary Cond
       if(Ny .eq. 1) then
        gammajm1=gamma(i,j)
        gammajp1=gamma(i,j)
!       elseif(j .eq. 1) then
       elseif(j .eq. 1 .or. grid(i,j-1) .le. -0.5) then
        gammajm1=gamma(i,j+1)
        gammajp1=gamma(i,j+1)
!       elseif(j .eq. Ny) then
       elseif(j .eq. Ny .or. grid(i,j+1) .le. -0.5) then
        gammajp1=gamma(i,j-1)
        gammajm1=gamma(i,j-1)
       else
        gammajp1=gamma(i,j+1)
        gammajm1=gamma(i,j-1)
       endif

!     Periodic Boundary
!
!       if(Ny .eq. 1) then
!        gammajm1=gamma(i,j)
!        gammajp1=gamma(i,j)
!       elseif(j .eq. 1) then
!        gammajm1=gamma(i,Ny)
!        gammajp1=gamma(i,j+1)
!       elseif(j .eq. Ny) then
!        gammajp1=gamma(i,1)
!        gammajm1=gamma(i,j-1)
!       else
!        gammajp1=gamma(i,j+1)
!        gammajm1=gamma(i,j-1)
!       endif


        gLapX=(gammaip1+gammaim1-2*gamma(i,j))/(dx**2)
        gLapY=(gammajp1+gammajm1-2*gamma(i,j))/(dy**2)
        gLaplace(i,j)=gLapX+gLapY


!
!        if(gammaip1 .eq. gamma(i,j)) then
!        thetai=1.d-10
!        else
!        thetai=(gamma(i,j)-gammaim1)/(gammaip1-gamma(i,j))+1.d-10
!        endif
!
!        if(gamma(i,j) .eq. gammaim1)then
!        thetaim1=1.d-10
!        else
!        thetaim1=(gammaim1-gammaim2)/(gamma(i,j)-gammaim1)+1.d-10
!        endif
!
!        psii=max(0.0,min(1.0,1.0/3.0+thetai/6.0,thetai))
!        psiim1=max(0.0,min(1.0,1.0/3.0+thetaim1/6.0,thetaim1))
!
!      xgradeC(i,j)=(1.0-psiim1+psii/thetai)*(-gammaim1+gamma(i,j))/(dx)


        ygradeC(i,j)=(gammajp1-gammajm1)/(2*dy)
        xgradeC(i,j)=(gammaip1-gammaim1)/(2*dx)
!        ygradeC(i,j)=(gammajp1-gamma(i,j))/(dy)
       enddo
      enddo

      return
      end



!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine UpdateDiscreteCells(t,h,Nx,Ny,Nc, gamma, ro,beta,
     . cells)

      implicit none

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      integer Nx, Ny,Nc, i, j, k, oldi, oldj
      double precision cells(Nc,2), t
      double precision gamma(Nx,Ny), rhogrid(Nx,Ny), betagrid(Nx,Ny)
      double precision beta(Nc), ro(Nc)
      double precision gLaplace(Nx,Ny),xgradeC(Nx,Ny),ygradeC(Nx,Ny)
      double precision grad, angle, h,speed, newx, newy
      integer grid(Nx,Ny)


      call functionLap(Nx,Ny,gamma,gLaplace,xgradeC,ygradeC)

      call GridifyingBetaRho(Nx,Ny,Nc,cells,grid,beta,ro,betagrid
     . ,rhogrid)


        speed=0.05/sqrt(dke0*0.024)
!        speed=0.0;

      do k= 1,Nc
        i=ceiling(cells(k,1)/dx)
        j=ceiling(cells(k,2)/dy)
        grad=sqrt(xgradeC(i,j)**2 +ygradeC(i,j)**2)
        if(rhogrid(i,j) .ge. 0.6 .and. grad .ge. 1.0
     .   .and. t/dk1 .gt. 50)then
!          write(6,*) 'Atempt move'
            angle=atan2(ygradeC(i,j),xgradeC(i,j))
            newx=cells(k,1)+h*speed*cos(angle)
            newy=cells(k,2)+h*speed*sin(angle)
            oldi=ceiling(cells(k,1)/dx)
            oldj=ceiling(cells(k,2)/dy)
            i=ceiling(newx/dx)
            j=ceiling(newy/dy)
            if (oldi .eq. i .and. oldj .eq. j)then
                cells(k,1)=newx
                cells(k,2)=newy
            else

                    if(i .gt.Nx)then
                        newx=newx-Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(i .lt.1)then
                        newx=newx+Nx*dx
                        i=ceiling(newx/dx)
                    endif

                    if(j .gt. Ny)then
                        newy=newy-Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(j .lt. 1)then
                        newy=newy+Ny*dy
                        j=ceiling(newy/dy)
                    endif

                    if(grid(i,j) .eq. 0) then
!                  write(6,*) 'Actual move'
                        cells(k,1)=newx
                        cells(k,2)=newy
                        grid(i,j)=1
                        grid(oldi,oldj)=0
                    endif
            endif


        endif
      enddo

      return
      end
