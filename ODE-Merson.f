      subroutine ODE(t,Nx,Ny,Nc,beta,gamma,ro,cells)
      
      implicit none
      double precision t
      integer Nx, Ny,Nc,i,j

      double precision beta(Nc),gamma(Nx,Ny),ro(Nc)
      double precision beta0(Nc),gamma0(Nx,Ny),ro0(Nc)
      double precision betak1(Nc),gammak1(Nx,Ny),rok1(Nc)
      double precision betak2(Nc),gammak2(Nx,Ny),rok2(Nc)
      double precision betak3(Nc),gammak3(Nx,Ny),rok3(Nc)
      double precision betak4(Nc),gammak4(Nx,Ny),rok4(Nc)
      double precision betak5(Nc),gammak5(Nx,Ny),rok5(Nc)
      double precision b1(Nc),g1(Nx,Ny),r1(Nc)
      double precision vdx(Nx,Ny), vdy(Nx,Ny)
      double precision cells(Nc,2)
      integer grid(Nx,Ny)

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      double precision tau,h,err
      integer iteration, index

      tau=0.d0
      h=dt


!      call flow(t,Nx,Ny,vdx,vdy)
      do j=1,Ny
       do i=1,Nx
        vdx(i,j)=0.0
        vdy(i,j)=0.0
       enddo
      enddo



 13   call FromCellToGrid(Nx,Ny,Nc,cells,grid)

      do j=1,Ny
         do i=1,Nx
            gamma0(i,j)=gamma(i,j)
         enddo
      enddo

      do i=1,Nc
         ro0(i)=ro(i)
         beta0(i)=beta(i)
      enddo

      iteration=0

      call rs(t,Nx,Ny,Nc,beta0,gamma0,ro0,betak1,gammak1,rok1,cells,vdx)
!     Runge-Kutta-Merson Method

 16   do j=1,Ny
         do i=1,Nx
            gamma(i,j)=gamma0(i,j)+h*gammak1(i,j)/3
         enddo
      enddo

      do i=1,Nc
        ro(i)=ro0(i)+h*rok1(i)/3
        beta(i)=beta0(i)+h*betak1(i)/3
      enddo


      call rs(t+h/3,Nx,Ny,Nc,beta,gamma,ro,betak2,gammak2,rok2,cells
     . ,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+gammak2(i,j))/6
       enddo
      enddo

      do i=1,Nc
        ro(i)=ro0(i)+h*(rok1(i)+rok2(i))/6
        beta(i)=beta0(i)+h*(betak1(i)+betak2(i))/6
      enddo

      call rs(t+h/3,Nx,Ny,Nc,beta,gamma,ro,betak3,gammak3,rok3,cells
     . ,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+3*gammak3(i,j))/8
       enddo
      enddo

      do i=1,Nc
        ro(i)=ro0(i)+h*(rok1(i)+3*rok3(i))/8
        beta(i)=beta0(i)+h*(betak1(i)+3*betak3(i))/8
      enddo


      call rs(t+h/2,Nx,Ny,Nc,beta,gamma,ro,betak4,gammak4,rok4,cells
     . ,vdx)

       do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)-3*gammak3(i,j)
     .   +4*gammak4(i,j))/2
       enddo
      enddo

      do i=1,Nc
        ro(i)=ro0(i)+h*(rok1(i)-3*rok3(i)
     .   +4*rok4(i))/2
        beta(i)=beta0(i)+h*(betak1(i)-3*betak3(i)
     .   +4*betak4(i))/2
      enddo



      call rs(t+h,Nx,Ny,Nc,beta,gamma,ro,betak5,gammak5,rok5,cells
     . ,vdx)

      do j=1,Ny
       do i=1,Nx
        gamma(i,j)=gamma0(i,j)+h*(gammak1(i,j)+4*gammak4(i,j)
     .   +gammak5(i,j))/6
       enddo
      enddo

      do i=1,Nc
        ro(i)=ro0(i)+h*(rok1(i)+4*rok4(i)
     .   +rok5(i))/6
        beta(i)=beta0(i)+h*(betak1(i)+4*betak4(i)
     .   +betak5(i))/6
      enddo



      do j=1,Ny
       do i=1,Nx
        g1(i,j)=gamma(i,j)-h*(gammak1(i,j)+gammak2(i,j)+gammak3(i,j)
     .   +gammak4(i,j)+gammak5(i,j))/5-gamma0(i,j)
       enddo
      enddo

      do i=1,Nc
        r1(i)=ro(i)-h*(rok1(i)+rok2(i)+rok2(i)+rok3(i)
     .   +rok4(i)+rok5(i))/5-ro0(i)
        b1(i)=beta(i)-h*(betak1(i)+betak2(i)+betak3(i)
     .   +betak4(i)+betak5(i))/5-beta0(i)
      enddo

      err=0.d0
      index=0
      do j=1,Ny
         do i=1,Nx
            err=max(abs(g1(i,j)),err)
            if ( gamma(i,j) .lt. 0)then
               index=1
            endif
         enddo
      enddo

      do i=1,Nc
         err=max(abs(b1(i)),abs(r1(i)))
         if (beta(i) .lt. 0 .or. ro(i) .lt. 0)then
            index=1
         endif
      enddo



      if (err .gt. tol .or. index .eq. 1) then
      h=h/2
      iteration=iteration+1

        if (iteration .gt. 2) then
            write(6,*) 't =',t,' index =',index, 'iteration=',iteration
        endif
        if (iteration .gt. 40) then
            write(6,*) 'Emergency Exit'
            call EXIT(0)
        endif

      go to 16
      endif


      t=t+h
      tau=tau+h

      call UpdateDiscreteCells(t,h,Nx,Ny,Nc, gamma, ro, beta,cells)

      h=dt

      if (tau + h .le. tout+tol*dt) then


       go to 13
      elseif(tau .lt. tout-tol*dt)then
         h = tout - tau

         go to 13
      endif


      return
      end




