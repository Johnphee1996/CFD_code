   program   Sca_adv
   !
   implicit none
   !
   integer nelem,npoin,i,tend,t,ie,ip,ngauss
   real*8 dx,CFL,dt,t1,t2,xg,t_end
   real*8,allocatable :: x(:),xc(:),y(:),temp(:)
   real*8,allocatable :: x_g(:),wg(:)
   !
   write (*,*) 'CFL = ?'
   read (*,*) CFL
   !
   write (*,*) 't_end = ?'
   read (*,*) t_end
   !
   call cpu_time(t1)
   !
   dx = 1.d-2
   nelem = 1.d0/dx
   npoin = nelem + 1
   dt = dx*CFL
   tend = t_end/dt + 1
   !
   write (*,*) 'tend = ', tend
   !
   allocate(x(npoin))
   allocate(xc(nelem))
   allocate(y(nelem))
   allocate(temp(nelem))
   !
   ngauss = 3
   allocate(x_g(ngauss))
   allocate(wg(ngauss))
   wg(1) = 5.d0/9.d0
   wg(2) = 8.d0/9.d0
   wg(3) = 5.d0/9.d0
   x_g(1) = -sqrt(3.d0/5.d0)
   x_g(2) = 0.d0
   x_g(3) = sqrt(3.d0/5.d0)
   !
   !---------------
   !Grid generation
   !---------------
   !
   do ip = 1, npoin
       x(ip) = dx*(ip - 1)
   end do
   !
   do ie = 1, nelem
       xc(ie) = 0.5d0*(x(ie) + x(ie+1))
   end do
   !
   !----------------
   !Initial solution
   !----------------
   !
   y = 0.d0
   !
   do ie = 1, nelem
       !
       do i = 1, ngauss
           !
           xg = 0.5d0*dx*x_g(i) + xc(ie)
           !
           if (ie .le. 60) then
               !
               y(ie) = y(ie) + 0.5d0*wg(i)*exp(-200.d0*(xg-0.3d0)**2.d0)
               !
           elseif(ie .le. 81) then
               !
               y(ie) = y(ie) + 0.5d0*wg(i)*1.d0
               !
           else
               !
               y(ie) = y(ie) + 0.5d0*wg(i)*0.d0
               !
           end if
           !
       end do
       !
   end do 
   !
   !--------------
   !Time evolution
   !--------------
   !
   do t = 1, tend
       !
       !construct the temporary variable
       !
       temp = y
       !
       do i = 1, nelem
           !
           if (i .eq. 1) then
               !
               y(i) = temp(i) - CFL*(temp(i) - temp(i+99)) !Periodic boundary condition
               !
           else
               !
               y(i) = temp(i) - CFL*(temp(i) - temp(i-1))
               !
           end if
           !
       end do
       !
   end do
   !
   !------
   !Output
   !------
   !
   open(1,file='outputp0.dat')
   !
   do ie = 1, nelem
       !
       write (1,*) xc(ie),y(ie)
       !
   end do
   !
   close(1)
   !
   call cpu_time(t2)
   !
   write (*,*) 'Total time =', t2-t1
   !
end program
       


