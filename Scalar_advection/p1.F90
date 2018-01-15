program   Sca_adv
!
implicit none
!
integer nelem,npoin,i,j,ie,ip,t,tend
integer nTVD,ngauss
real*8 dx,CFL,dt,t1,t2,xg,t_end
real*8 dcx,ul
real*8,allocatable :: x(:),xc(:)
real*8,allocatable :: y(:,:),temp(:,:),rhsel(:,:)
real*8,allocatable :: x_g(:),wg(:)
real*8,allocatable :: alpha(:),beta(:),flux(:)
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
dcx = 1.d-3*dx
nelem = 1.d0/dx
npoin = nelem + 1
dt = dx*CFL
tend = t_end/dt + 1
!
write (*,*) 'tend = ', tend
!
allocate(x(npoin))
allocate(flux(npoin))
allocate(xc(nelem))
allocate(y(2,nelem))
allocate(temp(2,nelem))
allocate(rhsel(2,nelem))
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
nTVD = 3
allocate(alpha(nTVD))
allocate(beta(nTVD))
alpha(1) = 0.d0
alpha(2) = 3.d0/4.d0
alpha(3) = 1.d0/3.d0
beta(1) = 1.d0
beta(2) = 1.d0/4.d0
beta(3) = 2.d0/3.d0
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
    !Even value of the element
    !
    do i = 1, ngauss
        !
        xg = 0.5d0*dx*x_g(i) + xc(ie)
        !
        if (ie .le. 60) then
            !
            y(1,ie) = y(1,ie) + 0.5d0*wg(i)*exp(-200.d0*(xg-0.3d0)**2.d0)
            !
        elseif(ie .le. 81) then
            !
            y(1,ie) = y(1,ie) + 0.5d0*wg(i)*1.d0
            !
        else
            !
            y(1,ie) = y(1,ie) + 0.5d0*wg(i)*0.d0
            !
        end if
        !
    end do
    !
    !Derivative value of the element
    !
    if (ie .le. 60) then 
        !
        y(2,ie) = (exp(-200.d0*(xc(ie)+dcx-0.3d0)**2) - exp(-200.d0*(xc(ie)-dcx-0.3d0)**2))&
                  &/2.d0/dcx
        !          
    elseif (ie .le. 81) then
        !
        y(2,ie) = (1.d0 - 1.d0)/2.d0/dcx
        !
    else
        !
        y(2,ie) = (0.d0 - 0.d0)/2.d0/dcx
        !
    end if
    !
end do
!
!----------------------------------
!Time evolution : TVD Runge Kutta 3
!----------------------------------
!    
do t = 1, tend
    !
    !construct the temporary variable
    !
    temp = y
    !
    do j = 1, nTVD
        !
        !get r.h.s
        !
        rhsel = 0.d0
        !
        do ip = 1, npoin
            !
            if (ip .eq. 1) then 
                !
                ul = temp(1,nelem) + temp(2,nelem)*dx/2.d0  ! Periodic Boundary Condition
                !
            else
                !
                ul = temp(1,ip-1) + temp(2,ip-1)*dx/2.d0
                !
            end if
            !
            flux(ip) = ul
            !
        end do
        !
        do ie = 1, nelem
            !
            !flux term
            !
            rhsel(1,ie) = rhsel(1,ie) - (flux(ie+1) - flux(ie))
            rhsel(2,ie) = rhsel(2,ie) - (flux(ie+1) + flux(ie))*dx/2.d0
            !
            !dB/dx*f domain integral
            !
            rhsel(2,ie) = rhsel(2,ie) + temp(1,ie)*dx
            !
        end do
        !
        !evolution
        !
        temp(1,:) = alpha(j)*y(1,:) + beta(j)*(temp(1,:) + CFL*rhsel(1,:))
        temp(2,:) = alpha(j)*y(2,:) + beta(j)*(temp(2,:) + CFL*12.d0/dx**2*rhsel(2,:))
        !
    end do
    !
    y = temp
    !
end do
!
!------
!Output
!------
!
open(1,file='outputp1.dat')
!
do ie = 1, nelem
    !
    write (1,*) xc(ie),y(1,ie),y(2,ie)
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
       


