program ste_exp
!
implicit none
!
integer nelem,npoin,ie,ip,ntime
integer i,itau,ngauss,iel,ier
!
real*8,allocatable :: x(:),xc(:),x_g(:),wg(:),flux1(:),flux2(:)
real*8,allocatable :: unkno(:,:),unkno_e(:,:),rhsel(:,:)
!
real*8 L,Lr0,Tr,nu,pi,dx,xg
real*8 ue,pe,t1,t2
real*8 b2,b3
real*8 ul,ur,pl,pr,ua,pa,du,dp
real*8 dt,lambda,dcx
real*8 res1,res2,res3,res1_0,res2_0,res3_0
real*8 l2_u,l2_p
real*8 tol_tau
!
write (*,*) 'Nelem = ?'
read (*,*) Nelem
!
write (*,*) 'dt = ?'
read (*,*) dt
!
call cpu_time(t1)
!
ntime = 1000000
!
npoin = nelem + 1
tol_tau = 1.d-11
!
ngauss = 3
allocate(x_g(ngauss))
allocate(wg(ngauss))
x_g(1) = -sqrt(3.d0/5.d0)
x_g(2) = 0.d0
x_g(3) = sqrt(3.d0/5.d0)
wg(1) = 5.d0/9.d0
wg(2) = 8.d0/9.d0
wg(3) = 5.d0/9.d0
!
allocate(x(npoin))
allocate(xc(nelem))
allocate(unkno(3,nelem))
allocate(unkno_e(2,nelem))
allocate(flux1(npoin))
allocate(flux2(npoin))
allocate(rhsel(3,nelem))
!
pi = 4.d0*atan(1.d0)
!
!---------------
!Grid generation
!---------------
!
L = 1.d0
dx = L/nelem
dcx = 1.d-3*dx
!
do ip = 1, npoin
    x(ip) = dx*(ip-1)
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
nu = 1.d0
Lr0 = 0.5d0*L/pi
Tr = Lr0*Lr0/nu
lambda = sqrt(nu/Tr)
!
unkno = 0.d0
!
do ie = 1, nelem
    !
    !Even value of the element
    !
    do i = 1, ngauss
        !
        xg = 0.5d0*dx*x_g(i) + xc(ie)
        unkno(1,ie) = unkno(1,ie) + 0.5d0*wg(i)*xg*(xg-1.d0)
        unkno(2,ie) = unkno(2,ie) + 0.5d0*wg(i)*(2.d0*xg - 1.d0)
        !
    end do
    !
    !Derivative value of the element
    !
    unkno(3,ie) = unkno(3,ie) + ((2.d0*(xc(ie)+dcx) - 1.d0) - &
                  & (2.d0*(xc(ie)-dcx) - 1.d0))*0.5d0/dx
    !
end do
!
!----------------
!Pseudo time loop
!----------------
!
do itau = 1, ntime
    !
    !get r.h.s
    !
    rhsel = 0.d0
    !
    do ip = 1, npoin
        !
        iel = ip - 1
        ier = ip
        !
        if (ip .eq. 1) then
            !
            ul = 0.d0
            ur = unkno(1,1) - unkno(2,1)*0.5d0*dx + unkno(3,1)*dx**2/12.d0
            !
            pr = unkno(2,1) - unkno(3,1)*0.5d0*dx
            pl = pr
            !
        elseif (ip .eq. npoin) then
            !
            ul = unkno(1,nelem) + unkno(2,nelem)*0.5d0*dx + unkno(3,nelem)*dx**2/12.d0
            ur = 0.d0
            !
            pl = unkno(2,nelem) + unkno(3,nelem)*0.5d0*dx
            pr = pl
            !
        else
            !
            ul = unkno(1,iel) + unkno(2,iel)*0.5d0*dx + unkno(3,iel)*dx**2/12.d0
            ur = unkno(1,ier) - unkno(2,ier)*0.5d0*dx + unkno(3,ier)*dx**2/12.d0
            !
            pl = unkno(2,iel) + unkno(3,iel)*0.5d0*dx   
            pr = unkno(2,ier) - unkno(3,ier)*0.5d0*dx
            !
        end if
        !
        ua = 0.5d0*(ul + ur)
        pa = 0.5d0*(pl + pr)
        !
        du = ur - ul
        dp = pr - pl
        !
        flux1(ip) = -nu*pa - 0.5d0*lambda*du
        flux2(ip) = -ua/Tr - 0.5d0*lambda*dp
        !
    end do
    !
    do ie = 1, nelem
        !
        !flux term
        !
        rhsel(1,ie) = rhsel(1,ie) - (flux1(ie+1) - flux1(ie))
        rhsel(2,ie) = rhsel(2,ie) - 0.5d0*dx*(flux1(ie+1) + flux1(ie))&
                      & - (flux2(ie+1) - flux2(ie))
        rhsel(3,ie) = rhsel(3,ie) - dx**2/12.d0*(flux1(ie+1) - flux1(ie))&
                      & - 0.5d0*dx*(flux2(ie+1) + flux2(ie))
        !
        do i = 1, ngauss
            !
            xg = 0.5d0*dx*x_g(i) + xc(ie)
            b2 = (xg - xc(ie))
            b3 = 0.5d0*(xg - xc(ie))**2 - dx**2/24.d0
            !
            !source term
            !
            rhsel(1,ie) = rhsel(1,ie) + 0.5d0*dx*wg(i)*nu*pi*pi*sin(pi*xg)
            rhsel(2,ie) = rhsel(2,ie) + 0.5d0*dx*wg(i)*b2*nu*pi*pi*sin(pi*xg) &
                          & - 0.5d0*dx*wg(i)*(unkno(2,ie) + unkno(3,ie)*b2)/Tr
            rhsel(3,ie) = rhsel(3,ie) + 0.5d0*dx*wg(i)*b3*nu*pi*pi*sin(pi*xg) &
                          & - 0.5d0*dx*wg(i)*b2*(unkno(2,ie) + unkno(3,ie)*b2)/Tr
            !
            !dCT/dx*F domain integral
            !
            rhsel(2,ie) = rhsel(2,ie) - 0.5d0*dx*wg(i)*nu*(unkno(2,ie) + unkno(3,ie)*b2)
            rhsel(3,ie) = rhsel(3,ie) - 0.5d0*dx*wg(i)*nu*b2*(unkno(2,ie) + unkno(3,ie)*b2) &
                          & - 0.5d0*dx*wg(i)*(unkno(1,ie) + unkno(2,ie)*b2 + unkno(3,ie)*b3)/Tr
            !
        end do
        !
    end do
    !
    !Solve & update
    !
    unkno(1,:) = unkno(1,:) + dt/dx*rhsel(1,:)
    unkno(2,:) = unkno(2,:) + dt/(dx + dx**3/12.d0)*rhsel(2,:)
    unkno(3,:) = unkno(3,:) + dt/(dx**3/12.d0 + dx**5/720.d0)*rhsel(3,:)
    !
    !check inner residual
    !
    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    !
    do ie = 1, nelem
        !
        res1 = res1 + rhsel(1,ie)*rhsel(1,ie)
        res2 = res2 + rhsel(2,ie)*rhsel(2,ie)
        res3 = res3 + rhsel(3,ie)*rhsel(3,ie)
        !
    end do
    !
    res1 = sqrt(res1)
    res2 = sqrt(res2)
    res3 = sqrt(res3)
    !
    if (itau .eq. 1) then
        !
        res1_0 = res1
        res2_0 = res2
        res3_0 = res3
        !
    end if
    !
    res1 = res1/res1_0
    res2 = res2/res2_0
    res3 = res3/res3_0
    !
    !Judge the convergence condition
    !
    if ((itau .eq. ntime) .or. (res1 .lt. tol_tau)) then
        !
        write(*,*) itau, 'Inner residual',res1
        exit
        !
    end if
    !
end do
!
!--------------------
!Compute the l2_error
!--------------------
!
do ie = 1, nelem
    !
    unkno_e(1,ie) = sin(pi*xc(ie)/L)
    unkno_e(2,ie) = pi/L*cos(pi*xc(ie)/L)
    !
end do
!
l2_u = 0.d0
l2_p = 0.d0
!
do ie = 1, nelem
    !
    do i = 1, ngauss
        !
        xg = 0.5d0*dx*x_g(i) + xc(ie)
        ue = sin(pi*xg/L)
        pe = pi/L*cos(pi*xg/L)
        b2 = xg - xc(ie)
        b3 = 0.5d0*(xg - xc(ie))**2 - dx**2/24.d0
        l2_u = l2_u + 0.5d0*dx*wg(i)*(unkno(1,ie) + unkno(2,ie)*b2 + unkno(3,ie)*b3 - ue)**2
        l2_p = l2_p + 0.5d0*dx*wg(i)*(unkno(2,ie) + unkno(3,ie)*b2 - pe)**2
        !
    end do
    !
end do
!
l2_u = sqrt(l2_u)
l2_p = sqrt(l2_p)
!
!------
!Output
!------
!
open(1,file='output4.dat')
    !
    do ie = 1, nelem
        !
        write(1,*) xc(ie),unkno(1,ie),unkno(2,ie),unkno_e(1,ie),unkno_e(2,ie)
        !
    end do
    !
close(1)
!
open(2,access='append',file='order4.dat')
    !
    write(*,*) 'L2_error'
    !
    write(*,*) 'Log(DOF^-1)         Log(e_u)         Log(e_p)'
    write(*,*) -log10(3.d0*nelem),log10(l2_u),log10(l2_p)
    write(2,*) -log10(3.d0*nelem),log10(l2_u),log10(l2_p)
    !
close(2)
!
call cpu_time(t2)
!
write(*,*) 'Total time = ', t2-t1
!
end program
      
      
      