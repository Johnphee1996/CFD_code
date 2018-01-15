program   uns_fohs_diff
!
implicit none
!
integer ie,ip,itime,iel,ier,ntime,i
integer nelem,npoin,itau,ngauss
real*8,allocatable :: x(:),xc(:),x_g(:),wg(:),flux1(:),flux2(:)
real*8,allocatable :: unkno(:,:),rhsel(:,:),unkno_n(:,:)
real*8,allocatable :: unkno_e(:,:)
!
real*8 L,Lr0,Tr,c,nu,pi,dx,xg
real*8 ue,pe,b2
real*8 ul,ur,pl,pr,ua,pa,du,dp
real*8 dt,dtau,tend
real*8 lambda
real*8 res1,res2,res1_0,res2_0
real*8 l2_u,l2_p
real*8 tol_tau
real*8 t1,t2
!
write (*,*) 'Nelem = ?'
read (*,*) nelem
!
write (*,*) 'dt = ?'
read (*,*) dt
!
write (*,*) 'dtau = ?'
!
read (*,*) dtau
!
write (*,*) 'ntime = ?'
read (*,*) ntime
!
call cpu_time(t1)
!
tend = dt*ntime
!
npoin = nelem + 1
!
tol_tau = 1.d-6
!
ngauss = 3
allocate(wg(ngauss)) 
allocate(x_g(ngauss))
wg(1) = 5.d0/9.d0
wg(2) = 8.d0/9.d0
wg(3) = 5.d0/9.d0
x_g(1) = -sqrt(3.d0/5.d0)
x_g(2) = 0.d0
x_g(3) = sqrt(3.d0/5.d0)
!
allocate(x(npoin))
allocate(xc(nelem))
allocate(unkno(2,nelem))
allocate(unkno_n(2,nelem))
allocate(rhsel(2,nelem))
allocate(flux1(npoin))
allocate(flux2(npoin))
allocate(unkno_e(2,nelem))
!
pi = 4.d0*atan(1.d0)
!
!---------------
!Grid generation
!---------------
!
L  = 1.d0
dx = L/nelem
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
nu = 0.06d0
 c = 50.d0
!
Lr0 = 0.5d0*L/pi
!
Tr = Lr0*Lr0/nu
!
lambda = sqrt(nu/Tr)
!
unkno = 0.d0
!
do ie = 1, nelem
    !
    do i = 1, ngauss
    !
    xg=0.5d0*dx*x_g(i)+xc(ie)
    unkno(1,ie) = unkno(1,ie) + c*sin(pi*xg/L)*0.5d0*wg(i)
    unkno(2,ie) = unkno(2,ie) + c*pi/L*cos(pi*xg/L)*0.5d0*wg(i)
    !
    end do
    !
end do
!
!------------------
!Physical time loop
!------------------
!
do itime = 1, ntime
    !
    !initial n-level solution
    !
    write (*,*) 'itime = ',itime
    unkno_n(:,:) = unkno(:,:)
    !
    !--------------
    !dual time loop
    !--------------
    !
    do itau = 1, 1000
        !
        !get r.h.s.
        !
        rhsel(:,:) = 0.d0
        !
        do ip = 1, npoin
            !
            iel = ip - 1
            ier = ip
            !
            if (ip .eq. 1) then
                !
                ul = 0.d0
                ur = unkno(1,1) - unkno(2,1)*0.5d0*dx
                !ul = -ur
                !
                pr = unkno(2,1)
                pl = pr
                !
            elseif (ip .eq. npoin) then
                !
                ul = unkno(1,nelem) + unkno(2,nelem)*0.5d0*dx
                ur = 0.d0
                !
                pl = unkno(2,nelem)
                pr = pl
                !
            else
                !
                ul = unkno(1,iel) + unkno(2,iel)*0.5d0*dx
                ur = unkno(1,ier) - unkno(2,ier)*0.5d0*dx
                !
                pl = unkno(2,iel)
                pr = unkno(2,ier)
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
            rhsel(2,ie) = rhsel(2,ie) - 0.5d0*dx*(flux1(ie+1) + flux1(ie)) &
                          & - (flux2(ie+1) - flux2(ie))
            !
            !source term
            !
            rhsel(2,ie) = rhsel(2,ie) - unkno(2,ie)/Tr*dx   
            !
            !dCT/dx*F domain integral
            !
            rhsel(2,ie) = rhsel(2,ie) - nu*unkno(2,ie)*dx
            !
            !physical time term
            !
            rhsel(1,ie) = rhsel(1,ie) - dx/dt*(unkno(1,ie) - unkno_n(1,ie))
            rhsel(2,ie) = rhsel(2,ie) - dx**3/12.d0/dt*(unkno(2,ie) - unkno_n(2,ie))
            !
        end do
        !
        !Solve & update
        !
        unkno(1,:) = unkno(1,:) + dtau/dx*rhsel(1,:)
        unkno(2,:) = unkno(2,:) + dtau/(dx+dx**3/12.d0)*rhsel(2,:)
        !
        !check inner residual
        !
        res1 = 0.d0
        res2 = 0.d0
        !
        do ie = 1, nelem
            !
            res1 = res1 + rhsel(1,ie)*rhsel(1,ie)
            res2 = res2 + rhsel(2,ie)*rhsel(2,ie)
            !
        end do
        !
        res1 = sqrt(res1)
        res2 = sqrt(res2)
        !
        if (itau .eq. 1) then
            !
            res1_0 = res1
            res2_0 = res2
            !
        end if
        !
        res1 = res1/res1_0
        res2 = res2/res2_0
        !
        !write (*,*) itau,'Inner residual',res1*res1_0,res1,res2*res2_0,res2
        !
        if ((itau .eq. 1000) .or. (res1 .lt. tol_tau)) then
            !
            write (*,*) itau,'Inner residual',res1
            exit
            !
        end if
        !
    end do
    !
end do
!
write(*,*) 't_end = ',tend
!
tend  = tend*pi*pi*nu
!
do ie = 1, nelem
    !
    unkno_e(1,ie) = c*exp(-tend)*sin(pi*xc(ie)/L)
    unkno_e(2,ie) = c*pi/L*exp(-tend)*cos(pi*xc(ie)/L)
    !
end do
!
!--------------------
!Compute the l2 error
!--------------------
!
l2_u = 0.d0
l2_p = 0.d0
!
do ie = 1, nelem
    !
    do i = 1, ngauss
    xg = 0.5d0*dx*x_g(i)+xc(ie)
    ue = c*exp(-tend)*sin(pi*xg/L)
    pe = c*pi/L*exp(-tend)*cos(pi*xg/L)
    b2 = xg - xc(ie)
    l2_u = l2_u + 0.5d0*dx*wg(i)*(unkno(1,ie) + unkno(2,ie)*b2 - ue)**2
    l2_p = l2_p + 0.5d0*dx*wg(i)*(unkno(2,ie) - pe)**2
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
open(1,file='output2.dat')
!
do ie = 1, nelem
    !
    write (1,*) xc(ie),unkno(1,ie),unkno(2,ie),unkno_e(1,ie),unkno_e(2,ie)
    !
end do
!
close(1)
!
open(2,access='append',file='order2.dat')
    !
    write(*,*) 'L2-error'
    !
    write(*,*) 'Log(DOF^-1)     Log(e_u)          Log(e_p)'
    write(*,*) -log10(2.d0*nelem),log10(l2_u),log10(l2_p)
    write(2,*) -log10(2.d0*nelem),log10(l2_u),log10(l2_p)
    !
close(2)
!
call cpu_time(t2)
!
write(*,*) 'Total time =', t2-t1
!
end program
