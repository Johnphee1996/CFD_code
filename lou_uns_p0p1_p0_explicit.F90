program   uns_fohs_diff
!
implicit none
!
integer ie,ip,ipl,ipr,itime,iel,ier,ntime,i
integer nelem,npoin,itau,ngauss
real*8,allocatable :: x(:),xc(:),x_g(:),wg(:)
real*8,allocatable :: unkno(:,:),rhsel(:,:),unkno_n(:,:)
real*8,allocatable :: unkno_e(:,:)
!
real*8 L,Lr0,Tr,c,nu,pi,dx,xg
real*8 b2l,b2r
real*8 ul,ur,pl,pr,ua,pa,du,dp,xl,xr
real*8 dt,dtau,tend
real*8 lambda,flux1,flux2
real*8 res1,res2,res1_0,res2_0
real*8 l2_u,l2_p
real*8 tol_tau
!

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
write (*,*) 'ntime =?'
read (*,*) ntime
write(*,*) 'ntime =',ntime
!
tend = dt*ntime
write (*,*) 't_end =',tend
npoin = nelem + 1
!
!
tol_tau = 1.d-6
ngauss=3
allocate(wg(ngauss)) 
allocate(x_g(ngauss))
wg(1)=5.d0/9.d0
wg(2)=8.d0/9.d0
wg(3)=5.d0/9.d0
x_g(1)=-sqrt(3.d0/5.d0)
x_g(2)=0.d0
x_g(3)=sqrt(3.d0/5.d0)
allocate(x(npoin))
allocate(xc(nelem))
allocate(unkno(2,nelem))
allocate(unkno_n(2,nelem))
allocate(rhsel(2,nelem))
!
!
pi = 4.d0*atan(1.d0)
!
!
!Grid generation
!
L  = 1.d0
dx = L/nelem
!
!
do ip = 1, npoin
    x(ip) = dx*(ip-1)
end do
!
do ie = 1, nelem
    xc(ie) = 0.5d0*dx + dx*(ie-1)
end do
!
!
!
!Initial solution
!
nu = 0.06d0
c  = 50
!
Lr0 = 0.5d0*L/pi
!
Tr = Lr0*Lr0/nu
!
lambda = sqrt(nu/Tr)
!
write (*,*) 'dt =',dt
write (*,*) 'dtau = ',dtau
!
unkno = 0.d0
!
do ie = 1, nelem
    !
    xl = x(ie)
    xr = x(ie+1)
    !
do i=1,ngauss
    xg=0.5d0*dx*x_g(i)+xc(ie)
    unkno(1,ie) = unkno(1,ie)+c*sin(pi*xg/L)*0.5d0*wg(i)
    unkno(2,ie) = unkno(2,ie)+c*pi/L*cos(pi*xg/L)*0.5d0*wg(i)
end do
    !
    !
end do
!
!
!Physical time loop
!
do itime = 1,ntime
    !
    !initial n-level solution
    !
    write (*,*) 'itime = ',itime
    unkno_n(:,:) = unkno(:,:)
    !
    !dual time loop
    !
    do itau = 1, 1000
        !
        !get r.h.s.
        !
        rhsel(:,:) = 0.d0
        !
        do ip = 1, npoin
            !
            iel = ip -1
            ier = ip
            !
            b2l = 0.5d0*dx
            b2r =-0.5d0*dx
            !
            if (ip .eq. 1) then
                !
                ul = 0.d0
                ur = unkno(1,1)+unkno(2,1)*b2r
                !ul = -ur
                !
                pl = unkno(2,1)
                pr = unkno(2,1)
                !
            elseif (ip .eq. npoin) then
                !
                ul = unkno(1,nelem)+unkno(2,nelem)*b2l
                ur = 0.d0
                !
                pl = unkno(2,nelem)
                pr = pl
            else
                !
                ul = unkno(1,iel)+unkno(2,iel)*b2l
                ur = unkno(1,ier)+unkno(2,ier)*b2r
                !
                pl = unkno(2,iel)
                pr = unkno(2,ier)
                !
            end if
            !
            ua = 0.5d0*(ul+ur)
            pa = 0.5d0*(pl+pr)
            !
            du = ur-ul
            dp = pr-pl
            !
            flux1 = -nu*pa-0.5d0*lambda*du
            flux2 = -ua/Tr-0.5d0*lambda*dp
            !
            !
            if (iel .ge. 1) then
                rhsel(1,iel) = rhsel(1,iel) - flux1
                rhsel(2,iel) = rhsel(2,iel) - flux1*b2l - flux2
            end if
            !
            if (ier .le. nelem) then
                rhsel(1,ier) = rhsel(1,ier) + flux1
                rhsel(2,ier) = rhsel(2,ier) + flux1*b2r + flux2
            end if
            !
        end do
        !
        !source term
        !
        do ie = 1, nelem
            !
            xl = x(ie)
            xr = x(ie+1)
            !
            rhsel(2,ie) = rhsel(2,ie) - unkno(2,ie)/Tr*dx
            !
        end do
        !
        !dCT/dx*F domain integral
        !
        do  ie = 1, nelem
        !
        rhsel(2,ie) = rhsel(2,ie) -nu*unkno(2,ie)*dx
        !
        end do
        !physical time term
        !
        do ie = 1, nelem
            !
            rhsel(1,ie) = rhsel(1,ie) - dx/dt*(unkno(1,ie)-unkno_n(1,ie))
            rhsel(2,ie) = rhsel(2,ie) - dx**3/12.d0/dt*(unkno(2,ie)-unkno_n(2,ie))
            !
        end do
        !
        !
        !Solve & update
        !
        do ie = 1, nelem
            !
            unkno(1,ie) = unkno(1,ie) + dtau/dx*rhsel(1,ie)
            unkno(2,ie) = unkno(2,ie) + dtau/(dx+dx**3/12.d0))*rhsel(2,ie)
            !
        end do
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
            res1_0 = res1
            res2_0 = res2
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
write(*,*) 't_end',tend
!
tend  = tend*pi*pi*nu
!
allocate(unkno_e(2,nelem))
!
do ie = 1, nelem
    !
    unkno_e(1,ie) = c*exp(-tend)*sin(pi*xc(ie)/L)
    unkno_e(2,ie) = c*pi/L*exp(-tend)*cos(pi*xc(ie)/L)
    !
end do
!
!l2 error
!
l2_u = 0.d0
l2_p = 0.d0
!
do ie = 1, nelem
    !
do i = 1,ngauss
    xg=0.5d0*dx*x_g(i)+xc(ie)
    l2_u = l2_u + dx/2*wg(i)*(unkno(1,ie)+unkno(2,ie)*(xg-xc(ie))-unkno_e(1,ie))**2
    l2_p = l2_p + dx/2*wg(i)*(unkno(2,ie)-unkno_e(2,ie))**2
    !
end do
 !
end do
!
l2_u = sqrt(l2_u)
l2_p = sqrt(l2_p)
!
write(*,*) 'L2-error'
!
write(*,*) 'Log(h)',log10(dx)
write(*,*) 'Log(e_u)',log10(l2_u)
write(*,*) 'Log(e_p)',log10(l2_p)
!
!
!Output
!
open(1,file='output.dat')
!
do ie = 1, nelem
    !
    write (1,*) xc(ie),unkno(1,ie),unkno(2,ie),unkno_e(1,ie),unkno_e(2,ie)
    !
end do
!
close(1)
!
end
