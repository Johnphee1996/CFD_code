  program   uns_fohs_diff
  !
  implicit none
  !
  integer ie,ip,ipl,ipr,itime,iel,ier,ntime
  integer nelem,npoin,itau
  real*8,allocatable :: x(:),xc(:)
  real*8,allocatable :: unkno(:,:),rhsel(:,:),unkno_n(:,:)
  real*8,allocatable :: unkno_e(:,:)
  !
  real*8 L,Lr0,Tr,c,nu,pi,dx
  real*8 ul,ur,pl,pr,ua,pa,du,dp,xl,xr
  real*8 dt,dtau,tend
  real*8 lambda,flux1,flux2
  real*8 res1,res2,res1_0,res2_0
  real*8 l2_u,l2_p
  real*8 tol_tau
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
  !
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
  do ie = 1, nelem
      !
      !xl = x(ie)
      !xr = x(ie+1)
      !ipl = ie
      !ipr = ipl + 1
      !unkno(1,ie) = c*L/pi*(cos(pi*xl/L)-cos(pi*xr/L))/dx
      unkno(1,ie) = c*sin(pi*xc(ie)/L)
      unkno(2,ie) = c*pi/L*cos(pi*xc(ie)/L)
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
              if (ip .eq. 1) then
                  !
                  ul = 0.d0
                  ur = unkno(1,1)
                  !ul = -ur
                  !
                  pl = unkno(2,1)
                  pr = unkno(2,1)
                  !
              elseif (ip .eq. npoin) then
                  !
                  ul = unkno(1,nelem)
                  ur = 0.d0
                  !
                  pl = unkno(2,nelem)
                  pr = pl
              else
                  !
                  ul = unkno(1,iel)
                  ur = unkno(1,ier)
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
              if (iel .ge. 1) then
                  rhsel(1,iel) = rhsel(1,iel) - flux1
                  rhsel(2,iel) = rhsel(2,iel) - flux2
              end if
              !
              if (ier .le. nelem) then
                  rhsel(1,ier) = rhsel(1,ier) + flux1
                  rhsel(2,ier) = rhsel(2,ier) + flux2
              end if
              !
          end do
          !
          !source term
          !
          do ie = 1, nelem
              !
              !xl = x(ie)
              !xr = x(ie+1)
              !
              !rhsel(1,ie) = rhsel(1,ie) + nu*pi*(cos(pi*xl)-cos(pi*xr))
              rhsel(2,ie) = rhsel(2,ie) - unkno(2,ie)/Tr*dx
              !
          end do
          !
          !physical time term
          !
          do ie = 1, nelem
              !
              rhsel(1,ie) = rhsel(1,ie) - dx/dt*(unkno(1,ie)-unkno_n(1,ie))
              !
          end do
          !
          !
          !Solve & update
          !
          do ie = 1, nelem
              !
              unkno(1:2,ie) = unkno(1:2,ie) + dtau/dx*rhsel(1:2,ie)
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
  ! 
  tend  = tend*pi*pi*nu
  !
  allocate(unkno_e(2,nelem))
  !
  do ie = 1, nelem
      !
      unkno_e(1,ie) = c*exp(-tend)*sin(pi*xc(ie))
      unkno_e(2,ie) = c*pi*exp(-tend)*cos(pi*xc(ie))
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
      l2_u = l2_u + (unkno(1,ie)-unkno_e(1,ie))**2
      l2_p = l2_p + (unkno(2,ie)-unkno_e(2,ie))**2
      !
  end do
  !
  l2_u = sqrt(l2_u/nelem)
  l2_p = sqrt(l2_p/nelem)
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
