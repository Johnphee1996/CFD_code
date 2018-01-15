  program   uns_fohs_diff
  !
  implicit none
  !
  integer ie,ip,ipl,ipr,itime,iel,ier,ntime,i
  integer nelem,npoin,itau,ngauss
  real*8,allocatable :: x(:),xc(:),x_g(:),wg(:),flux1(:),flux2(:)
  real*8,allocatable :: unkno(:,:),rhsel(:,:),unkno_n(:,:)
  real*8,allocatable :: unkno_e(:,:)
  real*8,allocatable :: delun(:,:)
  real*8,allocatable :: diago(:,:,:),lower(:,:,:),upper(:,:,:),diagp(:,:,:)
  !
  real*8 L,Lr0,Tr,c,nu,pi,dx,xg,det
  real*8 b2l,b2r,ue,pe
  real*8 ul,ur,pl,pr,ua,pa,du,dp,xl,xr
  real*8 dt,dtau,tend
  real*8 lambda
  real*8 res1,res2,res1_0,res2_0
  real*8 l2_u,l2_p
  real*8 tol_tau
  !
  
  !
  write (*,*) 'Nelem = ?'
  read (*,*) nelem
  !
  !write (*,*) 'dt = ?'
  !read (*,*) dt
  !
  write (*,*) 'dtau = ?'
  !
  read (*,*) dtau
  !
  !write (*,*) 'ntime =?'
  !read (*,*) ntime
  !write (*,*) 'ntime =',ntime
  !
  ntime = 1000000
  !tend = dtau*ntime
  !write (*,*) 't_end =',tend
  npoin = nelem + 1
  !
  !
  tol_tau = 1.d-11
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
  allocate(diago(2,2,nelem))
  allocate(upper(2,2,nelem-1))
  allocate(lower(2,2,2:nelem))
  allocate(diagp(2,2,nelem))
  allocate(delun(2,nelem))
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
      xc(ie) = 0.5d0*(x(ie)+x(ie+1))
  end do
  !
  !Initial solution
  !
  nu = 1.d0
  !
  Lr0 = 0.5d0*L/pi
  !
  Tr = Lr0*Lr0/nu
  !
  lambda = sqrt(nu/Tr)
  !
  !write (*,*) 'dt =',dt
  !
  unkno = 0.d0
  !
  do ie = 1, nelem
      !
      do i = 1, ngauss
          !
          xg = 0.5d0*dx*x_g(i) + xc(ie)
          unkno(1,ie) = unkno(1,ie) + 0.5d0*wg(i)*xg*(xg-1.d0)
          unkno(2,ie) = unkno(2,ie) + 0.5d0*wg(i)*(2.d0*xg - 1.d0)
          !
      end do
      !
      !
  end do
  !
  !
  !get l.h.s.
  !
  !
  !M/dtau
  !
  diago(1,1,:) = dx/dtau
  diago(2,2,:) = dx/dtau
  !
  !-dR/du
  !
  ! Source term
  !
  do ie = 1, nelem
      !
      diago(2,2,ie) = diago(2,2,ie) + dx/Tr
      !
  end do
  !
  !Boundary flux term
  !
  do ie = 1, nelem
      !
      !
      diago(1,1,ie) = diago(1,1,ie) + lambda
      diago(2,2,ie) = diago(2,2,ie) + lambda
      diago(1,2,ie) = 0
      diago(2,1,ie) = 0
      !
      !
      if (ie .gt. 1) then
          !
          !
          lower(1,1,ie) = -0.5d0*lambda 
          lower(1,2,ie) = 0.5d0*nu
          lower(2,1,ie) = 0.5d0/Tr
          lower(2,2,ie) = -0.5d0*lambda
          !
          !
      end if
      !
      !
      if (ie .lt. nelem) then
          !
          !
          upper(1,1,ie) = -0.5d0*lambda
          upper(1,2,ie) = -0.5d0*nu
          upper(2,1,ie) = -0.5d0/Tr
          upper(2,2,ie) = -0.5d0*lambda
          !
          !
      end if
      !
      !
  end do
  !
  !Boundary term
  !
  diago(1,2,1) = diago(1,2,1) + 0.5d0*nu
  diago(2,2,1) = diago(2,2,1) - 0.5d0*lambda
  diago(1,2,nelem) = diago(1,2,nelem) - 0.5d0*nu
  diago(2,2,nelem) = diago(2,2,nelem) - 0.5d0*lambda
  !
  !
  do ie = 1, nelem
      !
      !
      det = diago(1,1,ie)*diago(2,2,ie) - diago(1,2,ie)*diago(2,1,ie)
      det = 1.d0/det
      diagp(1,1,ie) = det*diago(2,2,ie)
      diagp(1,2,ie) = -det*diago(1,2,ie)
      diagp(2,1,ie) = -det*diago(2,1,ie)
      diagp(2,2,ie) = det*diago(1,1,ie)
      !
      !
  end do
  !
  ! 
  !
  !time loop
  !
  !do itime = 1, ntime
      !
      !initial n-level solution
      !
      !write (*,*) 'itime = ',itime
      !
      do itau = 1, ntime
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
                  pr = pl
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
              flux1(ip) = -nu*pa-0.5d0*lambda*du
              flux2(ip) = -ua/Tr-0.5d0*lambda*dp
              !
          end do
          !
          do ie = 1, nelem
          !flux term
              !
              rhsel(1,ie) = rhsel(1,ie) - (flux1(ie+1) - flux1(ie))
              rhsel(2,ie) = rhsel(2,ie) - (flux2(ie+1) - flux2(ie))
              !
          !source term
              !
              do i = 1, ngauss
                !
                xg = 0.5d0*dx*x_g(i) + xc(ie)
                rhsel(1,ie) = rhsel(1,ie) + 0.5d0*dx*wg(i)*nu*pi*pi*sin(pi*xg)
                !
              end do
              !
              rhsel(2,ie) = rhsel(2,ie) - unkno(2,ie)/Tr*dx 
              !
          !dCT/dx*F domain integral
              !
              !rhsel(2,ie) = rhsel(2,ie) -nu*unkno(2,ie)*dx
              !
          !physical time term
              !
              !rhsel(1,ie) = rhsel(1,ie) - dx/dt*(unkno(1,ie)-unkno_n(1,ie))
              !rhsel(2,ie) = rhsel(2,ie) - dx**3/12.d0/dt*(unkno(2,ie)-unkno_n(2,ie))
              !
          end do
          !
          !
          !Solve 
          !
          ! call jacobi(2,nelem,diagp,rhsel,delun)
           call lusgs(2,nelem,diago,diagp,lower,upper,rhsel,delun)
          !
          ! update
          !
          unkno(:,:) = unkno(:,:) + delun(:,:)
          !
          !unkno(:,:) = unkno(:,:) + dt/dx*rhsel(:,:)
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
          if ((itau .eq. ntime) .or. (res1 .lt. tol_tau)) then
              !
              write (*,*) itau,'Inner residual',res1
              exit
              !
          end if
          !
      end do
      !
  !end do
  !
  !write(*,*) 't_end', tend
  !
  !tend  = tend*pi*pi*nu
  !
  !
  !
  do ie = 1, nelem
      !
      unkno_e(1,ie) = sin(pi*xc(ie)/L)
      unkno_e(2,ie) = pi/L*cos(pi*xc(ie)/L)
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
      !    
      do i = 1, ngauss
          !
          !
          xg = 0.5d0*dx*x_g(i)+xc(ie)
          ue = sin(pi*xg/L)
          pe = pi/L*cos(pi*xg/L)
          l2_u = l2_u + dx/2.d0*wg(i)*(unkno(1,ie)-ue)**2
          l2_p = l2_p + dx/2.d0*wg(i)*(unkno(2,ie)-pe)**2
          !
          !
      end do
      !
      !
  end do
  !
  l2_u = sqrt(l2_u)
  l2_p = sqrt(l2_p)
  !
  open(2,access='append',file='orderste1.dat')
  write(*,*) 'L2-error'
      ! 
      !
      write(*,*) 'Log(DOF^-1)     Log(e_u)          Log(e_p)'
      write(*,*) -log10(2.d0*nelem),log10(l2_u),log10(l2_p)
      write(2,*) -log10(2.d0*nelem),log10(l2_u),log10(l2_p)
      !
      !
  close(2)
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

 
   subroutine jacobi(ndofe,nelem,diagp,rhsel,delun)
          !
          implicit none
          !
          integer ndofe,nelem
  	  real*8 diagp(ndofe,ndofe,nelem)
  	  real*8 rhsel(ndofe,nelem)
  	  real*8 delun(ndofe,nelem)
          integer ie,idofe,jdofe
          !
          delun = 0.d0
          !
          do ie = 1, nelem
              !
              do idofe = 1, ndofe
                  ! 
                  do jdofe = 1, ndofe
                      !
                      delun(idofe,ie) = delun(idofe,ie)+diagp(idofe,jdofe,ie)*rhsel(jdofe,ie)
                      ! 
                  end do
                  !
              end do
              !
          end do
          !
    end subroutine jacobi
        !
        subroutine lusgs(ndofe, nelem, diago, diagp, lower, upper, rhold, delun)
           !
           implicit none
           !
           integer ndofe, nelem
           real*8 diago(ndofe,ndofe,nelem),diagp(ndofe,ndofe,nelem)
           real*8 lower(ndofe,ndofe,2:nelem),upper(ndofe,ndofe,nelem-1)
           real*8 rhold(ndofe,nelem),delun(ndofe,nelem) 
           integer :: ie,ifa,iel,ier,i,j,ii,jj
           real*8, parameter :: c00 = 0.d0
           !
           !-------------------------
           !... lower (forward) sweep
           !-------------------------
           !
           ie = 1
           delun(1:ndofe,ie) = 0.0d0
           !
           do ii = 1, ndofe
               !
               do jj = 1, ndofe
                   !
                   delun(ii,ie) = delun(ii,ie) + diagp(ii,jj,ie)*rhold(jj,ie)
                   !
               end do
               !
           end do
           !
           do ifa = 2,  nelem
               !
               iel = ifa - 1
               ier = ifa
               !
               ! compute the new delun = lower * old delun
               !
               do ii = 1, ndofe
                   !
                   do jj = 1,ndofe
                       !
                       rhold(ii,ier) = rhold(ii,ier) - lower(ii,jj,ifa)*delun(jj,iel)
                       !
                   end do
                   !
               end do
               !
               delun(1:ndofe,ier) = 0.d0
               !
               do ii = 1, ndofe
                   !
                   do jj = 1, ndofe
                       !
                       delun(ii,ier) = delun(ii,ier) + diagp(ii,jj,ier)*rhold(jj,ier)
                       !
                   end do
                   !
               end do
               !
           end do
           !
           !
           !-----------------------------------------
           !
           !...  multiply the star by diagonal matrix                                              
           !
           !-----------------------------------------
           !
           !
           do ie = 1, nelem
               !
               rhold(1:ndofe,ie) = c00
               !
               do i = 1, ndofe
                   !
                   do j = 1, ndofe
                   !
                   rhold(i,ie) = rhold(i,ie) + diago(i,j,ie)*delun(j,ie)
                   !
                   end do
                   !
               end do
               !
           end do 
           !
           !-----------------------
           ! upper (backward) sweep
           !-----------------------
           !
           ie = nelem
           delun(1:ndofe,ie) = 0.0d0
           !
           do ii = 1, ndofe
               !
               do jj = 1, ndofe
                   !
                   delun(ii,ie) = delun(ii,ie) + diagp(ii,jj,ie)*rhold(jj,ie)
                   !
               end do
               !
           end do
           !
           do ifa = nelem-1, 1, -1
               !
               iel = ifa
               ier = ifa + 1
               !
               do ii = 1, ndofe
                   !
                   do jj = 1, ndofe
                       !
                       rhold(ii,iel) = rhold(ii,iel) - upper(ii,jj,ifa)*delun(jj,ier)
                       !
                   end do
                   !
               end do
               !
               delun(1:ndofe,iel) = 0.d0
               !
               do ii = 1, ndofe
                   !
                   do jj = 1, ndofe
                       !
                       delun(ii,iel) = delun(ii,iel) + diagp(ii,jj,iel)*rhold(jj,iel)
                       !
                   end do 
                   !    
               end do
               !
           end do
           !
           !
           return
        end subroutine lusgs
        !
        !
        subroutine sgs(nter, ndofe, nelem, nafac, diago, upper, lower, diagp, rhsel, delun, intfac)
        !
        implicit none
        !
        integer nter,ndofe,nelem,nafac
        real*8 diago(ndofe,ndofe,nelem),diagp(ndofe,ndofe,nelem)
        real*8 upper(ndofe,ndofe,nafac),lower(ndofe,ndofe,nafac)
        !
        real*8 rhsel(ndofe,nelem),delun(ndofe,nelem)
        integer intfac(2,nafac)
        !
        integer iter,iel,ier,ifa,ie,ii,jj
        !
        real*8 rhold(ndofe,nelem)
        real*8 ax(ndofe,nelem)
        real*8 res,resnom
        integer ip,jp
        !
        do iter = 1, nter
            !
            rhold = rhsel
            !
            !
            ! compute the new RHS : (L+D)*du = R-U*du
            !
            do ifa = 2, nafac - 1
                !
                iel = intfac(1,ifa)
                ier = intfac(2,ifa)
                !
                do ii =1, ndofe
                    do jj = 1, ndofe
                        !
                        rhold(ii,iel) = rhold(ii,iel) - upper(ii,jj,ifa)*delun(jj,ier)
                        !
                    end do
                end do
                !
            end do
            !
            !
          !
          !----------------------
          ! lower (forward) sweep
          !----------------------
          !
          ie = 1
          delun(1:ndofe,ie) = 0.0d0
          do ii = 1, ndofe
          do jj = 1, ndofe
            delun(ii,ie) = delun(ii,ie) + diagp(ii,jj,ie)*rhold(jj,ie)
          enddo
          enddo
          !
          do ifa = 2, nafac - 1
              !
              iel = intfac(1,ifa)
              ier = intfac(2,ifa)
              !
              ! compute the new delun = lower * old delun
              !
              do ii = 1, ndofe
                  do jj = 1,ndofe
                      !
                      rhold(ii,ier) = rhold(ii,ier) - lower(ii,jj,ifa)*delun(jj,iel)
                      !
                  end do
              end do
              !
              delun(1:ndofe,ier) = 0.d0
              !
              do ii = 1, ndofe
              do jj = 1, ndofe
                delun(ii,ier) = delun(ii,ier) + diagp(ii,jj,ier)*rhold(jj,ier)
              enddo
              enddo
                      
              !
          end do
          !
          ! renew the RHS
          !
           rhold = rhsel
           !
          !
          !  compute the new RHS = R - L*du
          !
          do ifa = 2, nafac - 1
              !
              iel = intfac(1,ifa)
              ier = intfac(2,ifa)
              !
              do ii =1, ndofe
                  do jj = 1, ndofe
                      !
                      rhold(ii,ier) = rhold(ii,ier) - lower(ii,jj,ifa)*delun(jj,iel)
                      !
                  end do
              end do
              !
          end do
          !
          !-----------------------
          ! upper (backward) sweep
          !-----------------------
          !
          ie = nelem
          delun(1:ndofe,ie) = 0.0d0
          do ii = 1, ndofe
          do jj = 1, ndofe
            delun(ii,ie) = delun(ii,ie) + diagp(ii,jj,ie)*rhold(jj,ie)
          enddo
          enddo
          !
          do ifa = nafac-1, 2, -1
              !
              iel = intfac(1,ifa)
              ier = intfac(2,ifa)
              !
              do ii = 1, ndofe
              do jj = 1, ndofe
                !
                rhold(ii,iel) = rhold(ii,iel) - upper(ii,jj,ifa)*delun(jj,ier)
                !
              enddo
              enddo
              !
              delun(1:ndofe,iel) = 0.d0
              !
              do ii = 1, ndofe
              do jj = 1, ndofe
                delun(ii,iel) = delun(ii,iel) + diagp(ii,jj,iel)*rhold(jj,iel)
              enddo     
              enddo
              !
          end do
          !
          !
        end do
        !
        end subroutine sgs
        !
        subroutine monitorsgs(ndofe,nelem,nafac,diago,upper,lower,diagp, &
                              rhsel,delun,intfac,res)
        implicit none
        !
        integer nter,ndofe,nelem,nafac
        real*8 diago(ndofe,ndofe,nelem),diagp(ndofe,ndofe,nelem)
        real*8 upper(ndofe,ndofe,nafac),lower(ndofe,ndofe,nafac)
        !
        real*8 rhsel(ndofe,nelem),delun(ndofe,nelem)
        integer intfac(2,nafac)
        !
        integer iter,iel,ier,ifa,ie,ii,jj
        !
        real*8 rhold(ndofe,nelem)
        real*8 ax(ndofe,nelem)
        real*8 res,resnom
        integer ip,jp
        !
        !
        !check res
        !
        ax = 0.d0
        !
        do ie = 1, nelem
            !
            do ip = 1, ndofe
                !
                do jp = 1, ndofe
                    !
                    ax(ip,ie) = ax(ip,ie) + diago(ip,jp,ie)*delun(jp,ie)
                end do
                !
            end do
            !
        end do
        !
        !
        do ifa = 2, nafac - 1
            !
            iel = intfac(1,ifa)
            ier = intfac(2,ifa)
            !
            do ip = 1, ndofe
                !
                do jp = 1, ndofe
                    !
                    ax(ip,iel) = ax(ip,iel) + upper(ip,jp,ifa)*delun(jp,ier)
                    ax(ip,ier) = ax(ip,ier) + lower(ip,jp,ifa)*delun(jp,iel)
                    !
                end do
                !
            end do
            !
        end do
        !
        res = 0.d0
        resnom = 0.d0
        !
        do ie = 1, nelem
            !
            do ip = 1, ndofe
                !
                res = res + (rhsel(ip,ie)-ax(ip,ie))*(rhsel(ip,ie)-ax(ip,ie))
                resnom = resnom + rhsel(ip,ie)*rhsel(ip,ie)
                !
            end do
            !
        end do
        !
        res = sqrt(res)
        resnom = sqrt(resnom)
        !
        res = res/resnom
        !
        !write (*,*) iter, res
        !if (res .le. 1E-10) flg=1
        !
        end subroutine monitorsgs
