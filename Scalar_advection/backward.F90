program   Sca_adv
!
implicit none
!
integer nelem,npoin,i,tend,t
real*8 dx,CFL,dt,t1,t2
real*8,allocatable :: x(:),y(:),temp(:)
!
write (*,*) 'dx = ?'
read (*,*) dx
!
write (*,*) 'CFL = ?'
read (*,*) CFL
!
call cpu_time(t1)
!
nelem = 5.d0/dx
npoin = nelem + 1
dt = dx*CFL
tend = 2.d0/dt + 1
allocate(x(npoin))
allocate(y(npoin))
allocate(temp(npoin))
!
!---------------
!Grid generation
!---------------
!
do i = 1, npoin
    x(i) = 1.d-2*(i-1)
end do
!
!----------------
!Initial solution
!----------------
!
do i = 1, npoin
    !
    if (i .lt. 61) then
        !
        y(i) = exp(-200.d0*((i-1)/100.d0 - 0.3d0)**2.d0)
        !
    elseif (i .lt. 82) then
        !
        y(i) = 1.d0
        !
    else
        !
        y(i) = 0.d0
        !
    end if
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
    do i = 1, npoin
        !
        if (i .eq. 1) then
            !
            y(i) = 0.d0
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
open(1,file='outputbac.dat')
!
do i = 1, npoin
    !
    write (1,*) x(i),y(i)
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

