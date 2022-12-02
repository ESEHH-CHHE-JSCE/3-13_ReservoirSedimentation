!*******************************************************************************
!
!     cip-sediment-1d-06.f90
!
!     cip-cls, staggerd, regular interval grid, equal width
!
!*******************************************************************************
    module common_variables
        integer, parameter :: k = kind(1.0d0)
!
        integer, parameter :: im = 500, nstep = 518400
        integer, parameter :: nx = im-1
!
        real(k), parameter :: dt = 0.500d0, dx = 5.0d0, width = 10.0d0
        real(k), parameter :: n_man=0.05d0, hth=0.001d0
        real(k), parameter :: ds = 0.001d0, lambda=0.4d0
!
        real(k), parameter :: pi = 3.141592654d0
        real(k), parameter :: gra = 9.80665d0
        real(k), parameter :: nuw = 1.d-6
        real(k), parameter :: dwat   = 1000.d0, dsand  = 2650.d0
        real(k), parameter :: sss = (dsand-dwat)/dwat
!
        real(k), parameter :: itr_sor = 100000, sigma_sor=0.0d0, err_sor = 1.0d-37
        real(k), parameter :: alpha1 = 1.0d0
!
        real(k), parameter :: tan_repose = -0.754d0
        real(k), parameter :: taumin = 1.0d-05, kappa = 0.4d0
!
    end module common_variables
!
!*******************************************************************************
!
    program cip_sediment_1d
    use common_variables
    implicit none
    integer istep
!
    real(k), dimension (0:im) :: yh, yhn, dh
    real(k), dimension (0:im) :: ym, ymn, yms1, yms2, yms3, yms4
    real(k), dimension (0:im) :: ysumm
    real(k), dimension (0:im) :: yu, yui
    real(k), dimension (0:im) :: yzb, yzbn, yzbs1
    real(k), dimension (0:im) :: taustb, tauste, taustbc, taustec
    real(k), dimension (0:im) :: qb, qs, cb, wscbot
    real(k), dimension (0:im) :: yc, ycn, ycs1, ycs2
    real(k), dimension (0:im) :: ych, ychn, ychs1, ychs2
    real(k), dimension (0:im) :: ysumch, ysumchn, ysumchs1, ysumchs2
    real(k) taustc
!
    integer i
!
!*** set initial conditions from here*******************************************
    print *,'k=',k

    call initl(ym, yh, yzb, yc, ych)
    call iwagakieq(taustc)
!
    open(10,file='test-time_1d_06_s.dat',status='unknown')
    open(50,file='test-time_1d_06_s_in.dat',status='unknown')
    open(60,file='test-time_1d_06_s_out.dat',status='unknown')
!
    call output(ym, yh, yzb, yc, 0)
!
!*** set initial conditions to here*********************************************
!
    do istep=1,nstep
        if(mod(istep,1000).eq.0) print *,'istep=',istep
!
        call bound(ym, yh, yzb, yc, istep)
!
!*** cal sediment from here*****************************************************
!
        call calvel(yu, ym, yh, 2)
        if (istep.eq.1) call cipsum(ysumch,ych)
!
        call shift(ych, ychs1)
        call shift(ysumch, ysumchs1)
        call shift(yzb, yzbs1)
!
!*** version : cip->non-advection from here*************************************
!
        call cipcls2_1d_r(ychs1,ysumchs1,yu,istep)
        call shift(ychs1, ychs2)
        call shift(ysumchs1, ysumchs2)
!
        call calcon(ycs2, ychs2, yh)
!
!*** version : cip->non-advection to here***************************************
!
        call bedstress(taustb, tauste, taustbc, taustec, ym, yu, yh)
!        call bedload_am(qb, taustb, tauste, taustc, yzb)
        call bedload_kp(qb, taustb, tauste, taustc, yzb)
        call suspendedload(qs, taustbc, taustec, taustc)
        call cal_wscbot(ycs2,yh,taustbc,wscbot)
        call update_ch(ychn, ysumchn, ychs2, ysumchs2, yh, qs, wscbot)
        call calbed(yzb, yzbs1, yzbn, qb, qs, wscbot)
!
!*** cal sediment to here*******************************************************
!
!*** cal flow from here*********************************************************
!
!*** version : cip->non-advection from here*************************************
!
        call calvel(yui, ym, yh, 1)
        call cipsum(ysumm,ym)
!
        call shift(ym, yms1)
!
        call cipcls2_1d_r(yms1,ysumm,yui,istep)
!
        call bound(yms1, yh, yzbn, yc, istep)
        call shift(yms1, yms2)
!
!*** version : cip->non-advection to here***************************************
!
        call visclc(yms2, yms2, yms3, yh)
        call bound(yms3, yh, yzbn, yc, istep)
!
        call drhs00(yms3, yms3, yms4, yh, yzbn)
        call bound(yms4, yh, yzbn, yc, istep)
!
        call dpreseq(yms4, yh, dh)
!
        call drhs01(yms2, yms4, ymn, yh, dh, yhn, yzbn)
        call bound(ymn, yhn, yzbn, yc, istep)
!
!*** modify sediment form here**************************************************
!
        call calcon(ycn, ychn, yhn)
!
!*** modify sediment to here****************************************************
!
        call shift(yhn, yh)
        call shift(ymn, ym)
        call shift(yzbn, yzb)
        call shift(ycn, yc)
        call shift(ychn, ych)
        call shift(ysumchn, ysumch)
!        calvel(yu, ym, yh, 1)
!
!*** cal flow to here***********************************************************
!
!        if(mod(istep,21600).eq.0) then
        if(mod(istep,7200).eq.0) then
            call output(ym, yh, yzb, yc, istep)
        end if
!
        if(mod(istep,120).eq.0) then
            do i=0,0
                write(50,220) istep*dt, i, qb(i), qs(i), wscbot(i), yc(i), taustbc(i), taustec(i), taustc, yzb(i+1), yh(i), ych(i), yu(i)
            end do
            do i=nx,nx
                write(60,220) istep*dt, i, qb(i), qs(i), wscbot(i), yc(i), taustbc(i), taustec(i), taustc, yzb(i), yh(i), ych(i), yu(i)
            end do
        end if
220 format(f13.5, i8, 11f13.5)
!
    end do
!
    close(10, status='keep')
    close(50, status='keep')
    close(60, status='keep')
!
    stop
    end program cip_sediment_1d
!
!*******************************************************************************
    subroutine initl(ym, yh, yzb, yc, ych)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(out) :: ym, yh, yzb, yc, ych
!
    integer i
!
    do i=0,im
        yzb(i)= 0.00d0+(real(im-1)-real(i))*dx*1.0d0/25.0d0
!
        yh(i)=50.0d0-yzb(i)
!
!        if(yh(i).le.0.43538d0) yh(i)=0.43538d0   ! for q=10.0m3/s w=10m i=1/25 n=0.05
        if(yh(i).le.1.407607d0) yh(i)=1.407607d0   ! for q=50.0m3/s w=10m i=1/25 n=0.05
!
        ym(i)= 10.0d0/width
!
        yc(i)= 0.0d0
        ych(i)= yc(i)*yh(i)
    end do
!
    return
    end subroutine initl
!
!*******************************************************************************
    subroutine bound(ym, yh, yzb, yc, istep)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(inout)  :: ym, yh, yzb, yc
    integer,      intent(in) :: istep
!
    integer i, n
    real(k) quan_in, h_out
!
    if (istep.lt.int(24.0d0*3600.d0/dt)) then
        n=0
    else if (istep.lt.int(48.0d0*3600.d0/dt)) then
        n=1
    else if (istep.lt.int(72.0d0*3600.d0/dt)) then
        n=2
    else if (istep.lt.int(96.0d0*3600.d0/dt)) then
        n=3
    else if (istep.lt.int(120.0d0*3600.d0/dt)) then
        n=4
    else if (istep.lt.int(144.0d0*3600.d0/dt)) then
        n=5
    else if (istep.lt.int(168.0d0*3600.d0/dt)) then
        n=6
    else if (istep.lt.int(192.0d0*3600.d0/dt)) then
        n=7
    else if (istep.lt.int(216.0d0*3600.d0/dt)) then
        n=8
    else
        n=9
    end if
!
    if (istep.le.int((6.0d0*3600.0d0+24.0d0*3600.0d0*n)/dt)) then
        quan_in=10.0d0+(100.0d0-10.0d0)/(6.0d0*3600.0d0)*(dt*istep-(24.0d0*3600.0d0)*n)
    else if (istep.le.int((12.0d0*3600.0d0+24.0d0*3600.0d0*n)/dt)) then
        quan_in=100.0d0-(100.0d0-60.0d0)/(6.0d0*3600.0d0)*(dt*istep-(6.0d0*3600.0d0)-(24.0d0*3600.0d0)*n)
    else if (istep.le.int((18.0d0*3600.0d0+24.0d0*3600.0d0*n)/dt)) then
        quan_in=60.0d0-(60.0d0-30.0d0)/(6.0d0*3600.0d0)*(dt*istep-(12.0d0*3600.0d0)-(24.0d0*3600.0d0)*n)
    else if (istep.le.int((24.0d0*3600.0d0+24.0d0*3600.0d0*n)/dt)) then
        quan_in=30.0d0-(30.0d0-10.0d0)/(6.0d0*3600.0d0)*(dt*istep-(18.0d0*3600.0d0)-(24.0d0*3600.0d0)*n)
    end if
!
    quan_in=50.0d0
!
    if(mod(istep,7200).eq.0) write(*,*) 'quan_in:', istep, n, quan_in

    h_out=50.0d0
!
    yh(0)=yh(2)
    yh(1)=yh(2)
!
    yh(nx)=h_out-yzb(nx)
    yh(im)=yh(nx)
!
    ym(1)=quan_in/width
    ym(0)=quan_in/width
!
    ym(nx)=ym(nx-1)
    ym(im)=ym(nx-1)
!
    return
    end subroutine bound
!
!*******************************************************************************
    subroutine output(ym, yh, yzb, yc, istep)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yh, yzb, yc
    integer,      intent(in) :: istep
!
    integer i
    real(k), dimension (0:im) :: yu
!
    write(10,*) 'time=', istep*dt

    call calvel(yu, ym, yh, 2)
!
    do i=1, nx
        write(10,210) real(i-1)*dx, yu(i), yh(i), (yh(i)+yzb(i)), yzb(i), yc(i)
    end do
210 format(8f13.5)
!
    return
    end subroutine output
!
!*******************************************************************************
    subroutine shift(yn,y)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: yn
    real(k), dimension (0:im), intent(out) :: y
!
    integer i
!
    do i=0,im
        y(i)=yn(i)
    end do
!
    return
    end subroutine shift
!
!*******************************************************************************
    subroutine calvel(yu, ym, yh, itp)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yh
    real(k), dimension (0:im), intent(out) :: yu
    integer, intent(in)  :: itp
!
    integer i
    real(k) yhi
!
    if (itp.eq.1) then
        do i=0,nx
            yhi=(yh(i)+yh(i+1))*0.5d0
            if (yhi.ge.hth) then
                yu(i)=ym(i)/yhi
            else
                yu(i)=0.0d0
            end if
        end do
            yu(im)=yu(nx)
    else if (itp.eq.2) then
        if (yh(0).ge.hth) then
            yu(0)= ym(0)/yh(0)
        else
            yu(0)= 0.0d0
        end if
        do i=1,im
            if (yh(i).ge.hth) then
                yu(i)= (ym(i-1)+ym(i))*0.5d0/yh(i)
            else
                yu(i)= 0.0d0
            end if
        end do
    end if
!
    return
    end subroutine calvel
!
!*******************************************************************************
    subroutine calcon(yc, ych, yh)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ych, yh
    real(k), dimension (0:im), intent(out) :: yc
!
    integer i
!
    do i=0,im
        if (yh(i).ge.hth) then
            yc(i) = ych(i)/yh(i)
        else
            yc(i)= 0.0d0
        end if
    end do
!
    return
    end subroutine calcon
!
!*******************************************************************************
    subroutine update_ch(ychn, ysumchn, ych, ysumch, yh, qs, wscbot)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ych, ysumch, yh, qs
    real(k), dimension (0:im), intent(inout) :: wscbot
    real(k), dimension (0:im), intent(out) :: ychn, ysumchn
!
    integer i
!
    do i=0,nx
        if (yh(i).ge.hth) then
            ychn(i) = ych(i)+(qs(i)-wscbot(i))*dt
        else
            ychn(i)= 0.0d0
        end if
    end do
!
    do i=0,nx
        if (ychn(i).lt.0.0d0) then
            wscbot(i)=ych(i)/dt+qs(i)
            ychn(i) = 0.0d0
!            write(*,*) 'ychn<0:', i, ych(i), qs(i), wscbot(i)
        end if
    end do
    ychn(im)=ychn(nx)
!
    do i=0,nx
        ysumchn(i) = ysumch(i) + ((qs(i)-wscbot(i)) +(qs(i+1)-wscbot(i+1)))*0.5d0*dt*dx
    end do
        ysumchn(im) = ysumchn(nx)
!
    return
    end subroutine update_ch
!
!*******************************************************************************
    subroutine visclc(ym, yms, ymn, yh)   ! ym:input, yms:before, ymn:after
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yms, yh
    real(k), dimension (0:im), intent(out) :: ymn
!
    integer i
    real(k) he, hn
!
    do i=1,nx
        he=(yh(i)+yh(i+1))/2.0d0
!
        if (he.ge.hth) then
            ymn(i)=yms(i)-gra*(n_man**2.0d0)/(he**(7.0d0/3.0d0))*ym(i)*abs(ym(i))*dt
        else
            ymn(i)=yms(i)
        end if
    end do
!
    return
    end subroutine visclc
!
!*******************************************************************************
    subroutine drhs00(ym, yms, ymn, yh, yzb)     ! for differential water depth : dh    ! ym:input, yms:before, ymn:after
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yms, yh, yzb
    real(k), dimension (0:im), intent(out) :: ymn
!
    integer i
!
    do i=1,nx
        if (yh(i).ge.hth.and.yh(i+1).ge.hth) then
            ymn(i)=yms(i)-(-(yh(i)+yzb(i))+(yh(i+1)+yzb(i+1)))*dt/dx*gra*(yh(i)+yh(i+1))*0.5d0
        else
            ymn(i)=0.0d0
        end if
    end do
!
    return
    end subroutine drhs00
!
!*******************************************************************************
    subroutine drhs01(ym, yms, ymn, yh, dh, yhn, yzb)     ! for differential water depth : dh ! ym:input, yms:before, ymn:after
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yms, yh, dh, yzb
    real(k), dimension (0:im), intent(out) :: ymn, yhn
!
    integer i
!
    do i=1,nx
        if (yh(i).ge.hth.and.yh(i+1).ge.hth) then
            ymn(i)=yms(i)-(-dh(i)+dh(i+1))*dt/dx*gra*(yh(i)+yh(i+1))*0.5d0
        else
            ymn(i)=0.0d0
        end if
    end do
!
    do i=2,nx
        yhn(i)=yh(i) - dt * ( (-ymn(i-1)+ymn(i))/dx )
    end do
!
    return
    end subroutine drhs01
!
!*******************************************************************************
    subroutine dpreseq(ym, yh, dhn)     ! for differential water depth : dh=h[n+1]-h[*]
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yh
    real(k), dimension (0:im), intent(out) :: dhn
!
    integer i, l, nn
    real(k) dx2, ap, aw, ae, ff, div
!
    real(k), dimension (:,:), allocatable :: a
    real(k), dimension (:), allocatable :: b, x
    integer allocatestatus, deallocatestatus
    integer n, n1, n2, m1, m2, itr, ier
    real(k) s, eps, ss
!
    real(k), dimension (:,:), allocatable :: a_tmp
    real(k), dimension (:), allocatable :: b_tmp
!
    n=nx
    m1=nx
    m2=m1+1
    n1=n+m2
    n2=-m2
    allocate (a(n2:n1,7),b(n),x(n2:n+m2), stat=allocatestatus)
    if (allocatestatus /= 0) stop "not enough memory"
    itr=itr_sor
    s=sigma_sor
    eps=err_sor
    a=0.0d0
    x=0.0d0
    b=0.0d0
!
    allocate (a_tmp(n2:n1,7),b_tmp(n), stat=allocatestatus)
    if (allocatestatus /= 0) stop "not enough memory"
    a_tmp=0.0d0
    b_tmp=0.0d0
!
    dx2=dx*dx
!
    do i=1,nx
        l=i
!
        aw=0.0d0
        ae=0.0d0
!
!-------------------------------------------------------------------------------
        if (i==1) then              ! west  dh(1)=dh(2)
            ap=-1.0d0/dx2*yh(i)
            ae= 1.0d0/dx2*yh(i)
            ff= 0.0d0
        else if (i==nx) then        ! east  dh(i)=0.0d0
            ap=-1.0d0/dx2*yh(i)
            ff= 0.0d0
        else if (i>=2.and.i<=(nx-1)) then
            div= (-ym(i-1)+ym(i))/dx
            if (yh(i)<=hth) then
                ap=-(1.0d0+1.0d0)/dx2*yh(i)
                ff=0.0d0
            else
                ap=-1.0d0/dx2*(yh(i-1)+2.0d0*yh(i)+yh(i+1))*0.5d0-1.0d0/(gra*dt*dt)
                aw=1.0d0/dx2*(yh(i-1)+yh(i))*0.5d0 !- (-yh(i-1)+yh(i+1))/yh(i)/dx2*0.25d0
                ae=1.0d0/dx2*(yh(i)+yh(i+1))*0.5d0 !+ (-yh(i-1)+yh(i+1))/yh(i)/dx2*0.25d0
                ff=div/dt/gra
            end if
        else
            write(*,*) "error on preseq2d!"
        end if
!-------------------------------------------------------------------------------
        a(l,1)=0.0d0
        a(l,2)=0.0d0
        a(l,3)=aw
        a(l,4)=ap
        a(l,5)=ae
        a(l,6)=0.0d0
        a(l,7)=0.0d0
        b(l)=ff
    end do
!
    i=1
    do l=1, n
!        write(30,'(2i5,8f18.6)') i, l, a(l,1), a(l,2), a(l,3), a(l,4), a(l,5), a(l,6), a(l,7), b(l)
!
        b_tmp(l)=b(l)
        do nn=1, 7
            a_tmp(l,nn)=a(l,nn)
        end do
!
!        if (i<nx) then
        i=i+1
!        end if
    end do
!
   call ilucgs9(a,n,n1,n2,m1,m2,b,eps,itr,s,x,ier)
!
    i=1
    do l=1, n
        dhn(i)=x(l)
!
        ss=0.0d0
        if (a_tmp(l,2)/=0.0d0) ss=ss+a_tmp(l,2)*x(l-nx)
        if (a_tmp(l,3)/=0.0d0) ss=ss+a_tmp(l,3)*x(l-1)
        if (a_tmp(l,4)/=0.0d0) ss=ss+a_tmp(l,4)*x(l)
        if (a_tmp(l,5)/=0.0d0) ss=ss+a_tmp(l,5)*x(l+1)
        if (a_tmp(l,6)/=0.0d0) ss=ss+a_tmp(l,6)*x(l+nx)
        if (abs(b_tmp(l)-ss)>1.0d-12) then
        write(*,*) 'dpreseq1:', i, l, b_tmp(l)-ss
        end if
!
        i=i+1
    end do
!
    deallocate (a, b, x, stat=deallocatestatus)
    if (deallocatestatus /= 0) stop "not open memory"
!
    if (eps>1.0d-12) then
        write(*,*) 'loop too much in ilucgs err =',eps
        write(*,*) 'ilucgs itteration =',itr
    end if
!
    return
    end subroutine dpreseq
!
!*******************************************************************************
      subroutine cipcls2_1d_r(f,r,u,istep)     ! cip-csl2, rational function
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(inout) :: f, r, u
    integer istep
!
    real(k) xmpm, dxmm, dxmkm, skm, alphakm, betakm, kappakm, kaikm, fm, rm, zz
    real(k), dimension (0:im) :: fn, rn, xx, dm
    integer i, j, km, kmup, kmup2, isg
    integer, dimension (0:im) :: kmi
!
    do i=0, im
        xx(i)=(real(i)-(real(im)-1.0d0)/2.0d0)*dx
        kmi(i)=-1
    end do
!
    do i=1,nx
        zz=1.0d0
        isg=-sign(zz,u(i))
!
        xmpm=xx(i)-u(i)*dt
!
        do j=1,nx
            if ((((xmpm-xx(j))*(xmpm-xx(j+isg)))<0.).or.(xmpm==xx(j))) kmi(i)=j
        end do
!
        km=kmi(i)
        kmup=km+isg
        kmup2=km-int((1-isg)/2)
!
        if (kmup<0.or.kmup>im.or.kmi(i)==-1) then
            dm(i) = 0.0d0
            fn(i) = f(i)
        else
            dxmm =xmpm-xx(km)
            dxmkm=xx(kmup)-xx(km)
!
            skm= real(isg)*r(kmup2)/dxmkm
            if ((f(kmup)-skm).ne.0.0d0) then
                betakm=( abs( (skm-f(km))/(f(kmup)-skm) ) - 1.0d0 )/dxmkm
                if (((skm-f(km))/(f(kmup)-skm))>=0.0d0) then
                    alphakm=1.0d0
                else
                    alphakm=0.0d0
                end if
            else
                betakm=-1.0/dxmkm
                if ((skm-f(km))>=0.0d0) then
                    alphakm=1.0d0
                else
                    alphakm=0.0d0
                end if
            end if
!
            kappakm = ( f(km) - skm + ( f(kmup)-skm )*( 1.0d0+alphakm*betakm*dxmkm ) )/(dxmkm*dxmkm)
            kaikm   = skm*alphakm*betakm + ( skm-f(km) )/dxmkm - kappakm*dxmkm
!
            dm(i) = ( ( kappakm*dxmm + kaikm )*dxmm + f(km) )*dxmm/( 1.0d0+alphakm*betakm*dxmm )
!
            fm    = ( ( 3.0d0*kappakm*dxmm + 2.0d0*kaikm )*dxmm + f(km) - alphakm*betakm*dm(i) )/( 1.0d0+alphakm*betakm*dxmm )
            fn(i) = ( 1.0d0 - (u(km)-u(kmup))/(xx(km)-xx(kmup))*dt )*fm
        end if
    end do
!
    do i=1,nx-1
        if (kmi(i)==-1.or.kmi(i+1)==-1) then
            rn(i) =r(i)
        else
            rm=0.0d0
            do j=0, im
                if ((j-kmi(i))*(j-(kmi(i+1)-1))<=0.) rm=rm+r(j)
            end do
!
           rn(i) = rm + ( dm(i+1)-dm(i) )
        end if
    end do
!
    do i=1,nx
        f(i) = fn(i)
    end do
!
    do i=1,nx-1
        r(i) = rn(i)
    end do
!
    return
    end subroutine cipcls2_1d_r
!
!*******************************************************************************
    subroutine cipsum(ysum,y)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: y
    real(k), dimension (0:im), intent(out) :: ysum
!
    integer i
!
    do i=0,nx
        ysum(i) = (y(i)+y(i+1))*0.5d0*dx
    end do
        ysum(im) = ysum(nx)
!
    return
    end subroutine cipsum
!
!*******************************************************************************
      subroutine ilucgs9(a,n,n1,n2,m1,m2,b,tol,itr,s,x,ier)
!*******************************************************************************
!  incomplete lu decomposition conjugate gradient squared method for   
!      the finite difference method.                                   
!                                                                      
!  parameters: same as ilubcg9 routine.                                
!                                                                      
!  copyright	tsutomu oguni      sep. 1 1997      ver. 1             
!*******************************************************************************
    use common_variables
    implicit none
    integer n, n1, n2, m1, m2, itr, ier
    real(k) tol, s
    real(k), dimension (n2:n1,7) :: a
    real(k), dimension (n) :: b
    real(k), dimension (n2:n+m2) :: x
!
    real(k), dimension (:), allocatable :: w, p, q, r, d, r0, e, h
    integer allocatestatus, deallocatestatus
    real(k) th, ss, c1, c2, c3, x1, x2, y, alpha, beta, res
    integer i, kk
!
    allocate (w(n2:n+m2), p(n2:n+m2), q(n2:n+m2), r(n2:n+m2), d(n2:n+m2), r0(n), e(n), h(n), stat=allocatestatus)
    if (allocatestatus /= 0) stop "not enough memory"
!
    ier = 0
!    if (n1 < n+m2 .or. m1 <= 1 .or. m1 >= m2 .or. m2 > n .or. n2 > -m2 .or. s < 0.0) then
    if (n1 < n+m2 .or. m1 <= 1 .or. m1 >= m2 .or. n2 > -m2 .or. s < 0.0) then
        write(*,*) '(subr. ilucgs9) invalid agument. ',n,n1,n2,m1,m2,s
        ier = 2
!        return
        go to 100
    end if
!
    th = 1.0d0
    if (s > 0.0 .and.s < 1.0) then
        th = s
        s = 1.0d0
    end if
    do i=1-m2,0
        d(i) = 0.0d0
        x(i) = 0.0d0
        p(i) = 0.0d0
        q(i) = 0.0d0
        r(i) = 0.0d0
        w(i) = 0.0d0
        w(i+n+m2) = 0.0d0
        p(i+n+m2) = 0.0d0
        q(i+n+m2) = 0.0d0
        r(i+n+m2) = 0.0d0
        x(i+n+m2) = 0.0d0
    end do
    d(1:n) = 0.0d0
!   incomplete cholesky decomposition
    if (s /= 0.0) then
        do i=1,n
            ss=s*a(i,4)-a(i,3)*(a(i-1,5)+(a(i-1,6)+a(i-1,7))*th)*d(i-1) &
                -a(i,2)*(a(i-m1,6)+(a(i-m1,5)+a(i-m1,7))*th)*d(i-m1) &
                -a(i,1)*(a(i-m2,7)+(a(i-m2,6)+a(i-m2,5))*th)*d(i-m2)
            d(i) = 1.0d0 / ss
        end do
    else
        do i=1,n
            ss = a(i,4) - a(i,3)*a(i-1,5)*d(i-1) - a(i,2)*a(i-m1,6)*d(i-m1) - a(i,1)*a(i-m2,7)*d(i-m2)
            d(i) = 1.0d0 / ss
        end do
    end if
    do i=1,n
       q(i) = a(i,1)*x(i-m2)+a(i,2)*x(i-m1)+a(i,3)*x(i-1)+a(i,4)*x(i)+a(i,5)*x(i+1)+a(i,6)*x(i+m1)+a(i,7)*x(i+m2)
    end do
    r(1:n) = b(1:n) - q(1:n)
    do i=1,n
        r(i) = d(i)*(r(i)-a(i,3)*r(i-1)-a(i,2)*r(i-m1)-a(i,1)*r(i-m2))
    end do
    do i=n,1,-1
        r(i) = r(i)-d(i)*(a(i,5)*r(i+1)+a(i,6)*r(i+m1)+a(i,7)*r(i+m2))
    end do
    do i=1,n
        r0(i) = r(i)
        p(i) = r(i)
        e(i) = r(i)
    end do
    c1 = dot_product(r(1:n),r(1:n))
!   iteration phase
!    do k=1,itr
    do kk=1,itr
        do i=1,n
            q(i) = a(i,1)*p(i-m2)+a(i,2)*p(i-m1)+a(i,3)*p(i-1)+a(i,4)*p(i)+a(i,5)*p(i+1)+a(i,6)*p(i+m1)+a(i,7)*p(i+m2)
        end do
        do i=1,n
            q(i)=d(i)*(q(i)-a(i,3)*q(i-1)-a(i,2)*q(i-m1)-a(i,1)*q(i-m2))
        end do
        do i=n,1,-1
            q(i)=q(i)-d(i)*(a(i,5)*q(i+1)+a(i,6)*q(i+m1)+a(i,7)*q(i+m2))
        end do
        c2 = dot_product(q(1:n),r0(1:n))
        if (c2 == 0.0) then
                ier = 5
                tol = res
!                itr = k
                itr = kk
!                return
            go to 100
        end if
        alpha = c1 / c2
        c3 = 0.0d0
        x1 = 0.0d0
        x2 = 0.0d0
        h(1:n) = e(1:n) - alpha * q(1:n)
        w(1:n) = e(1:n) + h(1:n)
        do i=1,n
        q(i) = a(i,1)*w(i-m2)+a(i,2)*w(i-m1)+a(i,3)*w(i-1)+a(i,4)*w(i)+a(i,5)*w(i+1)+a(i,6)*w(i+m1)+a(i,7)*w(i+m2)
        end do
        do i=1,n
            q(i)=d(i)*(q(i)-a(i,3)*q(i-1)-a(i,2)*q(i-m1)-a(i,1)*q(i-m2))
        end do
        do i=n,1,-1
            q(i)=q(i)-d(i)*(a(i,5)*q(i+1)+a(i,6)*q(i+m1)+a(i,7)*q(i+m2))
        end do
        do i=1,n
            y = x(i)
            r(i) = r(i) - alpha * q(i)
            x(i) = x(i) + alpha * w(i)
            c3 = c3 + r(i) * r0(i)
            x1 = x1 + y * y
            x2 = x2 + (x(i) - y)**2
        end do
        if (x1 /= 0.0) then
            res = dsqrt(x2 / x1) 
            if (res <= tol) then
!                itr = k
                itr = kk
                ier = 0
                tol = res
                if (th /= 1.0d0) s = th
!                return
                go to 100
            end if
        end if 
        beta = c3 / c1
        c1 = c3
        do i=1,n
            e(i) = r(i) + beta * h(i)
            p(i) = e(i) + beta * (h(i) + beta * p(i))
        end do
!
    end do
    ier = 1
    write(*,*) '(subr. ilucgs9) no convergence. '

    stop
    
    tol = res
    if (th /= 1.0d0) s = th
!    return
    go to 100
!
100 deallocate (w, p, q, r, d, r0, e, h, stat=deallocatestatus)
    if (deallocatestatus /= 0) stop "not open memory"
    return
!
    end subroutine ilucgs9
!
!*******************************************************************************
    subroutine bedstress(taustb, tauste, taustbc, taustec, ym, yu, yh)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: ym, yu, yh
    real(k), dimension (0:im), intent(out) :: taustb, tauste, taustbc, taustec
!
    integer i
    real(k) ui, hi
!
    do i=0,im
        if (yh(i).ge.hth) then
            taustbc(i)=(n_man**2.0d0)*yu(i)*abs(yu(i))/(sss*ds*yh(i)**(1.0d0/3.0d0))
            if (yu(i).gt.0.d0) then
                taustec(i)= ((yu(i)/(6.0d0+2.5d0*log(yh(i)/ds/(1.0d0+2.0d0*taustbc(i)))))**2.0d0)/(sss*gra*ds)
            else if (yu(i).lt.0.d0) then
                taustec(i)=-((yu(i)/(6.0d0+2.5d0*log(yh(i)/ds/(1.0d0-2.0d0*taustbc(i)))))**2.0d0)/(sss*gra*ds)
            else
                taustec(i)=0.0d0
            end if
        else
            taustbc(i)=0.0d0
            taustec(i)=0.0d0
        end if
    end do
!
    do i=0,nx
        hi=(yh(i)+yh(i+1))*0.5d0
        if (hi.ge.hth) then
            ui=ym(i)/hi
            taustb(i)=(n_man**2.0d0)*(ui)*abs(ui)/(sss*ds*hi**(1.0d0/3.0d0))
            if (ui.gt.0.d0) then
                tauste(i)= ((ui/(6.0d0+2.5d0*log(hi/ds/(1.0d0+2.0d0*taustb(i)))))**2.0d0)/(sss*gra*ds)
            else if (ui.lt.0.d0) then
                tauste(i)=-((ui/(6.0d0+2.5d0*log(hi/ds/(1.0d0-2.0d0*taustb(i)))))**2.0d0)/(sss*gra*ds)
            else
                tauste(i)=0.0d0
            end if
        else
            taustb(i)=0.0d0
            tauste(i)=0.0d0
        end if
    end do
        taustb(im)=taustbc(im)
        tauste(im)=taustec(im)
!
    return
    end subroutine bedstress
!
!*******************************************************************************
    subroutine iwagakieq(taustc)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), intent(out)     :: taustc
!
    real(k) ustc
!
    if(ds*100.d0.ge.0.3030d0) then
        ustc=sqrt(80.9d0*ds*100.d0)*0.01d0
    elseif(ds*100.d0.ge.0.1180d0) then
        ustc=sqrt(134.6d0*(ds*100.d0)**(31.d0/22.d0))*0.01d0
    elseif(ds*100.d0.ge.0.0565d0) then
        ustc=sqrt(55.0d0*ds*100.d0)*0.01d0
    elseif(ds*100.d0.ge.0.0065d0) then
        ustc=sqrt(8.41d0*(ds*100.d0)**(11.d0/32.d0))*0.01d0
    else
        ustc=sqrt(226.d0*ds*100.d0)*0.01d0
    end if
!
    taustc=ustc**2.0d0/(sss*gra*ds)
!
    return
    end subroutine iwagakieq
!
!*******************************************************************************
    subroutine bedload_am(qb, taustb, tauste, taustc, yzb) ! ashida-michiue
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: taustb, tauste, yzb
    real(k), dimension (0:im), intent(out) :: qb
    real(k), intent(in)      :: taustc
!
    integer i
    real(k) qblimit
!
    qb=0.d0
!
    do i=0,nx
        if (abs(taustb(i)).le.taustc) cycle
!
        if (taustb(i).gt.0.d0) then
            qb(i)= 17.d0*abs(tauste(i))**(3.d0/2.d0)*(1.0d0-taustc/abs(taustb(i)))*(1.0d0-sqrt(taustc/abs(taustb(i))))*sqrt(sss*gra*ds**3.d0)
!
            qblimit=(1.d0-lambda)*dx/dt*(yzb(i)-yzb(i+1))*0.5d0
            if ((qb(i).gt.qblimit).and.(qblimit.ge.0.d0)) qb(i)=qblimit
        else if (taustb(i).lt.0.d0) then
            qb(i)=-17.d0*abs(tauste(i))**(3.d0/2.d0)*(1.0d0-taustc/abs(taustb(i)))*(1.0d0-sqrt(taustc/abs(taustb(i))))*sqrt(sss*gra*ds**3.d0)
!
            qblimit=(1.d0-lambda)*dx/dt*(yzb(i)-yzb(i+1))*0.5d0
            if ((qb(i).lt.qblimit).and.(qblimit.le.0.d0)) qb(i)=qblimit
        else
            qb(i)=0.d0
        end if
    end do
        qb(0)=qb(1)
        qb(im)=qb(nx)
!
    return
    end subroutine bedload_am
!
!*******************************************************************************
    subroutine bedload_kp(qb, taustb, tauste, taustc, yzb) ! kovacs-parker
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: taustb, tauste, yzb
    real(k), dimension (0:im), intent(out) :: qb
    real(k), intent(in)      :: taustc
!
    integer i
    real(k) kc, mus, tantx, qblimit
!
    mus=1.d0
    qb=0.d0
!
    do i=1,nx
        if (abs(taustb(i)).le.taustc) cycle
!
        if (taustb(i).gt.0.d0) then
            tantx=(-yzb(i)+yzb(i+1))/dx
            mus=-tan_repose
            kc=1.d0+tantx/mus
            if (kc.gt.1.9999d0) then
                kc=1.9999d0
            else if (kc.le.0.0001d0) then
                kc=0.0001d0
            end if
            qb(i)= 17.d0/kc*tauste(i)**(3.d0/2.d0)*(1.0d0-kc*taustc/taustb(i))*(1.0d0-sqrt(kc*taustc/taustb(i)))*sqrt(sss*gra*ds**3.d0)
!
            qblimit=(1.d0-lambda)*dx/dt*(yzb(i)-yzb(i+1))*0.5d0
            if ((qb(i).gt.qblimit).and.(qblimit.ge.0.d0)) qb(i)=qblimit
        else if (taustb(i).lt.0.d0) then
            tantx=-(-yzb(i)+yzb(i+1))/dx
            mus=tan_repose
            kc=1.d0+tantx/mus
            if (kc.gt.1.9999d0) then
                kc=1.9999d0
            else if (kc.le.0.0001d0) then
                kc=0.0001d0
            end if
            qb(i)=-17.d0/kc*abs(tauste(i))**(3.d0/2.d0)*(1.0d0-kc*taustc/abs(taustb(i)))*(1.0d0-sqrt(kc*taustc/abs(taustb(i))))*sqrt(sss*gra*ds**3.d0)
!
            qblimit=(1.d0-lambda)*dx/dt*(yzb(i)-yzb(i+1))*0.5d0
            if ((qb(i).lt.qblimit).and.(qblimit.le.0.d0)) qb(i)=qblimit
        end if
    end do
        qb(0)=qb(1)
        qb(im)=qb(nx)
!
    return
    end subroutine bedload_kp
!
!*******************************************************************************
    subroutine suspendedload(qs, taustbc, taustec, taustc)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: taustbc, taustec
    real(k), dimension (0:im), intent(out) :: qs
    real(k),      intent(in) :: taustc
!
    real(k) qsus
    integer i, icon, iten
!
    qs=0.0d0
!
    do i=0,im
!
        if (taustbc(i).lt.taustc) cycle
!
        call itakura_kishi(1.0d0,ds,taustec(i),taustbc(i),taustc,qsus,icon,iten)
!
        qs(i)=qsus
!
        if ( icon /=0 .or. qsus < 0.d0 ) then
            write(*,*) i,qsus,taustbc(i),taustc
            stop
        end if
    end do
!
    return
    end subroutine suspendedload
!
!*******************************************************************************
    subroutine itakura_kishi(fraction,diameter,taust_effective,taust,taust_c,qs,icon,iten)
!*******************************************************************************
    use common_variables
    implicit none
!
!    real,     intent(in)  :: fraction
    real(k),     intent(in)  :: fraction
    real(k),     intent(in)  :: diameter
    real(k),     intent(in)  :: taust_effective, taust, taust_c
    real(k),     intent(out) :: qs
    integer,     intent(out) :: icon
    integer,     intent(in)  :: iten
!
!
    real(k)                  :: conk,alpst,bst0
    real(k)                  :: bst
    real(k)                  :: taust_cdi
    real(k)                  :: omega
    real(k)                  :: ss
    real(k)                  :: ws
!
    data conk,alpst,bst0/0.008d0,0.14d0,0.143d0/
!
    icon =0
    qs=0.0
!    if(taust < taust_c) return
!
    ss=(dsand-dwat)/dsand
!
    bst=bst0
!
    call omega_for_itakura_kishi(taust_effective,bst,omega)
!
    call rubey(diameter,ws)
!
    qs=fraction*conk*(alpst*ss*gra*diameter                           &
               /sqrt(taust_effective*sss*gra*diameter)*omega          &
               -ws)
!
    qs = max(qs, 0.d0)
!    if ( qs < 0.0 ) icon=-1
!
    return
    end subroutine itakura_kishi
!
!*******************************************************************************
    subroutine omega_for_itakura_kishi(taust_effective,bst,omega)
!*******************************************************************************
    use common_variables
    implicit none
!
    real(k),     intent(in)  :: taust_effective
    real(k),     intent(in)  :: bst
    real(k),     intent(out) :: omega
!
    real(k)                  :: a1,a2,a3
    real(k)                  :: eta0
    real(k)                  :: sqrt2,x,bb,aa
    real(k)                  :: numerator,denominator
    real(k)                  :: a_prime
    real(k)                  :: tmp1, tmp2, gsi, dgsi
    integer                  :: count
!
    data a1,a2,a3/0.4361836d0,-0.1201676d0,0.937298d0/
    data eta0/0.5d0/
!
    sqrt2=sqrt(2.d0)
!
    a_prime=bst/taust_effective-1.d0/eta0
!
    if(a_prime.ge.0.) x=a_prime*sqrt2
    if(a_prime.lt.0.) x=-a_prime*sqrt2
!
    bb=1.d0/(1.d0+0.33627d0*x)
    aa=a1*bb+a2*bb**2+a3*bb**3
    numerator=1./sqrt(pi)*exp(-a_prime**2)
    if(a_prime.ge.0.d0) denominator=numerator*aa/sqrt2
    if(a_prime.lt.0.d0) denominator=1.d0-numerator*aa/sqrt2
    numerator=0.5d0*numerator
!
    if ((taust_effective.eq.0.d0).or.(denominator.eq.0.d0)) then
        omega=0.d0
    else
        omega=taust_effective/bst                         &
              *numerator                                  &
              /denominator                                &
              +taust_effective/bst/eta0                   &
              -1.d0
    end if
!
    return
    end subroutine omega_for_itakura_kishi
!
!*******************************************************************************
    subroutine cal_wscbot(cb,dep,taustbc,wscbot)
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: cb, dep, taustbc
    real(k), dimension (0:im), intent(out) :: wscbot
!
    integer    :: i
    real(k)    :: ws, cbot, taubc
!
    wscbot=0.d0
!
    do i=0, im
        if (dep(i).lt.hth) cycle
!
        call rubey(ds,ws)
!
        taubc=dwat*sss*gra*ds*taustbc(i)

        call bottom_concentration1(cb(i),taubc,dep(i),ds,cbot)
!
        wscbot(i)=ws*cbot
!
    end do
!
    return
    end subroutine cal_wscbot
!
!*******************************************************************************
    subroutine bottom_concentration1(cb,tau,dep,diameter,cbot)
!*******************************************************************************
    use common_variables
    implicit none
    real(k),     intent(in)  :: cb, tau, dep, diameter
    real(k),     intent(out) :: cbot
!
    real(k)                  :: beta, keddy, ws, alf

    cbot=cb
!
    alf=0.d0
!
    if (cb<=0.d0) return
!      
    if( tau < taumin ) return
!
    keddy=kappa*sqrt(tau/dwat)*dep/6.d0
    call rubey(diameter,ws)
    beta=ws*dep/keddy
!
    if ( beta > 20.0d0 ) then
        alf=beta
    else if ( beta <= 0.d0 ) then
        alf=0.0d0
    else if ( beta < 1.d-6 ) then
        alf=1.0d0
    else
        alf=beta/(1.0d0-exp(-beta))
    end if
!
    cbot=cb*alf
!
    return
    end subroutine bottom_concentration1
!
!*******************************************************************************
    subroutine rubey(diameter,ws)
!*******************************************************************************
    use common_variables
    implicit none
!
    real(k),     intent(in)  :: diameter
    real(k),     intent(out) :: ws
!
!    if(diameter.gt.0.001) then
!        ws=32.8d0*sqrt(diameter*100.)*0.01
!    else
        ws=sqrt(2.0d0/3.0d0*sss*gra*diameter+(6.0d0*nuw/diameter)**2.0d0)-6.0d0*nuw/diameter
!    end if
!
    return
    end subroutine rubey
!
!*******************************************************************************
    subroutine calbed(yzb, yzbs, yzbn, qb, qs, wscbot)   ! yzb:input, yzbs:before, yzbn:after
!*******************************************************************************
    use common_variables
    implicit none
    real(k), dimension (0:im), intent(in)  :: yzb, yzbs, qb, qs, wscbot
    real(k), dimension (0:im), intent(out) :: yzbn
!
    integer i
!
    do i=1,nx
        yzbn(i)=yzbs(i)-( (-qb(i-1)+qb(i))/dx+qs(i)-wscbot(i))/(1.d0-lambda)*dt
    end do
!
    return
    end subroutine calbed
!
!*******************************************************************************
