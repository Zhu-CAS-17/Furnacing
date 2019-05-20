c all variables are in international unit
c change radi_c
c change tion id_ele
c conduc
c rmu_c
c th_conduc
c matr_ele
c add_con
c tine grid_mag
c l:746 set reflection boundary condition
c l:3819, 602 set no radiation in chamber
c l:3420 set no radiation in SiC powder
c l:1628 set coil position
c tine bound: set radiation source term and latent heat source term

c radi_ini: coil_xx0,coil_yy0
c need change th_conduc, heat_cap, conduc, rho_cal
c parameters:
c cden crys_len pres_sys coil_xx0 z_pres
c cal_ex1 for gr_rate_r in sterling and imech systems

      program main
	include 'sc.h'
c for sticking model
c      call sticking
c	call cal_gr_advan
c for calculation of growth rate for JCG paper
c gr_2300k10k gr_2300k100k need all the following subroutines
c      call cal_gr_ster
c      call cal_gr_ster_pres
c      call cal_gr_ster_gradient
c fixed pressure
c      call cal_gr_ster_gradient_t
c      stop


c	open(1,file='alfa_gra')
c	write(*,*)'Temperature (K), conductivity (W/m/K), 
c     1    diffusivity (m2/s)'
c	do tkel=250,3000,250
c	alfa=cal_diffu(tkel)
c	enddo
c      stop

      open(unit=11, file='line1')
      rewind 11
      open(unit=12, file='tmax_t.dat')
      rewind 12

c initial field choice: init=0 freezing initial
c                           =1 restart initial

      init=1
c interstep iteration number
c EACH TIME STEP NEED TO BE CONVERGENCE. 
      last=10
	iter=0
	nfmax=7
c timestep number

c For test purpose, please change lastt=2. 

      lastt=1


	crys_len=0.00
        ndt=0

c e.g. 'input.dat': 1 0.05 100
	open(1,file='input.dat')
	read(1,*)init,dtm,lastt
      read(1,*)idesign,pres_sys,cden0,crys_len
      read(1,*)coil_move,felt_incr

      read(1,*)icontrol
      read(1,*)t_ref
      p_int=0

        read(1,*)(relax(i),i=1,6)
        read(1,*)(nsolve(i),i=1,6)
        read(1,*)(ntimes(i),i=1,6)

	close(1)
      cden=cden0

c in old version, pres_sys is in Torr, we changed it to Pa
c      pres_sys=pres_sys*133.3

	darcy=1d10

c for Sterling furnace
        open(11,file='geom.dat')
        do i=1,10
        read(11,*)
        enddo
        read(11,*)id_total
        do i=1,id_total
        read(11,*)
        read(11,*)xx_0(i),xx_1(i),yy_0(i),yy_1(i),id_el0(i)
        if(id_el0(i).eq.10)read(11,*)air_len0(i)
        enddo
        close(11)
	open(99,file='grid.dat')
        read(99,*)
      read(99,*)nxx
      do i=1,nxx
      read(99,*)iddx(i),xxd(i)
      enddo
        read(99,*)
      read(99,*)nyy
      do i=1,nyy
      read(99,*)iddy(i),yyd(i)
      enddo
	close(99)

        freq=10600

        open(1,file='coil.dat')
        read(1,*)freq
        read(1,*)coil_xx0,coil_yy0,coil_dxx,width_coil,xx_r,yy_r
        close(1)

        open(1,file='radi.dat')
        read(1,*)
        read(1,*)rr_c_t,zb_c_t,zt_c_t
        read(1,*)
        read(1,*)rr_c_l,zb_c_l,zt_c_l
        read(1,*)
        read(1,*)rc,zb0,zt0,z_pres,z_pres2
        close(1)

	open(1,file='latent.dat')
	read(1,*)
	read(1,*)gr_lat,zb_lat,zt_lat
	close(1)

        omega=freq*2*3.14159


        icom_mag=10
        if(itest.eq.1)icom_mag=50
        iout_dat=20


c coordinate choice: mode=0 rectangular
c                        =1 cylinder coordinate
      mode=1


c      call grid
c      call start

	write(*,*)'for general'
c generate xc, yc
      if(idesign.eq.1)then
      call grid_mag_design2
      else if(idesign.eq.2)then
	call grid_mag_design2
	else if(idesign.eq.3)then
	call grid_mag_design3
        else if(idesign.eq.4)then
        call grid_mag_design4
        else if(idesign.eq.5)then
        call grid_mag_design5
        else if(idesign.eq.6)then
        call grid_mag_design4
	endif
c for x,y,hksi,heta,ak1,ae1...
	call setup1
c call restart to read u, v, w
	write(*,*)'start'
	call start
c initiate computational domain for radition, used in gamsor and boun
	write(*,*)'radi_ini'
	call radi_ini 
      tend=tinit+dtm*dfloat(lastt)
      tpass=tinit+dtm*ndt
      print 18, tinit
	print 19, tend
18	format(8x,' start time of process : ',e10.4,
     1    'nondimensional time')
19	format(8x,'   end time of process : ',e10.4,
     1    'nondimensional time')

c --- time loop starts --------------

c lsolve: for u, v, p, pc, t
c	do 41 k=1,4
c41	lsolve(k)=.true.
      lsolve(5)=.true.
c	lsolve(6)=.true.
	do i=1,l1
	do j=1,m1
	if(id_ele(i,j).eq.0)then
c	write(*,*)'less than 0',i,j,t(i,j)
	t(i,j)=293
	endif
	enddo
	enddo
c loop
10	continue



      do i=1,l1
      if(x(i,2).gt.zb_c_t)then
      t_top=t(i-1,2)
      goto 91
      endif
      enddo
91    continue

      do i=l1,1,-1
      if(x(i,2).lt.zb_c_l)then
      t_bot=t(i+1,2)
      goto 93
      endif
      enddo
93    continue


        if(mod(ndt,icom_mag).eq.0)then
	goto 90
	open(11,file='mag.dat')
	read(11,*)l1o,m1o
	if(l1o.ne.l1.or.m1o.ne.m1)stop 'err in mag.dat'
	read(11,88)((qth(i,j),i=1,l1),j=1,m1)
	read(11,88)((fc(i,j),i=1,l1),j=1,m1)
	close(11)
	do i=1,l1
	do j=1,m1
	fco(i,j)=fc(i,j)
	enddo
	enddo
88	format(6(1pe12.5,1x))
      goto 102
90	continue


      if(icontrol.ne.0)then
      if(icontrol.eq.1)then
      cden=cden0+(t_ref-t_top)*400/483
      endif
      if(icontrol.eq.2)then
c PI control
c Tau=5000s
      p_int=p_int+(t_ref-t_top)*dtm
      cden=cden0+0.6*(t_ref-t_top)*400/483+1/(0.5*5000)*p_int
      endif
      write(*,*)'cden changed to',cden,'t_top',t_top
      endif

	call cal_elemag
102    continue
	endif

	call varmax
      lconv=.true.
	do i=1,nfmax
	lcon(i)=.false.
	enddo

      call dense
	write(*,*)'setup1'
c for x,y,hksi,heta,ak1,ae1...
      call setup1

c for ruksi,rueta

      if(tpass.eq.0.0) call initl
	
	write(*,*)'setup2'
c get source term
	do i=1,l1
	do j=1,m1
	pc(i,j)=0
	enddo
	enddo

20	continue
	write(*,*)'ndt',ndt,'tpass=',tpass
	write(*,*)'iter=',iter
	if(iter.eq.0)then
	write(*,*)lcon(6)
	write(*,*)ccen(i_rt,j_rr-1)
	endif
	call cal_ex1
      call setup2
	if(iter.eq.0)then
	write(*,*)'con_rb',ccen(i_rb-1,2),ccen(i_rb,2)
	write(*,*)'ccen',ccen(i_rt,j_rr-1)
	write(*,*)'stream'
	write(*,*)'p1=',p(1,1),p(1,2),p(1,3)
	endif
      call stream
	if(iter.eq.0)then
	write(*,*)'output'
	write(*,*)'p2=',p(1,1),p(1,2),p(1,3)
	write(*,*)'p3=',p(1,1),p(1,2),p(1,3)
	endif
      iter=iter+1
      if(iter.gt.6.and.lconv) go to 60
        if(itest.eq.1)goto 60
c	if(iter.ge.2)goto 60
      if(iter.le.last) go to 20


      print 42,last
   42 format(/5x,'not convergence after',i4,' iterations'//)
60	continue
      ndt=ndt+1
      tpass=tpass+dtm
	if(tpass.ge.tend) lstop=.true.
	call output
      if(lstop) go to 80

      write(11,62)tpass,fluxm1
c      write(12,62)tpass,fluxm2
   62 format(4x,f12.4,3x,f12.4)
c not move grid
c      call movegrid

c update fo values of u, v, t, ccen
      call newstep
c      call movegrid
      go to 10
   80 continue

      write(*,*)'temperature at z=0.1 m'

	xi0=0.1
	i0=interid(xi0,x,l1)
	do j=1,m1
	write(*,*)y(i0,j),t(i0,j)
	enddo

      write(*,*)'temperature at center'
        open(22,file='profvt.dat')
        do i=1,l1
        write(*,*)x(i,1),t(i,1)
	if(t(i,1).ne.293)then
        write(22,*)x(i,1),t(i,1)
	endif
        enddo
        close(22)

        write(*,*)'temperature at z=0.178m (design2 4)'
        write(*,*)'or z=0.175m (design3)'
        open(22,file='profht.dat')
	if(idesign.eq.3)then
        xi0=0.175
      else if(idesign.eq.1.or.idesign.eq.2)then
	xi0=0.178
        else if(idesign.eq.4)then
	xi0=0.178
        else if(idesign.eq.5)then
	xi0=0.178
        else if(idesign.eq.6)then
      xi0=0.1765

	endif
        i0=interid(xi0,x,l1)
        do j=1,m1
        write(*,*)y(i0,j),t(i0,j)
	if(t(i0,j).ne.293)then
        write(22,*)y(i0,j),t(i0,j)
	endif
        enddo
        close(22)
        open(1,file='cham_t.dat')
        write(1,*)'Top and bottom temperature inside growth chamber'
        write(1,*)t(i_rt+1,1),t(i_rb-1,1),tmax
        close(1)

      stop
      end



	function rho_cal(i,j)
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: griphe foam
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11:copper 12: bottom graphite 13: crystal
c 14: gas with apparent conductivity 50 W/m/K

	include 'sc.h'
	dimension rho_d(12)
	data rho_d/2200,2200,1896,0.448,2200,2200,2200,2200,1900,0.448,
     1    8960,950/

c graphite density : 2200 kg/m3
c SiC density: 3160 kg/m3, porosity: 0.4
c air density: 0.448 kg/m3, 1atm and 800C
	rho_cal=0.448
	id=id_ele(i,j)
	if(id.ne.0.and.id.le.12)then
	rho_cal=rho_d(id)
	endif

c	if(id.eq.3.or.id.eq.4)rho_cal=rho_cal+ccen(i,j)*0.040097

        if(id.eq.0.or.id.eq.4.or.id.eq.10.or.id.eq.14)then
	tt0=t(i,j)
	if(tt0.le.293)tt0=293
	rho_cal=pres_sys*0.039948/8.3145/tt0
	endif

	if(id.eq.13)rho_cal=3160

	return
	end

	function gam_cal(i,j)
	include 'sc.h'
	real*8 gam_mat(4)
	data gam_mat/20,5,50,0.05/
	id=id_mat(i,j)
c for  conductivity
c 1: graphite 2: SiC 3: sintered SiC 4: Argon
	if(nf.eq.5.and.id.ge.1.and.id.le.4)then
	gam_cal=gam_mat(id)
	endif

	return
	end

	function id_mat(i,j)
	include "sc.h"
c 1: graphit
c 2: SiC powder
c 3: Sic Xtal
c 4: air
	xxx=x(i,j)
	yyy=y(i,j)
	id_mat=1
	if((xxx.ge.0.5*0.0254.and.xxx.le.3*0.0254)
     1   .and.yyy.le.1.25*0.0254)then
	id_mat=2
	endif
	x_xtal=0.01
	if((xxx.ge.(5.5-x_xtal)*0.0254
     1   .and.xxx.le.5.5*0.0254)
     1  .and.yyy.le.1*0.0254)then
	id_mat=3
	endif
	if((xxx.ge.3*0.0254.and.xxx.le.(5.5-x_xtal)*0.0254)
     1   .and.yyy.le.1*0.0254)then
	id_mat=4
	endif
	return
	end

      subroutine user
	include 'sc.h'

   95 format(e12.5)
   96 format(6e12.5)
   88 format(6(1pe12.5,1x))
   87 format(1x,i2)
   86 format(1x,'title=','"zhanghui"')
   84 format(1x,'variables ="x","y","u","v","T(K)",
     1    "concentration","partial pressure","parti","st"')
   83 format(1x,'zone t= "zone 1" , i=',i3,' , j=',i3,' , f= block')
c  t(solid)=t(melting): lst=m1;lsolve(ns)=.false. zhanghui 2/28/91
c  if steady rdtm=0.0 
	

      entry grid

      mode=0
      l0=72
      m0=68
      lst3=55
      lst=m0/2

      xlen=(0.5+2.5+2.5+7./16.)*0.0254
      ylen=3.25/2.*0.0254

      l1=l0
      m1=m0
      l2=l1-1
      m2=m1-1

      ib=9
      bx(1)=0.0  
      by(1)=0.0
      bx(2)=0.015  
      by(2)=0.3
      bx(3)=0.03
      by(3)=0.5
      bx(4)=0.046  
      by(4)=0.6
      bx(5)=0.07  
      by(5)=0.7
      bx(6)=0.10
      by(6)=0.8
      bx(7)=0.14
      by(7)=0.9
      bx(8)=0.18
      by(8)=0.95
      bx(9)=0.25  
      by(9)=1.0

      inum=9
      do 40 k=1,inum
      xin(k)=by(k)
   40 yin(k)=bx(k)

      if(init.eq.0)call initgrid

      return

c --------------------------------------------------------------------

      entry start

      lstop=.false.
      lortho=.false.
c 1:u 2: v 3:p 4:pc 5:t 6:c

c      ntimes(1)=5
c      ntimes(2)=5
c      ntimes(5)=6
c      ntimes(6)=3
c      nsolve(1)=3
c      nsolve(2)=3
c      nsolve(3)=3
c      nsolve(4)=3
c      nsolve(5)=5
c      nsolve(6)=10
c        if(itest.eq.1)then
c        ntimes(5)=2
c        nsolve(5)=2
c        endif
c for u, v, p, pc, t
c      relax(1)=0.05
c      relax(2)=0.05
c      relax(3)=1
c        relax(4)=1
c
c      relax(5)=1.
c      if(itest.eq.1)relax(5)=1
c      relax(6)=1.

c --- set/read initial variable field -------

      re=1.0
      pr=0.45
      grash=70

      ra=grash*pr
      rre=1.0/re
      rpr=1.0/pr

c graphite heat capacity
      rdtm=1./dtm
	write(*,*)'dtm=',dtm
c	rdtm=0

      stes=1.32d-2
      stel=5.1d-3

      fksl=1.0
      frsl=1.0
      fcsl=1.0

      do 5 i=1,l1
      do 5 j=1,m1
	nf=5
	id_bound(i,j)=id_bound_f(i,j)
      rho(i,j)=rho_cal(i,j)
      gam(i,j)=th_conduc(i,j)
      con(i,j)=0.0
	a_mag(1,i,j)=0
	a_mag(2,i,j)=0
	void(i,j)=1
	if(id_ele(i,j).eq.3)void(i,j)=0.4
	ap(i,j)=0.0
	apc(i,j)=0.0
	fc(i,j)=0
5	continue

      do 6 i=1,l1
      do 6 j=1,m1
      do k=1,nfmax
      f(i,j,k)=0.0
	enddo
	t(i,j)=293
	id=id_ele(i,j)
c	if(id.eq.3)t(i,j)=2260
c	if(id.eq.4)t(i,j)=2260
      i_powder(i,j)=0

6	continue


      if(init.ne.0) call restart

      do 11 i=1,l1
      do 11 j=1,m1
	ruksi(i,j)=0
	rueta(i,j)=0
      do k=1,nfmax
c      if(nf.eq.3.or.nf.eq.4) go to 11
      fo(i,j,k)=f(i,j,k)
	fn(i,j,k)=f(i,j,k)
	enddo
	rho_o(i,j)=rho(i,j)
   11 continue

c for magnetic potential
	do i=1,l1
	do j=1,m1
	fco(i,j)=fc(i,j)
	enddo
	enddo
      return
	entry varmax

	tmin=1d6
	tmax=0
	do i=1,l1
	do j=1,m1
	if(tmin.gt.t(i,j))tmin=t(i,j)
	if(tmax.lt.t(i,j))tmax=t(i,j)
	enddo
	enddo
	write(*,*)'tmax=',tmax,'tmin=',tmin
	return

      entry dense
      return
	end


	subroutine bound
	include 'sc.h'
c set bi

	write(*,*)nf,'u=',u(i_rt,2),u(i_rt+1,2),ruksi(i_rt+1,2),
     1    ruksi(i_rt,2)

      do  i=1,l1
      ruksi(i,1)=0.0
	ruksi(i,m1)=0.0
	enddo
      do j=1,m1
      rueta(1,j)=0.0
	rueta(l1,j)=0.0
	enddo

      do i=1,l1
      rueta(i,2)=0.0
	rueta(i,m1)=0.0
	enddo
      do j=1,m1
      ruksi(2,j)=0.0
	ruksi(l1,j)=0.0
	enddo

	if(nf.ne.1)goto 20
	call boun_b
	do i=1,l1
	gam(i,1)=0.0
	ajm(i,2)=0
	f(i,1,nf)=f(i,2,nf)
	do j=2,m2
	id=id_bound(i,j)
	if(id.ne.0)then
	ajm(i,j)=0
	goto 105
	endif
	enddo
105	continue
	enddo
	mm=j_rr
	do j=2,mm-1
c top
	i=i_rt
	u(i+1,j)=growthsp(j)
	enddo
	write(*,*)nf,'ugr=',u(i_rt,2),u(i_rt+1,2),ruksi(i_rt+1,2),
     1    ruksi(i_rt,2)
20	continue
	if(nf.ne.2)goto 300
	call boun_b
300	continue
	if(nf.ne.3)goto 400
	call boun_b
	do i=1,l1
	gam(i,1)=0.0
	ajm(i,2)=0
	f(i,1,nf)=f(i,2,nf)
	enddo
400	continue
	if(nf.ne.4)goto 500
	call boun_b
	do i=1,l1
	gam(i,1)=0.0
	ajm(i,2)=0
	f(i,1,nf)=f(i,2,nf)
	enddo
500	continue
	if(nf.ne.5)goto 600
	qmax=0
	volt=0
	do i=2,l2
	do j=2,m2
	con(i,j)=con(i,j)+qth(i,j)*vol(i,j)
	if(id_bound(i,j).ne.0)then
	qmax=qmax+qth(i,j)*vol(i,j)
	volt=volt+vol(i,j)
	endif
	enddo
	enddo
	write(*,*)'qmax=',qmax
	write(*,*)'volt=',volt,y(2,m1)*y(2,m1)*pai*(x(l1,2)-x(1,2))


c add radiation heat source

	call radi_c

c out_en: energy radiation into ambient
        out_en=0
c rr: radius, zb: bottom position, bt: top position
	call radi_cyl(rr_c_t,zb_c_t,zt_c_t,out_en)

	call radi_cyl(rr_c_l,zb_c_l,zt_c_l,out_en)

        write(*,*)'out_en from top and bottom holes',out_en

c for latent heat

      if(dtm.le.2)then
       relax_b=1
	else
       relax_b=2/dtm
	endif
	do j=2,j_rr-1
	i=i_rt
	ap(i,j)=ap(i,j)/relax_b
	con(i,j)=con(i,j)+(1.0-relax_b)*ap(i,j)*t(i,j)
	heat_flux=gr_lat/1000./3600.*3166.*10.58d6
	con(i,j)=con(i,j)+ak1(i,j)*heat_flux
	enddo

	do j=2,j_rr-1
	do i=2,l1
	if(x(i,j).ge.zb_lat.and.x(i,j).le.zt_lat)then
	ap(i,j)=ap(i,j)/relax_b
	con(i,j)=con(i,j)+(1.0-relax_b)*ap(i,j)*t(i,j)
	heat_flux=gr_lat/1000./3600.*3166.*10.58d6
	con(i,j)=con(i,j)-vol(i,j)*heat_flux/(zt_lat-zb_lat)
	endif
	enddo
	enddo

	bi=5
c We used sig=0.9 before for Sterling calculation, however the value may be
c  too large, choose sig=0.8 since 05.02.03
        sig=0.8
	ste=5.67051d-8


	tgas=800

        if(dtm.le.2)then
        relax_b=1
	else
        relax_b=2/dtm
	endif


	do j=2,m2
c lower
	do i=2,l2
	id=id_ele(i,j)
	if(id.ne.0.and.id.ne.4)then

c make sure that at every boundary element the boundary conditions are
c effective radiation and convective heat transfer
c	aim(i,j)=0
	id_e=id_ele(i,j)
	if(r(i,j).gt.rr_c_l)then
c foam
c	bi=5
	bi=0
	fview=1
	else
c	bi=5
	bi=0
	fview=0.011
	fview=0
	goto 100
	endif

	ifst=i
	radi_b=ak1(i,j)*(fview*ste*sig*(t(i,j)**4-293.0**4)
     1   +bi*(t(i,j)-tgas))
c	con(i,j)=con(i,j)-radi_b
	ap(i,j)=ap(i,j)/relax_b
	con(i,j)=con(i,j)+(1.0-relax_b)*ap(i,j)*t(i,j)

	con(i,j)=con(i,j)+ak1(i,j)*(fview*ste*sig*293.0**4+bi*tgas)
	ap(i,j)=ap(i,j)+ak1(i,j)*(fview*ste*sig*t(i,j)**3+bi)
	out_en=out_en+radi_b
c	write(*,*)con(i,j)
	goto 100
	endif
	enddo
c upper 
100	continue
	tgas=800

	do i=l2,2,-1
	id=id_ele(i,j)
	if(id.ne.0.and.id.ne.4)then

c	aip(i,j)=0

	id_e=id_ele(i,j)
	if(r(i,j).gt.rr_c_t)then
c	bi=5
	bi=0
	fview=1
	else
c	bi=5
	bi=0
	fview=0.052
	fview=0
	goto 200
	write(*,*)id,y(i,j),t(i,j)
	endif

	ilst=i
	radi_b=ak1(i+1,j)*(fview*ste*sig*(t(i,j)**4-293.0**4)
     1   +bi*(t(i,j)-tgas))
c	con(i,j)=con(i,j)-radi_b
	ap(i,j)=ap(i,j)/relax_b
	con(i,j)=con(i,j)+(1.0-relax_b)*ap(i,j)*t(i,j)

	con(i,j)=con(i,j)+ak1(i+1,j)*(fview*ste*sig*293.0**4+bi*tgas)
	ap(i,j)=ap(i,j)+ak1(i+1,j)*(fview*ste*sig*t(i,j)**3+bi)
	out_en=out_en+radi_b
	goto 200
	endif
	enddo
200	continue
	if(j.eq.m1/2)then
	write(*,*)'ifst=',ifst,ilst
	write(*,*)t(ifst,j),t(ilst,j)
	endif
	enddo

	write(*,*)'out energy through upper and lower walls=',out_en
c wall
	bi=10
	bi=0
	fview=1
c tamb: temperature of copper coil
	tamb=293
	tgas=800
c emissivity and radius
c sig1, r1: felt cylinder surface
c sig2, r2: coil
	sig1=0.8
	r1=0.125

        sig2=0.3
c        sig2=0.12
c sig2, 0.12, is underestimated.
c for radiation condition without reflection sig2=1

	r2=0.2007625

      if(idesign.eq.5)then
      r1=0.165
      r2=0.240
      sig1=0.8
      sig2=0.3

      endif

c tamb: ambient temperature

	sqt=0
	env_o=0
	jj0=0
	do i=2,l2

	sqt=sqt+ae1(i,m1)
	do j=m2,2,-1
	id=id_ele(i,j)
	if(id.ne.0.and.id.ne.4)then

c	ajp(i,j)=0

	env_o=env_o+ajm(i,j+jj0)
     1   *(t(i,j-1+jj0)-t(i,j+jj0))

c	write(*,*)ajm(i,j),i,j
c	radi_b=ae1(i,j+1)*(fview*ste*sig*(t(i,j)**4-tamb**4)
c formulation p.286,Thermal Radiation Heat Transfer, R. Siegel, J. Howell
        sig_eff0=1/(1/sig1+r1/r2*(1/sig2-1))


        sig_eff=sig_eff0
	xx0=x(i,j)
	del_xx=(r2-r1)*0.6
	del_xxo=(r2-r1)*0.3
	if(xx0.ge.0.00245-del_xxo.and.xx0.le.(0.00245+del_xx))then
	sig_eff=sig1+(xx0-(0.00245-del_xxo))/(del_xx+del_xxo)
     1     *(sig_eff0-sig1)
	endif
	if(xx0.ge.0.212-del_xx.and.xx0.le.0.212+del_xxo)then
	sig_eff=sig_eff0+(xx0-(0.212-del_xx))/(del_xx+del_xxo)
     1     *(sig1-sig_eff0)
	endif

	if(xx0.le.0.00245-del_xxo.or.xx0.ge.0.212+del_xxo)sig_eff=sig1

	radi_b=ae1(i,j+1)*
     1  (ste*(t(i,j)**4-tamb**4)*sig_eff
     1   +bi*(t(i,j)-tgas))
c	con(i,j)=con(i,j)-radi_b

	ap(i,j)=ap(i,j)/relax_b
	con(i,j)=con(i,j)+(1.0-relax_b)*ap(i,j)*t(i,j)

c	con(i,j)=con(i,j)+ae1(i,j+1)*(fview*ste*sig*tamb**4+bi*tgas)
	con(i,j)=con(i,j)+ae1(i,j+1)*
     1    (ste*tamb**4*sig_eff+bi*tgas)
	ap(i,j)=ap(i,j)+ae1(i,j+1)*
     1   (ste*t(i,j)**3*sig_eff+bi)
	out_en=out_en+radi_b

c	if(mod(ndt,5).eq.0.and.iter.eq.0)then
c	write(*,*)'out bound',x(i,j),t(i,j)
c	endif

	goto 530
	endif
	enddo
530	continue
	enddo
	write(*,*)'sqt=',sqt,y(2,m1)*(x(l1,2)-x(1,2))*2*pai
	write(*,*)'out energy through sidewall',env_o
	write(*,*)'out energy=',out_en

	do i=1,l1
	gam(i,1)=0.0
	ajm(i,2)=0
	f(i,1,nf)=f(i,2,nf)
	do j=2,m2
	id=id_bound(i,j)
	if(id.ne.0)then
	ajm(i,j)=0
	goto 111
	endif
	enddo
111	continue

	enddo
600	continue

	if(nf.ne.6)goto 700
c concentration
	call boun_b
	do i=1,l1
	gam(i,1)=0.0
	ajm(i,2)=0
	f(i,1,nf)=f(i,2,nf)
	write(*,*)'con',i,nf,f(i,1,nf),f(i,2,nf)
	do j=2,m2
	id=id_bound(i,j)
	if(id.ne.0)then
	ajm(i,j)=0
	goto 601
	endif
	enddo
601	continue
	enddo

	mm=j_rr
	do j=2,mm-1
c	growthsp(j)=0
c top
	i=i_rt
c J=fcl1+fpl1*f(i+1,j,nf)=-aim(i+1,j)*(f(i+1,j,nf)-f(i,j,nf))/ak1(i,j)
c flux kg/s/m2 for SiC2, ccen: mol/m3, growthsp: m/s, Mt=0.0521kg/mol
	fcl1=0.0521*growthsp(j)*ccen(i,j)
c	fcl1=0
	fpl1=0
	aimak1=aim(i+1,j)/ak1(i,j)
	aimak1=aip(i,j)/ak1(i,j)
c	write(*,*)'aim',aim(i+1,j)
	con(i,j)=con(i,j)-aip(i,j)*fcl1/(aimak1+fpl1)
	ap(i,j)=ap(i,j)-aip(i,j)*aimak1/(aimak1+fpl1)
	ap(i,j)=ap(i,j)+aip(i,j)
	aip(i,j)=0
	f(i+1,j,nf)=-fcl1/(aimak1+fpl1)+aimak1/(aimak1+fpl1)*f(i,j,nf)
	enddo

c	write(*,*)((con(i,j),i=1,l1),j=m1/2,m1/2)
700	continue
	return
	end

	subroutine bound_sig
	include 'sc.h'

      if(iter.ne.0) go to 63
      do 51 i=1,l1
      ruksi(i,1)=0.0
   51 ruksi(i,m1)=0.0
      do 52 j=1,m1
      rueta(1,j)=0.0
   52 rueta(l1,j)=0.0

      do 61 i=1,l1
      rueta(i,2)=0.0
   61 rueta(i,m1)=0.0
      do 62 j=1,m1
      ruksi(2,j)=0.0
   62 ruksi(l1,j)=0.0
   63 continue


      do 8 j=1,m1
      t(1,j)=2450
8	t(l1,j)=2300
	xin(1)=0
	yin(1)=2450
	xin(2)=(0.5+1.5)*0.0254
	yin(2)=2450
	xin(3)=(0.5+2.5)*0.0254
	yin(3)=2420
	xin(4)=(0.5+2.5+0.8)*0.0254
	yin(4)=2390
	xin(5)=(0.5+2.5+1.6)*0.0254
	yin(5)=2360
	xin(6)=(0.5+2.5+2.5)*0.0254
	yin(6)=2330
	xin(7)=(0.5+2.5+2.5+7./16.)*0.0254
	yin(7)=2300
	inum=7
	do i=1,l1
	call lineinp(x(i,m1),t(i,m1))
	enddo
	if(nf.eq.5)then
	do i=1,l1
	gam(i,1)=0.0
	enddo
	endif
      return
	end

c ------------------------------------------------------------------

      subroutine output
	include 'sc.h'
   91 format(///11x,'positions of the nodes in cartesian coordinates')
   92 format(///11x,'the nodes in cylinder coordinates (x,y)-(z,r)')
   93 format(11x,47(1h*),//)
   95 format(e12.5)
   96 format(6e12.5)
   88 format(6(1pe12.5,1x))
   87 format(1x,i2)
   86 format(1x,'title=','"zhanghui"')
   84 format(1x,'variables ="x","y","u","v","T(K)",
     1    "concentration","partial pressure","parti","p","pc","st"')
   83 format(1x,'zone t= "zone 1" , i=',i3,' , j=',i3,' , f= block')

        if(mod(ndt,iout_dat).eq.0)goto 80
      if(lstop) go to 80
      if(tpass.eq.0.and.mode.eq.0) print 91
      if(tpass.eq.0.and.mode.eq.1) print 92
      if(tpass.eq.0) return

      if(iter.eq.0) print 70,tpass
   70 format(////5x,' time passed : ',e10.4,' nondimensional time'/)

      if(iter.eq.0) print 71
   71 format(//' iter',5x,'smax',7x,'ssum',5x,'u(5,5)',4x,
     &'v(5,5)',4x,'t(5,5)',4x,'t(5,5)'/)
      print 72,iter,smax,ssum,u(5,5),v(5,5),t(5,5),t(5,5)
   72 format(i4,2x,10e11.4)
	write(*,*)'l1=',l1,'m1=',m1
c     if(iter.eq.0) print 71
c  71 format(//' iter',5x,'smax',7x,'ssum',5x,'res(1)',4x,
c    &'res(2)',4x,'res(5)',4x,'res(7)'/)
c     print 72,iter,smax,ssum,res(1),res(2),res(5),res(7)
c  72 format(i4,2x,10e11.4)
c     endif
      return

   80 print 81
   81 format(///)

        write(12,*)tpass,t(i_rt+1,1),t(i_rb-1,1),
     1   tmax,t_top,t_bot
      lpgrid=.false.
      if(mode.eq.0.and.lpgrid) print 91
      if(mode.eq.1.and.lpgrid) print 92
      if(lpgrid) print 93
      if(lpgrid) call print
      lpgrid=.false.
      call print


	call write_output
	call output_user
      call setup1

c   add this part to calculate Reynolds stress
      open(unit=10, file='plot.plt')
      rewind 10
      write(10,*)'title=','"zhanghui"'
      write(10,*)'variables ="x","y","u","v","T(K)"'
c       write(10,*)'variables ="x","y","u","v","T(K)",
c     1    "concentration","partial pressure","parti","p","pc","st"')
      write(10,83)l1,m1
      write(10,88)((y(i,j),i=1,l1),j=1,m1)
      write(10,88)((x(i,j),i=1,l1),j=1,m1)
      write(10,88)((v(i,j),i=1,l1),j=1,m1)
      write(10,88)((u(i,j),i=1,l1),j=1,m1)
      write(10,88)((t(i,j),i=1,l1),j=1,m1)
c      write(10,88)((ccen(i,j),i=1,l1),j=1,m1)
c      write(10,88)((ccen(i,j)*8.31*t(i,j),i=1,l1),j=1,m1)
c      write(10,88)((parti_p(i,j),i=1,l1),j=1,m1)
c      write(10,88)((p(i,j),i=1,l1),j=1,m1)
c      write(10,88)((pc(i,j),i=1,l1),j=1,m1)
c      write(10,88)((f(i,j,10),i=1,l1),j=1,m1)

	write(*,*)'con80',ccen(80,1)
	if(idesign.eq.3)then
	write(10,1001)
      write(10,*)7
c heater
	write(10,*)5
	write(10,*)0, -0.005
	write(10,*)0.050,-0.005
	write(10,*)0.050,0.195
      write(10,*)0.00, 0.195
      write(10,*)0.00,0.187
c heater
	write(10,*)4
	write(10,*)0,0.005
	write(10,*)0.035,0.005
	write(10,*)0.035,0.187
      write(10,*)0.00,0.187
c foam
        write(10,*)8
	write(10,*)0,-0.005
	write(10,*)0.0075,-0.005
        write(10,*)0.0075,-0.105
        write(10,*)0.145,-0.105
        write(10,*)0.145,0.255
        write(10,*)0.008,0.255
        write(10,*)0.008,0.205
        write(10,*)0,0.205
c crucible
	write(10,*)4
	write(10,*)0,0.029
	write(10,*)0.03,0.029
        write(10,*)0.03,0.173
        write(10,*)0,0.173
c crystal interface
	write(10,*)2
        write(10,*)0,0.165-crys_len
        write(10,*)0.03,0.165-crys_len
c argon and charge interface
	write(10,*)2
        write(10,*)0,0.163
        write(10,*)0.03,0.163
c argon and sic powder surface
	write(10,*)2
	write(10,*)0,0.113
	write(10,*)0.03,0.113

c for stick
c      write(10,*)6
c      write(10,*)0,0.184
c      write(10,*)0.015,0.184
c      write(10,*)0.015,0.193
c      write(10,*)0.008,0.199
c      write(10,*)0.008,0.224
c      write(10,*)0,0.224

	else if(idesign.eq.2)then	


	write(10,1001)
      write(10,*)7
c heater
	write(10,*)5
	write(10,*)0, -0.007
      write(10,*)0.075,-0.007
      write(10,*)0.075,0.207
	write(10,*)0.023, 0.207
	write(10,*)0.023,0.195
c heater
	write(10,*)4
	write(10,*)0,0.005
	write(10,*)0.060,0.005
	write(10,*)0.060,0.195
	write(10,*)0.023,0.195
c insulation
	write(10,*)12
	write(10,*)0,-0.007
	write(10,*)0.0075,-0.007
	write(10,*)0.0075,-0.055
	write(10,*)0.1225,-0.055
	write(10,*)0.1225,0
	write(10,*)0.125,0
      write(10,*)0.125,0.255+felt_incr
      write(10,*)0.105,0.255+felt_incr
      write(10,*)0.105,0.270+felt_incr
      write(10,*)0.008,0.270+felt_incr
      write(10,*)0.008,zb_c_t
      write(10,*)0,zb_c_t
c crucible
	write(10,*)4
	write(10,*)0,0.058
      write(10,*)0.0455,0.058
      write(10,*)0.0455,0.178
	write(10,*)0,0.178
c crystal interface
	write(10,*)2
	write(10,*)0,0.178-crys_len
	write(10,*)0.025,0.178-crys_len
c argon/charge interface
	write(10,*)2
	write(10,*)0,0.140
      write(10,*)0.0455,0.140

c for stick
	write(10,*)6
	write(10,*)0,0.193
	write(10,*)0.020,0.193
	write(10,*)0.020,0.203
	write(10,*)0.008,0.213
      write(10,*)0.008,zb_c_t
      write(10,*)0,zb_c_t

        else if(idesign.eq.4)then 


	write(10,1001)
      write(10,*)14
c bottom heater
      write(10,*)6
	write(10,*)0, -0.007
      write(10,*)0.075,-0.007
      write(10,*)0.075,0
      write(10,*)0.060,0
      write(10,*)0.060,0.005
      write(10,*)0.055,0.005

c side heater
      write(10,*)5
      write(10,*)0.060,0
      write(10,*)0.075,0
      write(10,*)0.075,z_incr+0.200
      write(10,*)0.060,z_incr+0.200
      write(10,*)0.060,0

c top heater
      write(10,*)7
      write(10,*)0.075,z_incr+0.200
      write(10,*)0.075,z_incr+0.207
      write(10,*)tos_d+0.002,z_incr+0.207
      write(10,*)tos_d+0.002,z_incr+0.195
      write(10,*)0.060,z_incr+0.195
      write(10,*)0.060,z_incr+0.200
      write(10,*)0.075,z_incr+0.200

c bottom foam
      write(10,*)7
	write(10,*)0,-0.007
	write(10,*)0.0075,-0.007
	write(10,*)0.0075,-0.055
	write(10,*)0.1225,-0.055
	write(10,*)0.1225,0
      write(10,*)0.075,0
      write(10,*)0.075,-0.007

c side felt
      write(10,*)5
      write(10,*)0.105,0
	write(10,*)0.125,0
        write(10,*)0.125,z_incr+0.255
        write(10,*)0.105,z_incr+0.255
      write(10,*)0.105,0
c side foam
      write(10,*)5
      write(10,*)0.085,0
      write(10,*)0.105,0
      write(10,*)0.105,z_incr+0.255
      write(10,*)0.085,z_incr+0.255
      write(10,*)0.085,0

c top foam
      write(10,*)7
      write(10,*)0.105,z_incr+0.255
        write(10,*)0.105,z_incr+0.270
        write(10,*)tos_a,z_incr+0.270
      write(10,*)tos_a,z_incr+0.240
      write(10,*)0.085,z_incr+0.240
      write(10,*)0.085,z_incr+0.255
      write(10,*)0.105,z_incr+0.255

c top felt
      write(10,*)3
        write(10,*)tos_a,z_incr+0.232
      write(10,*)0.085,z_incr+0.232
      write(10,*)0.085,0

c bottom crucible
      write(10,*)6
      write(10,*)0,0.040
      write(10,*)0.051,0.040
      write(10,*)0.051,z_incr+0.178
      write(10,*)0.0455,z_incr+0.178
        write(10,*)0.0455,0.044
      write(10,*)0,0.044
c top lid
      write(10,*)6
      write(10,*)0,z_incr+0.190
      write(10,*)0.051,z_incr+0.190
      write(10,*)0.051,z_incr+0.153
      write(10,*)0.055,z_incr+0.153
      write(10,*)0.055,z_incr+0.193
      write(10,*)0,z_incr+0.193
c top lid
      write(10,*)4
      write(10,*)0,z_incr+0.178
      write(10,*)0.051,z_incr+0.178
      write(10,*)0.051,z_incr+0.180
      write(10,*)0,z_incr+0.180

c crystal interface
	write(10,*)2
        write(10,*)0,z_incr+0.178-crys_len
        write(10,*)0.025,z_incr+0.178-crys_len
c argon/charge interface
	write(10,*)2
	write(10,*)0,0.140
        write(10,*)0.0455,0.140

c for stick
	write(10,*)6
        write(10,*)0,z_incr+0.193
      write(10,*)tos_d,z_incr+0.193
      write(10,*)tos_d,z_incr+0.193+tos_h2
      write(10,*)tos_a,z_incr+0.193+tos_h-tos_h1
      write(10,*)tos_a,z_incr+0.193+tos_h
      write(10,*)0,z_incr+0.193+tos_h

c pedestal
	write(10,1001)
      write(10,*)1
      
      write(10,*)3+5*8+4
      write(10,*)0,0.004
      write(10,*)0.055,0.004
      write(10,*)0.055,0.006
      do id=0,7
      write(10,*)0.010,0.005+0.001+id*0.004
      write(10,*)0.010,0.005+0.001+0.002+id*0.004
      write(10,*)0.055,0.005+0.001+0.002+id*0.004
      write(10,*)0.055,0.005+0.001+0.004+id*0.004
      write(10,*)0.01,0.005+0.001+0.004+id*0.004
      enddo
      write(10,*)0.01,0.040
      write(10,*)0.055,0.040
      write(10,*)0.055,0.042
      write(10,*)0.051,0.042

  

        else if(idesign.eq.5)then 


	write(10,1001)
      write(10,*)14
c bottom heater
      write(10,*)6
	write(10,*)0, -0.007
        write(10,*)0.095,-0.007
      write(10,*)0.095,0
      write(10,*)0.075,0
      write(10,*)0.075,0.005
      write(10,*)0.068,0.005

c side heater
      write(10,*)5
      write(10,*)0.075,0
      write(10,*)0.095,0
      write(10,*)0.095,z_incr+0.20
      write(10,*)0.075,z_incr+0.20
      write(10,*)0.075,0

c top heater
      write(10,*)7
      write(10,*)0.095,z_incr+0.20
        write(10,*)0.095,z_incr+0.207
      write(10,*)tos_d+0.002,z_incr+0.207
      write(10,*)tos_d+0.002,z_incr+0.195
      write(10,*)0.075,z_incr+0.195
      write(10,*)0.075,z_incr+0.20
      write(10,*)0.095,z_incr+0.20

c insulation bottom foam
      write(10,*)6
	write(10,*)0,-0.007
      write(10,*)0.0075,-0.007
      write(10,*)0.0075,-0.055
      write(10,*)0.1625,-0.055
      write(10,*)0.1625,0
      write(10,*)0.095,0

c side felt
      write(10,*)5
      write(10,*)0.145,0
      write(10,*)0.165,0
      write(10,*)0.165,z_incr+0.255
      write(10,*)0.145,z_incr+0.255
      write(10,*)0.145,0

c top foam
      write(10,*)7
      write(10,*)0.145,z_incr+0.255
      write(10,*)0.145,z_incr+0.270
        write(10,*)0.009,z_incr+0.270
        write(10,*)0.009,z_incr+0.240
        write(10,*)0.125,z_incr+0.240
        write(10,*)0.125,z_incr+0.255
        write(10,*)0.145,z_incr+0.255

c side foam
      write(10,*)5
      write(10,*)0.125,0
      write(10,*)0.145,0
      write(10,*)0.145,z_incr+0.255
      write(10,*)0.125,z_incr+0.255
      write(10,*)0.125,0

c inside felt
      write(10,*)7
        write(10,*)0.009,z_incr+0.231
        write(10,*)0.009,z_incr+0.232
        write(10,*)0.125,z_incr+0.232
        write(10,*)0.125,0
        write(10,*)0.095,0
        write(10,*)0.095,z_incr+0.207
        write(10,*)tos_d,z_incr+0.207

c crucible
      write(10,*)6
      write(10,*)0,0.044
        write(10,*)0.060,0.044
        write(10,*)0.060,z_incr+0.178
      write(10,*)0.064,z_incr+0.178
      write(10,*)0.064,0.040
      write(10,*)0,0.040
c top lid
      write(10,*)6
      write(10,*)0.0,z_incr+0.190
      write(10,*)0.064,z_incr+0.190
      write(10,*)0.064,z_incr+0.153
      write(10,*)0.068,z_incr+0.153
      write(10,*)0.068,z_incr+0.193
      write(10,*)0,z_incr+0.193
c lid
      write(10,*)4
      write(10,*)0,z_incr+0.178
      write(10,*)0.064,z_incr+0.178
      write(10,*)0.064,z_incr+0.180
      write(10,*)0,z_incr+0.180

c crystal interface
      write(10,*)3
        write(10,*)0,z_incr+0.178-crys_len
      write(10,*)0.05,z_incr+0.178-crys_len
      write(10,*)0.05,z_incr+0.190
c argon/charge interface
	write(10,*)2
	write(10,*)0,0.140
        write(10,*)0.060,0.140

c for stick
	write(10,*)6
        write(10,*)0.0,z_incr+0.193
      write(10,*)tos_d,z_incr+0.193
      write(10,*)tos_d,z_incr+0.193+tos_h2
      write(10,*)tos_a,z_incr+0.193+tos_h-tos_h1
      write(10,*)tos_a,z_incr+0.193+tos_h
      write(10,*)0.0,z_incr+0.193+tos_h

	write(10,1001)
      write(10,*)1

c pedestal
      
      write(10,*)3+5*8+4
      write(10,*)0,0.004
      write(10,*)0.068,0.004
      write(10,*)0.068,0.006
      do id=0,7
      write(10,*)0.010,0.005+0.001+id*0.004
      write(10,*)0.010,0.005+0.001+0.002+id*0.004
      write(10,*)0.068,0.005+0.001+0.002+id*0.004
      write(10,*)0.068,0.005+0.001+0.004+id*0.004
      write(10,*)0.01,0.005+0.001+0.004+id*0.004
      enddo
      write(10,*)0.01,0.040
      write(10,*)0.068,0.040
      write(10,*)0.068,0.042
      write(10,*)0.064,0.042

        else if(idesign.eq.6)then 


	write(10,1001)
      write(10,*)9
c bottom heater
      write(10,*)6
	write(10,*)0, -0.007
      write(10,*)0.075,-0.007
      write(10,*)0.075,0
      write(10,*)0.060,0
      write(10,*)0.060,0.005
      write(10,*)0.0,0.005

c side heater
      write(10,*)5
      write(10,*)0.060,0
      write(10,*)0.075,0
      write(10,*)0.075,z_incr+0.200
      write(10,*)0.060,z_incr+0.200
      write(10,*)0.060,0

c top heater
      write(10,*)6
      write(10,*)0.0,z_incr+0.195
      write(10,*)0.060,z_incr+0.195
      write(10,*)0.060,z_incr+0.200
      write(10,*)0.075,z_incr+0.200
      write(10,*)0.075,z_incr+0.207
      write(10,*)0.0,z_incr+0.207

c bottom foam
      write(10,*)7
	write(10,*)0,-0.007
	write(10,*)0.0075,-0.007
	write(10,*)0.0075,-0.055
	write(10,*)0.1225,-0.055
	write(10,*)0.1225,0
      write(10,*)0.075,0
      write(10,*)0.075,-0.007

c side felt
      write(10,*)5
      write(10,*)0.105,0
      write(10,*)0.1225,0
        write(10,*)0.1225,z_incr+0.255
        write(10,*)0.105,z_incr+0.255
      write(10,*)0.105,0
c side foam
      write(10,*)5
      write(10,*)0.085,0
      write(10,*)0.105,0
      write(10,*)0.105,z_incr+0.255
      write(10,*)0.085,z_incr+0.255
      write(10,*)0.085,0

c top foam
      write(10,*)7
      write(10,*)0.105,z_incr+0.255
        write(10,*)0.105,z_incr+0.270
        write(10,*)0.0225,z_incr+0.270
      write(10,*)0.0225,z_incr+0.240
      write(10,*)0.085,z_incr+0.240
      write(10,*)0.085,z_incr+0.255
      write(10,*)0.105,z_incr+0.255

c top felt
      write(10,*)5
        write(10,*)0.015,z_incr+0.240
      write(10,*)0.085,z_incr+0.240
      write(10,*)0.085,z_incr+0.207
      write(10,*)0.015,z_incr+0.207
      write(10,*)0.015,z_incr+0.240
c argon/charge interface
	write(10,*)2
      write(10,*)0,0.1315
      write(10,*)0.046,0.1315


	write(10,1001)
      write(10,*)5

c bottom crucible
      write(10,*)6
      write(10,*)0,0.0365
      write(10,*)0.051,0.040
      write(10,*)0.051,z_incr+0.1765
      write(10,*)0.046,z_incr+0.1765
        write(10,*)0.046,0.0415
      write(10,*)0,0.0415
c top crucible cover
      write(10,*)6
      write(10,*)0,z_incr+0.183
      write(10,*)0.051,z_incr+0.183
      write(10,*)0.051,z_incr+0.1365
      write(10,*)0.055,z_incr+0.1365
      write(10,*)0.055,z_incr+0.186
      write(10,*)0,z_incr+0.186
c seed lid
      write(10,*)4
      write(10,*)0,z_incr+0.1765
      write(10,*)0.051,z_incr+0.1765
      write(10,*)0.051,z_incr+0.178
      write(10,*)0,z_incr+0.178

c crystal interface
	write(10,*)2
      write(10,*)0,z_incr+0.1765-crys_len
      write(10,*)0.025,z_incr+0.1765-crys_len



c for stick
	write(10,*)6
      write(10,*)0,z_incr+0.186
      write(10,*)tos_d,z_incr+0.186
      write(10,*)tos_d,z_incr+0.186+tos_h2
      write(10,*)tos_a,z_incr+0.186+tos_h-tos_h1
      write(10,*)tos_a,z_incr+0.186+tos_h
      write(10,*)0,z_incr+0.186+tos_h

c pedestal
	write(10,1001)
      write(10,*)1
      
      write(10,*)3+5*8+4
      write(10,*)0.0,0.004
      write(10,*)0.055,0.004
      write(10,*)0.055,0.006
      do id=0,7
      write(10,*)0.010,0.005+0.001+id*0.004
      write(10,*)0.010,0.005+0.001+0.002+id*0.004
      write(10,*)0.055,0.005+0.001+0.002+id*0.004
      write(10,*)0.055,0.005+0.001+0.004+id*0.004
      write(10,*)0.01,0.005+0.001+0.004+id*0.004
      enddo
      write(10,*)0.01,0.040
      write(10,*)0.055,0.040
      write(10,*)0.055,0.042
      write(10,*)0.051,0.042

	endif




c for coils
c coil_xx0, coil_yy0, coil_dxx from radi_ini

      if(idesign.eq.1.or.idesign.eq.2.or.idesign.eq.3.or.
     1   idesign.eq.4.or.idesign.eq.6)then
	write(10,1001)
      write(10,*)10


	do 700 id=1,5
      xx0=coil_xx0+coil_dxx*(id-1)
      yy0=coil_yy0
c        xx_r=0.01905
c        yy_r=0.0047625
	write(10,*)5
	write(10,*)yy0-yy_r,xx0-xx_r
	write(10,*)yy0+yy_r,xx0-xx_r
	write(10,*)yy0+yy_r,xx0+xx_r
	write(10,*)yy0-yy_r,xx0+xx_r
	write(10,*)yy0-yy_r,xx0-xx_r
c        xx_r=0.01905-width_coil
c        yy_r=0.0047625-width_coil
        xx_r1=xx_r-width_coil
        yy_r1=yy_r-width_coil
	write(10,*)5
        write(10,*)yy0-yy_r1,xx0-xx_r1
        write(10,*)yy0+yy_r1,xx0-xx_r1
        write(10,*)yy0+yy_r1,xx0+xx_r1
        write(10,*)yy0-yy_r1,xx0+xx_r1
        write(10,*)yy0-yy_r1,xx0-xx_r1
700	continue
1001	format(1x,'GEOMETRY X=0 Y=0 T=LINE CS=GRID C=BLACK')

      endif

      if(idesign.eq.5)then
      do 800 id=1,10
      xx0=coil_xx0+coil_dxx*(id-1)
      yy0=coil_yy0
	write(10,111)yy0,xx0
      write(10,*)0.25*0.0254
	write(10,111)yy0,xx0
      write(10,*)0.25*0.0254-width_coil
800   continue
111	format('GEOMETRY X=',f7.3,' Y=',f7.3,
     1     ' T=CIRCLE, CS=GRID, C=BLACK')

      endif

	close(10)
c	write(*,*)((t(i,j),i=1,l1),j=m1/2,m1/2)
      return
	end



      subroutine write_output
      include 'sc.h'

        goto 100
      open(unit=9, file='output', form='unformatted')
      write(9) tpass,l1,m1
      write(9)((x(i,j),i=1,l1),j=1,m1)
      write(9)((y(i,j),i=1,l1),j=1,m1)
      write(9)((u(i,j),i=1,l1),j=1,m1)
      write(9)((v(i,j),i=1,l1),j=1,m1)
      write(9)((p(i,j)/(1+cdif),i=1,l1),j=1,m1)
      write(9)((t(i,j),i=1,l1),j=1,m1)
      write(9)((ccen(i,j),i=1,l1),j=1,m1)
      write(9)((ruksi(i,j),i=1,l1),j=1,m1)
      write(9)((rueta(i,j),i=1,l1),j=1,m1)
c stream
      write(9)((f(i,j,10),i=1,l1),j=1,m1)
      write(9)((qth(i,j),i=1,l1),j=1,m1)
c for magnetic field
      write(9)((fc(i,j),i=1,l1),j=1,m1)
	close(9)
        goto 200

100   open(unit=9, file='output')
      write(9,*) tpass,l1,m1
      write(9,196)((x(i,j),i=1,l1),j=1,m1)
      write(9,196)((y(i,j),i=1,l1),j=1,m1)
      write(9,196)((u(i,j),i=1,l1),j=1,m1)
      write(9,196)((v(i,j),i=1,l1),j=1,m1)
      write(9,196)((p(i,j)/(1+cdif),i=1,l1),j=1,m1)
      write(9,196)((t(i,j),i=1,l1),j=1,m1)
      write(9,196)((ccen(i,j),i=1,l1),j=1,m1)
      write(9,196)((ruksi(i,j),i=1,l1),j=1,m1)
      write(9,196)((rueta(i,j),i=1,l1),j=1,m1)
c stream
      write(9,196)((f(i,j,10),i=1,l1),j=1,m1)
      write(9,196)((qth(i,j),i=1,l1),j=1,m1)
c for magnetic field

c (fc,f(1,1,8)
      write(9,196)((f(i,j,8),i=1,l1),j=1,m1)
      write(9,196)((f(i,j,9),i=1,l1),j=1,m1)
	close(9)

196     format(6e20.12)
200     continue
      return
      end



	subroutine boun_b
	include 'sc.h'
c left boundary
	do 1005 i=2,l2
	jb=m2
c search for solid part
1001	do j=jb,2,-1
	if(id_bound(i,j).ne.0)then
	jb=j
	goto 1002
	endif
	enddo
	goto 1005
1002	do j=jb,2,-1
	if(id_bound(i,j).eq.0)then
	jb=j
	goto 1003
	endif
	enddo
	goto 1005
1003	do j=jb,1,-1
	if(id_bound(i,j).eq.0)then
	je=j
	else
	goto 1004
	endif
	enddo
1004	write(*,*)i,j,k,'boun_b'
	if(nf.eq.1)then
	do j1=jb,je,-1
	u(i,j1)=0
	enddo
	else if(nf.eq.2)then
	do j1=jb,je,-1
	v(i,j1)=0
	enddo
	else if(nf.eq.3)then
	do j1=jb,je,-1
	p(i,j1)=0
	enddo
	if(jb+1.le.m2)then
	ajm(i,jb+1)=0
	endif
	p(i,jb)=p(i,jb+1)
	else if(nf.eq.4)then
	do j1=jb,je,-1
	pc(i,j1)=0
	enddo
	if(jb+1.le.m2)then
	ajm(i,jb+1)=0
	endif
	pc(i,jb)=pc(i,jb+1)
	endif

	if(nf.eq.6)then
	do j1=jb,je,-1
	ccen(i,j1)=0
	enddo

	if(jb+1.le.m2)then
	ajm(i,jb+1)=0
	endif
	ccen(i,jb)=ccen(i,jb+1)


	endif

	goto 1001
1005	continue

c right boundary
	do 2005 i=2,l2
	jb=2
2001	do j=jb,m2
	if(id_bound(i,j).ne.0)then
	jb=j
	goto 2002
	endif
	enddo
	goto 2005
2002	do j=jb,m2
	if(id_bound(i,j).eq.0)then
	jb=j
	goto 2003
	endif
	enddo
	goto 2005
2003	do j=jb,m1
	if(id_bound(i,j).eq.0)then
	je=j
	else
	goto 2004
	endif
	enddo
2004	continue
c	write(*,*)i,j,'boun_b'
	if(nf.eq.1)then
	do j1=jb,je
	u(i,j1)=0
	enddo
	else if(nf.eq.2)then
	do j1=jb,je
	v(i,j1)=0
	enddo
	else if(nf.eq.3)then
c	write(*,*)'p in boun_b',p(i,jb-1),p(i,jb)
	if(jb-1.ge.2)then
c	ajp(i,jb-1)=0
	endif
	do j1=jb,je
c	p(i,j1)=0
	enddo
c	p(i,jb)=p(i,jb-1)
	p(i,jb)=p(i,jb-1)+(p(i,jb-1)-p(i,jb-2))
     1     /0.5/(heta(i,jb-1)+heta(i,jb-2))
     1    *0.25*(heta(i,jb)+heta(i,jb-1))
	else if(nf.eq.4)then
	do j1=jb,je
c	pc(i,j1)=0
	enddo
	if(jb-1.ge.2)then
c	ajp(i,jb-1)=0
	endif
c	pc(i,jb)=pc(i,jb-1)
	pc(i,jb)=pc(i,jb-1)+(pc(i,jb-1)-pc(i,jb-2))
     1     /0.5/(heta(i,jb-1)+heta(i,jb-2))
     1    *0.25*(heta(i,jb)+heta(i,jb-1))
	endif


	if(nf.eq.6)then
	if(jb-1.ge.2)then
	ajp(i,jb-1)=0
	endif

	do j1=jb,je
	ccen(i,j1)=0
	enddo
	ccen(i,jb)=ccen(i,jb-1)

	endif

	goto 2001
2005	continue
	write(*,*)'crystal speed'

c upper boundary
	do 3005 j=2,m2
	ib=2
3001	do i=ib,l2
	if(id_bound(i,j).ne.0)then
	ib=i
	goto 3002
	endif
	enddo
	goto 3005
3002	do i=ib,l2
	if(id_bound(i,j).eq.0)then
	ib=i
	goto 3003
	endif
	enddo
	goto 3005
3003	do i=ib,l1
	if(id_bound(i,j).eq.0)then
	ie=i
	else
	goto 3004
	endif
	enddo
3004	continue
	if(nf.eq.1)then
	u(ib,j)=0
	else if(nf.eq.2)then
	v(ib,j)=0
	else if(nf.eq.3)then
	aip(ib-1,j)=0
	p(ib,j)=p(ib-1,j)
	else if(nf.eq.4)then
	if(ib-1.ge.2)then
	aip(ib-1,j)=0	
	pc(ib,j)=pc(ib-1,j)
	endif
	else if(nf.eq.6)then
	if(ib-1.ge.2)then
c	aip(ib-1,j)=0
c	ccen(ib,j)=ccen(ib-1,j)
	endif
	endif

	goto 3001
3005	continue

c lower boundary
        do 4005 j=2,m2
        ib=l2
4001	do i=ib,2,-1
	if(id_bound(i,j).ne.0)then
	ib=i
	goto 4002
	endif
	enddo
	goto 4005
4002	do i=ib,2,-1
	if(id_bound(i,j).eq.0)then
	ib=i
	goto 4003
	endif
	enddo
	goto 4005
4003	do i=ib,1,-1
	if(id_bound(i,j).eq.0)then
	ie=i
	else
	goto 4004
	endif
	enddo
4004	continue
	if(nf.eq.1)then
	do i1=ib,ie,-1
	u(i1,j)=0
	enddo
	else if(nf.eq.2)then
	do i1=ib,ie,-1
	v(i1,j)=0
	enddo
	else if(nf.eq.3)then
	do i1=ib,ie,-1
c	p(i1,j)=0
	enddo
	if(ib+1.le.l2)then
c	aim(ib+1,j)=0
	endif
c	p(ib,j)=p(ib+1,j)
	p(ib,j)=p(ib+1,j)+(p(ib+1,j)-p(ib+2,j))
     1   /0.5/(hksi(ib+1,j)+hksi(ib+2,j))
     1   *0.25*(hksi(ib,j)+hksi(ib+1,j))
	else if(nf.eq.4)then
	do i1=ib,ie,-1
c	pc(i1,j)=0
	enddo
	if(ib+1.le.l2)then
c	aim(ib+1,j)=0
	endif
c	pc(ib,j)=pc(ib+1,j)
	pc(ib,j)=pc(ib+1,j)+(pc(ib+1,j)-pc(ib+2,j))
     1   /0.5/(hksi(ib+1,j)+hksi(ib+2,j))
     1   *0.25*(hksi(ib,j)+hksi(ib+1,j))
	else if(nf.eq.6)then
	if(ib+1.le.l2)then
	aim(ib+1,j)=0
	endif
	ccen(ib,j)=ccen(ib+1,j)
	endif

	goto 4001
4005	continue
	return
      end


      subroutine gamsor
	include 'sc.h'
      do 101 i=1,l1
      do 101 j=1,m1


      gam(i,j)=viscosity(i,j)
c	if ((nf.eq.1).or.(nf.eq.2)) then
c      if(i.ge.lst3) gam(i,j)=1.d10
c	endif
101	continue   

      if(nf.ne.1) go to 201

      return


  201 if(nf.ne.2) go to 301
      return


301	if(nf.ne.5) goto 601
c for temperature
      do 310 i=1,l1
      do 310 j=1,m1

c	gam(i,j)=gam_cal(i,j)
	gam(i,j)=th_conduc(i,j)
310	continue
	return

601	if(nf.ne.6)goto 701
c concentration
	do i=1,l1
	do j=1,m1
c molecular weight of SiC2, Mt=0.0521 kg/mol
	gam(i,j)=0.009*0.0521

c	gam(i,j)=0
	enddo
	enddo
	return

701	continue
	if(nf.eq.10)then
	do i=1,l1
	do j=1,m1
	gam(i,j)=1
	enddo
	enddo
	endif
      return
      end


	subroutine radi_ini
c this subroutine is called after reading input.dat
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt

	include "sc.h"

	z_incr=0

        coil_xx0=coil_xx0+coil_move

c rr_c_t, zb_c_t, zt_c_t for top hole
c rr_c_l, zb_c_l, zt_c_l for lower hole
c these are used in subroutine bound


      if(idesign.eq.1.or.idesign.eq.2)then
      tos_a=0.008
      tos_d=0.025
	else if(idesign.eq.3)then
      tos_a=0.008
      tos_d=0.025
        else if(idesign.eq.4)then
      tos_a=0.008
      tos_d=0.025
      tos_h1=0.018
      tos_h2=0.010
      tos_h=0.038
        else if(idesign.eq.5)then
      tos_a=0.008
      tos_d=0.050
      tos_h1=0.018
      tos_h2=0.010
      tos_h=0.038
        else if(idesign.eq.6)then
      tos_a=0.0125
      tos_d=0.0225
      tos_h1=0.018
      tos_h2=0.010
      tos_h=0.034

      endif


c ar: area of each cell
c only used in radi_c and cal_ex1



	do i=l1,1,-1
	if(x(i,2).lt.zb0)then
	i_rb=i+1
	if(id_ele(i,2).eq.4)then
	write(*,*)id_ele(i,2),x(i,2),r(i,2),i
	write(*,*)(r(l1/2,j),j=1,m1)
	stop 'err in radi_c zb0'
	endif
	goto 20
	endif
	enddo
20	do i=1,l1
	if(x(i,2).gt.zt0)then
	i_rt=i-1
	if(id_ele(i,2).eq.4)then
	write(*,*)id_ele(i,2),x(i,2),r(i,2)
	stop 'err in radi_c zt0'
	endif
	goto 40
	endif
	enddo
40	do j=1,m1
	if(y(i_rt,j).gt.rc)then
	j_rr=j-1
	if(id_ele(i_rt,j).eq.4)then
	stop 'err in radi_c rc'
	endif
	goto 60
	endif
	enddo
60	return
	end

c for top and bottom holes

	subroutine radi_cyl(rr_c,zb_c,zt_c,out_en)

	include 'sc.h'
c determine rr_c: radius of cylindrical enclosure, 
c           zb_c: coordinate of bottom, 
c           zt_c: coordinate of top
c ib_c: axial coordinate of the bottom point in the enclosure
c it_c: axial coordinate of the top point in the enclosure.
c jr_c: radial coordinate of the rightest point
c incre: 1 means in normal position, -1: upside down

	if(zb_c.le.zt_c)then
	do i=l1,1,-1
	if(x(i,2).lt.zb_c)then
	ib_c=i+1

	goto 20
	endif
	enddo
20	do i=1,l1
	if(x(i,2).gt.zt_c)then
	it_c=i-1

	goto 40
	endif
	enddo
40	nn=it_c-ib_c+1
c nn: number of points in vertical direction
	incre=1
	else
	do i=l1,1,-1
	if(x(i,2).lt.zb_c)then
	ib_c=i
	goto 50
	endif
	enddo
50	do i=1,l1
	if(x(i,2).gt.zt_c)then
	it_c=i
	goto 55
	endif
	enddo
55	nn=abs(it_c-ib_c)+1
	incre=-1
	endif

	do j=1,m1
	if(y(it_c,j).gt.rr_c)then
	jr_c=j-1
	goto 60
	endif
	enddo
60	mm=jr_c

c bottom line
	do j=1,mm
	rb(j)=r(ib_c,j)
	zb(j)=0
	tb(j)=t(ib_c-incre,j)
	enddo

c side wall
	do i=1,nn
	rw(i)=rr_c
	zw(i)=x((i-1)*incre+ib_c,jr_c)-zb_c
	tw(i)=t((i-1)*incre+ib_c,jr_c+1)
	enddo

c top line
	do j=1,mm
	rs(j)=r(it_c,j)
	zs(j)=zt_c-zb_c
c set as environmental temperature
        tsurf(j)=300
c        tsurf(j)=t(it_c+incre,j)
	enddo
        eb=0.8
        ew=0.8
c es environment, 
        es=1.
c	write(*,*)'radi_cy',mm,nn,ib_c,it_c,jr_c,zb_c,zt_c,rr_c
c	write(*,*)tb(1)
c	pause

c    ra0,za0,rb0,zb0,rc0,zc0,rd0,zd0

c eb: base, ew: side wall, es: top surface
c for top and bottom holes
	call radi_g(mm,nn,mm,rb,zb,tb,rw,zw,tw,rs,zs,tsurf,0d0,0d0,
     1   rr_c,0d0,rr_c,zt_c-zb_c,0d0,zt_c-zb_c,eb,ew,es,qb,qw,qsurf,
     1   pb,pw,ps,arb,arw,ars)

c	call radi(mm,nn,tb,rc,zt0-zb0,zw,tw,rs,zs,tsurf,eb,ew,es,
c     1  qb,qw,qsurf,iprint)
c	write(*,*)(qb(i),i=1,mm)
c	write(*,*)(qw(i),i=1,nn)
c	write(*,*)(qsurf(i),i=1,mm)

        if(dtm.le.2)then
	relax_b=1
	else
        relax_b=2/dtm
	endif

	do i=2,nn-1
c sidewall from i_rb+1 to i_rt-1
	i1=(i-1)*incre+ib_c
	j1=jr_c+1
	ss=qw(i)*arw(i)
	if(i.eq.2)ss=ss+qw(1)*arw(1)
	if(i.eq.nn-1)ss=ss+qw(nn)*arw(nn)
	ap(i1,j1)=ap(i1,j1)/relax_b
	pp=pw(i)*arw(i)
	if(i.eq.2)pp=pp+pw(1)*arw(1)
	if(i.eq.nn-1)pp=pp+pw(nn)*arw(nn)

	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)
	con(i1,j1)=con(i1,j1)+ss
	ap(i1,j1)=ap(i1,j1)-pp


	enddo

c        goto 30
c open boundary at the top

	do j=2,mm-1
c top
	i1=it_c+incre
	j1=j
	ss=qsurf(j)*ars(j)
	if(j.eq.2)ss=ss+qsurf(1)*ars(1)
	if(j.eq.mm-1)ss=ss+qsurf(mm)*ars(mm)
	pp=ps(j)*ars(j)
	if(j.eq.2)pp=pp+ps(1)*ars(1)
	if(j.eq.mm-1)pp=pp+ps(mm)*ars(mm)

c	ap(i1,j1)=ap(i1,j1)/relax_b
c	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)
c	con(i1,j1)=con(i1,j1)+ss
c	ap(i1,j1)=ap(i1,j1)-pp

	out_en=out_en+ss+pp*t(i1,j1)


	enddo

30	continue

	do j=2,mm-1
c bottom

	i1=ib_c-incre
	j1=j
	ss=qb(j)*arb(j)
	if(j.eq.2)ss=ss+qb(1)*arb(1)
	if(j.eq.mm-1)ss=ss+qb(mm)*arb(mm)
	pp=pb(j)*arb(j)
	if(j.eq.2)pp=pp+pb(1)*arb(1)
	if(j.eq.mm-1)pp=pp+pb(mm)*arb(mm)

	ap(i1,j1)=ap(i1,j1)/relax_b
	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)

	con(i1,j1)=con(i1,j1)+ss
	ap(i1,j1)=ap(i1,j1)-pp


	enddo


	return
	end

	subroutine radi_c
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt
	include "sc.h"
c ar: area of each cell
c determin rc, zb0, zt0
	call radi_ini
40	nn=i_rt-i_rb+1

60	mm=j_rr

	do j=1,mm
	rb(j)=r(i_rb,j)
	zb(j)=0
	tb(j)=t(i_rb-1,j)
	enddo

	do i=1,nn
	rw(i)=rc
	zw(i)=x(i+i_rb-1,j_rr)-zb0
	tw(i)=t(i+i_rb-1,j_rr+1)
	enddo

	do j=1,mm
	rs(j)=r(i_rt,j)
	zs(j)=zt0-zb0
        tsurf(j)=t(i_rt+1,j)

	enddo
	eb=0.8
	ew=0.8
        es=0.8
c	write(*,*)mm,nn,i_rb,i_rt,j_rr,zb0,zt0
c	write(*,*)tb(1)
c	pause

c for growth chamber
	call radi_g(mm,nn,mm,rb,zb,tb,rw,zw,tw,rs,zs,tsurf,0d0,0d0,
     1   rc,0d0,rc,zt0-zb0,0d0,zt0-zb0,eb,ew,es,qb,qw,qsurf,
     1   pb,pw,ps,arb,arw,ars)

c testing part
c	call radi_g(1,nn,mm,rb,zb,tb,rw,zw,tw,rs,zs,tsurf,0d0,0d0,
c     1   rc,0d0,rc,zt0-zb0,0d0,zt0-zb0,eb,ew,es,qb,qw,qsurf,
c     1   pb,pw,ps,arb,arw,ars)
c	write(*,*)(qb(i),i=1,mm)
c	write(*,*)(qw(i),i=1,nn)
c	write(*,*)(qsurf(i),i=1,mm)
c	call radi(mm,nn,tb,rc,zt0-zb0,zw,tw,rs,zs,tsurf,eb,ew,es,
c     1  qb,qw,qsurf,iprint)
c	write(*,*)(qb(i),i=1,mm)
c	write(*,*)(qw(i),i=1,nn)
c	write(*,*)(qsurf(i),i=1,mm)

        if(dtm.le.2)then
	relax_b=1
	else
        relax_b=2/dtm
	endif

	do i=2,nn-1
c sidewall from i_rb+1,i_rt-1

	i1=i+i_rb-1

	j1=j_rr+1
	ss=qw(i)*arw(i)
	if(i.eq.2)ss=ss+qw(1)*arw(1)
	if(i.eq.nn-1)ss=ss+qw(nn)*arw(nn)
	pp=pw(i)*arw(i)
	if(i.eq.2)pp=pp+pw(1)*arw(1)
	if(i.eq.nn-1)pp=pp+pw(nn)*arw(nn)

	ap(i1,j1)=ap(i1,j1)/relax_b
	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)

	con(i1,j1)=con(i1,j1)+ss
	ap(i1,j1)=ap(i1,j1)-pp
c	write(*,*)con(i1,j1),qw(i),ae1(i1,j1)
c	pause
	enddo

	do j=2,mm-1
c top

	i1=i_rt+1
	j1=j
	ss=qsurf(j)*ars(j)
	if(j.eq.2)ss=ss+qsurf(1)*ars(1)
	if(j.eq.mm-1)ss=ss+qsurf(mm)*ars(mm)
	pp=ps(j)*ars(j)
	if(j.eq.2)pp=pp+ps(1)*ars(1)
	if(j.eq.mm-1)pp=pp+ps(mm)*ars(mm)

	ap(i1,j1)=ap(i1,j1)/relax_b
	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)

	con(i1,j1)=con(i1,j1)+ss
	ap(i1,j1)=ap(i1,j1)-pp
c	write(*,*)con(i1,j1),t(i1,j1)
	enddo

	do j=2,mm-1
c bottom
	i1=i_rb-1
	j1=j
	ss=qb(j)*arb(j)
	if(j.eq.2)ss=ss+qb(1)*arb(1)
	if(j.eq.mm-1)ss=ss+qb(mm)*arb(mm)
	pp=pb(j)*arb(j)
	if(j.eq.2)pp=pp+pb(1)*arb(1)
	if(j.eq.mm-1)pp=pp+pb(mm)*arb(mm)

	ap(i1,j1)=ap(i1,j1)/relax_b
	con(i1,j1)=con(i1,j1)+(1.0-relax_b)*ap(i1,j1)*t(i1,j1)

	con(i1,j1)=con(i1,j1)+ss
	ap(i1,j1)=ap(i1,j1)-pp
c	write(*,*)ak1(i1,j1),ak1(i1+1,j1)
c	pause
	enddo

c	if(mod(ndt,5).eq.0.and.iter.eq.0)then
c	write(*,*)'center t(i_rb-1->i_rt+1)',(t(i,2),i=i_rb-1,i_rt+1)
c	write(*,*)'wall t(i_rb-1->i_rt+1)',
c     1     (t(i,j_rr+1),i=i_rb-1,i_rt+1)
c	write(*,*)'upper t(1->j_rr+1)',(t(i_rt+1,j),j=1,j_rr+1)
c	write(*,*)'lower t(1->j_rr+1)',(t(i_rb-1,j),j=1,j_rr+1)
c	endif
	return
	end


	subroutine radi_g(m1,m2,m3,rb,zb,tb,rw,zw,tw,rs,zs,ts,
     1  ra0,za0,rb0,zb0,rc0,zc0,rd0,zd0,eb,ew,es,qb,qw,qs,
     1  pb,pw,ps,arb,arw,ars)

	parameter (nnn=500)
	implicit real*8(a-h,o-z)
	dimension r(nnn),z(nnn),t(nnn),em(nnn),ar(nnn),
     1    f(201,201),a(201,201),b(201,201),c(201,201),
     1    q(nnn),p(nnn),
     1    rb(nnn),zb(nnn),tb(nnn),qb(nnn),pb(nnn),arb(nnn),
     1    rw(nnn),zw(nnn),tw(nnn),qw(nnn),pw(nnn),arw(nnn),
     1    rs(nnn),zs(nnn),ts(nnn),qs(nnn),ps(nnn),ars(nnn)

	do i=1,m1
	r(i)=rb(i)
	z(i)=zb(i)
	t(i)=tb(i)
	em(i)=eb
	enddo
	do i=1,m2
	r(i+m1)=rw(i)
	z(i+m1)=zw(i)
	t(i+m1)=tw(i)
	em(i+m1)=ew
	enddo
	do i=1,m3
	r(i+m1+m2)=rs(m3-i+1)
	z(i+m1+m2)=zs(m3-i+1)
	t(i+m1+m2)=ts(m3-i+1)
	em(i+m1+m2)=es
	enddo
	m=m1+m2+m3
	if(m.ge.500)stop 'err radi_g'

c	write(*,*)m1,m2,m3
c	write(*,*)(r(i),i=1,m)
c	write(*,*)(z(i),i=1,m)
c	write(*,*)(t(i),i=1,m)
c	write(*,*)(em(i),i=1,m)

c calculate vier factor f
	call view_2ring(m1,m2,m3,m,ra0,za0,rb0,zb0,rc0,zc0,rd0,zd0,
     1   r,z,ar,f)
	ste=5.67d-8
	do i=1,m
	do j=1,m
	a(i,j)=-f(i,j)*(1-em(j))/em(j)

	b(i,j)=f(i,j)
	if(i.eq.j)then
	a(i,j)=a(i,j)+1/em(i)
	b(i,j)=b(i,j)-1
	endif
c	c(i,j)=a(i,j)
	enddo
	enddo

	call minv(a,m)
c	do i=1,m
c	do j=1,m
c	sum=0
c	do k=1,m
c	sum=sum+a(i,k)*c(k,j)
c	enddo
c	write(*,*)'ac=',i,j,sum
c	enddo
c	enddo

	do i=1,m
	sum=0
	do j=1,m
	c(i,j)=0
	do k=1,m
	c(i,j)=c(i,j)+a(i,k)*b(k,j)
	enddo
	sum=sum+c(i,j)
	enddo
c	write(*,*)i,sum
	enddo
	sum=0
	do i=1,m
	p(i)=0
	q(i)=0

c If we set p(i), then it may cause numerical instability

	do j=1,m
c	if(j.eq.i)then
c	p(i)=c(i,j)*ste*t(j)**3
c	else
	q(i)=q(i)+c(i,j)*ste*t(j)**4
c	endif
	enddo
	sum=sum+q(i)*ar(i)

	enddo


	do i=1,m1
	qb(i)=q(i)
	pb(i)=p(i)
	arb(i)=ar(i)

	enddo
	do i=1,m2
	qw(i)=q(i+m1)
	pw(i)=p(i+m1)
	arw(i)=ar(i+m1)

	enddo
	do i=1,m3
	qs(m3-i+1)=q(i+m1+m2)
	ps(m3-i+1)=p(i+m1+m2)
	ars(m3-i+1)=ar(i+m1+m2)
c	write(*,*)'qs=',qs(m3-i+1),ts(m3-i+1)
	enddo
	return
	end


	subroutine radi(m,n,tb,rc,l,z,tw,rs,zs,tsurf,
     1   eb,ew,es,qb,qw,qsurf,iprint)
      implicit real*8(a-h,o-z)
      external disk, minv


c    A fortran program for radiative exchange between surfaces of
c    a cylindrical crystal growth furnace.

       dimension fwb(100),fbw(100),fww(100,100)
     1,fws(100,100),fsw(100,100),fss(100,100),fbs(100)
     2,fsb(100),aw(100),as(100),
     1  z(100),tw(100),qw(100),
     1  rs(100),zs(100),tsurf(100),qsurf(100),
     1  r(100),zm(100),ts(100),qs(100),
     1  a(201,201)
     4,b(201,201),c(201,201),tfwb(100),tfbw(100),tfww(100,100)
     5,tfws(100,100),tfsw(100,100),tfss(100,100),tfbs(100)
     6,tfsb(100),ASURF(100),awmo(100)

c      common/areaview/m,n,tb,eb,es,ew,unit
      real*8 l
      LOGICAL fsupdate
      integer unit
c      save aw,as,Tfbw,Tfww,Tfsw,Tfbs,Tfws,Tfss

c       namelist / input / n,m,rc,z,rs,zs,tw,tsurf,tb,eb,es,ew
c      open(unit=2,file='cyl.dat',status='old')
c      read(2,input)

      pi=3.141592654
      sigma=5.67d-8
c  sigma=1.0
      nt=1+n+m

	write(*,*)m,n,tb
c	pause
      fsupdate=.true.
      unit=6

      do i=1,m
        r(i)=rs(m+1-i)
        zm(i)=zs(m+1-i)
        ts(i)=tsurf(m+1-i)
      end do

c      l=z(n)+(z(n)-z(n-1))/2.
c      l=z(n)

      if(fsupdate)then




c     Calculate view factors between the base and the side wall
	if(iprint.eq.1)then
      WRITE(UNIT,1001)
	endif
 1001 FORMAT(' ',3X,'I',9X,'Fbw',14X,'Fwb')

      do i=1,n
        if(i.eq.1)then
          z1=(z(2)+z(1))/2.
          aw(i)=2.*pi*rc*z1
          call disk(rc,rc,z1,ff)
          fbw(1)=1-ff
          fwb(1)=rc/2./z1*fbw(1)
        else
          z1=(z(i)+z(i-1))/2.
           if(i.eq.n)then
            z2=l
           else
          z2=(z(i+1)+z(i))/2.
           end if
          aw(i)=2.*pi*rc*(z2-z1)
          call disk(rc,rc,z1,ff1)
          call disk(rc,rc,z2,ff2)
          fbw(i)=ff1-ff2
          fwb(i)=rc/2./(z2-z1)*fbw(i)
        end if
C        write(6,*)i,' fbw=',fbw(i),' fwb=',fwb(i)
	if(iprint.eq.1)then
         WRITE(UNIT,1002)I,FBW(I),FWB(I)
	endif
 1002    FORMAT(' ',I4,2(2X,F14.7))

      end do

c     Calculate view factors between the base and the surface

	if(iprint.eq.1)then
      WRITE(UNIT,1003)
	endif
 1003 FORMAT(' ',3X,'I',9X,'Fbs',14X,'Fsb')

      do i=1,m-1
        if(i.eq.1)then
          r1=rc
          z1=L
        else
          r1=(r(i)+r(i-1))/2.
          z1=(zm(i)+zm(i-1))/2.
        end if
       z2=(zm(i)+zm(i+1))/2.
       r2=(r(i)+r(i+1))/2.
       as(i)=pi*(r1+r2)*dsqrt((r1-r2)**2+(z1-z2)**2)
c       write(6,*)'rc=',rc,' r1=',r1,' r2=',r2,' z1=',z1,' z2=',z2
       call disk(rc,r1,z1,ff1)
       call disk(rc,r2,z2,ff2)
       fbs(i)=ff1-ff2
c       write(6,*)'ff1= ',ff1,'  ff2=',ff2
       fsb(i)=rc**2/(r1+r2)/dsqrt((r1-r2)**2+(z1-z2)**2)*fbs(i)
C        write(6,*)i,' fbs=',fbs(i),' fsb=',fsb(i)
	if(iprint.eq.1)then
         WRITE(UNIT,1002)I,FBS(I),FSB(I)
	endif

      end do

      z1=zm(m)
      r1=(r(m)+r(m-1))/2.
      as(m)=pi*r1**2 
cout      as(m)=pi*r1*dsqrt(r1**2+(zm(m)-zm(m-1))**2)
      call disk(rc,r1,z1,ff1)
      fbs(m)=ff1
      fsb(m)=pi*rc**2/as(m)*fbs(m) 
cout      fsb(m)=rc**2/r1**2*fbs(m) 
C        write(6,*)m,' fbs=',fbs(m),' fsb=',fsb(m)
	if(iprint.eq.1)then
	WRITE(UNIT,1002)M,FBS(M),FSB(M) 
	endif

c      check view factors from the base

      sum=0.
      do i=1,n
        sum=sum+fbw(i)
      end do
      do i=1,m
        sum=sum+fbs(i)
      end do
	if(iprint.eq.1)then
      write(UNIT,*)'SUM FOR BASE =',sum
	endif

c     Calculate view factors between walls
	if(iprint.eq.1)then
      WRITE(UNIT,1004)
	endif
 1004 FORMAT(' ',3X,'I',3X,'J',9X,'Fww(i,j)',10X,'Fww(j,i)')

      do i=1,n
       do j=1,n
       if(i.gt.j)go to 1
         if(i.eq.j)then
          if(i.eq.1.or.i.eq.n)then
           if(i.eq.1)zz=(z(2)+z(1))/2.
           if(i.eq.n)zz=l-(z(n-1)+z(n))/2.
          else
           z1=(z(i)+z(i-1))/2.
           z2=(z(i+1)+z(i))/2.
           zz=z2-z1
          end if
          call disk(rc,rc,zz,f1)
          fww(i,i)=1.-rc/zz*(1.-f1)
         else
             if(i.eq.1)then
          z1=0.
             else
          z1=(z(i)+z(i-1))/2.
             end if
          z2=(z(i)+z(i+1))/2.
          z3=(z(j)+z(j-1))/2.
           if(j.eq.n)then
            z4=l
           else
          z4=(z(j)+z(j+1))/2.
           end if
          z2z3=z3-z2
          z2z4=z4-z2
          z1z3=z3-z1
          z1z4=z4-z1
           if(abs(z2z3).le.0.000001)then
             f23=1.
           else
          call disk(rc,rc,z2z3,f23)
           end if
          call disk(rc,rc,z2z4,f24)
          call disk(rc,rc,z1z3,f13)
          call disk(rc,rc,z1z4,f14)
          fww(j,i)=rc/2./(z4-z3)*(f23-f24-f13+f14)
          fww(i,j)=(z4-z3)/(z2-z1)*fww(j,i)
         end if
C        write(6,*)i,j,' fww=',fww(i,j),' fww=',fww(j,i)
	if(iprint.eq.1)then
         WRITE(UNIT,1005)I,J,FWW(I,J),FWW(J,I)
	endif
 1005    FORMAT(' ',2(1X,I4),2(2X,F14.7))

   1    continue
        end do
       end do

c      Calculate view factors between wall and surface
	if(iprint.eq.1)then
      WRITE(UNIT,1006)
	endif
 1006 FORMAT(' ',3X,'I',3X,'J',9X,'Fws(i,j)',10X,'Fsw(j,i)')

      do i=1,n
       do j=1,m
        if(i.eq.1)then
         z1=0.
        else
         z1=(z(i)+z(i-1))/2.
        end if
        if(i.eq.n)then
          z2=l
        else
        z2=(z(i)+z(i+1))/2.
        end if
        if(j.eq.1)then
         z3=l
         r3=rc
        else
         z3=(zm(j)+zm(j-1))/2.
         r3=(r(j)+r(j-1))/2.
        end if
          if(j.lt.m)then
        z4=(zm(j)+zm(j+1))/2.
        r4=(r(j)+r(j+1))/2.
          else
        z4=z3
        r4=0.
          end if
        z2z3=z3-z2
        z2z4=z4-z2
        z1z3=z3-z1
        z1z4=z4-z1
          if(abs(z2z3).le.0.00001)then
        f23=(r3/rc)**2
          else
        call disk(rc,r3,z2z3,f23)
          end if
          if(j.eq.m)then
             f24=0.
          else
          if(dabs(z2z4).le.0.0000001)then
             f24=(r4/rc)**2
          else
        call disk(rc,r4,z2z4,f24)
          end if
          end if
        call disk(rc,r3,z1z3,f13)
          if(j.eq.m)then
              f14=0.
          else
        call disk(rc,r4,z1z4,f14)
          end if
        fws(i,j)=rc/2./(z2-z1)*(f23-f24-f13+f14)
        fsw(j,i)=2.*rc*(z2-z1)/(r3+r4)/dsqrt((r3-r4)**2+(z3-z4)**2)
     1 *fws(i,j)
C        write(6,*)i,j,' fws=',fws(i,j),' fsw=',fsw(j,i)
	if(iprint.eq.1)then
         WRITE(UNIT,1005)I,J,FWS(I,J),FSW(J,I)     
	endif

       end do
      end do

c     Check view factors from the wall

      do i=1,n
        sum=fwb(i)
       do j=1,n
        sum=sum+fww(i,j)
       end do

       do j=1,m
        sum=sum+fws(i,j)
       end do
	if(iprint.eq.1)then
        write(UNIT,*) i, ' SUM FOR WALL =',sum
	endif
      end do

c     Calculate view factors between two surface rings
	if(iprint.eq.1)then
      WRITE(UNIT,1007)
	endif
 1007 FORMAT(' ',3X,'I',3X,'J',9X,'Fss(i,j)',10X,'Fss(j,i)')

      do i=1,m
       do j=1,m
        if(i.gt.j)go to 2
          if(i.eq.j)then
           if(i.eq.1)then
            r1=rc
            r2=(r(1)+r(2))/2.
            zz=(zm(2)+zm(1))/2.-l
            go to 3
           end if
           if(i.eq.m)then
            fss(i,i)=0.0
cout            zz=(zm(m)-zm(m-1))/2.
cout            r1=(r(m)+r(m-1))/2.
cout            fss(i,i)=1-r1/dsqrt(r1**2+zz**2)
            go to 2
           end if
           zz=(zm(i+1)-zm(i-1))/2.
           r1=(r(i-1)+r(i))/2.
           r2=(r(i+1)+r(i))/2.
    3      if(dabs(zz).le.0.000001)then
           fss(i,i)=0.0
           else
           call disk(r1,r2,zz,f12)
           fss(i,i)=1.-1./(r1+r2)/dsqrt(zz**2+(r2-r1)**2)*
     1 (r1**2+r2**2-2.*r1**2*f12)
           end if
          go to 2
        end if 
        if(i.eq.1)then
          z1=l
          r1=rc
        else
          z1=(zm(i)+zm(i-1))/2.
          r1=(r(i)+r(i-1))/2.
        end if
          z3=(zm(j)+zm(j-1))/2.
          r3=(r(j)+r(j-1))/2.
          z2=(zm(i)+zm(i+1))/2.
          r2=(r(i)+r(i+1))/2.
           if(j.eq.m)then
             z4=z3
             r4=0.
           else
          z4=(zm(j)+zm(j+1))/2.
          r4=(r(j)+r(j+1))/2.
            end if
        z2z3=z3-z2
        z2z4=z4-z2
        z1z3=z3-z1
        z1z4=z4-z1  
         if(dabs(z2z3).le.0.000001)then
        f23=1.
         else
        call disk(r2,r3,z2z3,f23)
         end if
         if(dabs(z2z4).le.0.0000001)then
        f24=1.
         else
        call disk(r2,r4,z2z4,f24)
         end if
         if(dabs(z1z3).le.0.0000001)then
        f13=1.
         else
        call disk(r1,r3,z1z3,f13)
         end if
         if(dabs(z1z4).le.0.0000001)then
        f14=1.
         else
        call disk(r1,r4,z1z4,f14)
         end if
          aa=f13-f14
          bb=f23-f24
C         write(6,*)aa,bb

c        if(i.eq.m-1.and.j.eq.m)then
c          fss(j,i)=1.-r1**2/r3**2*f13
c          fss(i,j)=r3**2/(r1+r3)/dsqrt((r1-r3)**2+(z1-z3)**2)
c     1   *fss(j,i)
c          go to 2
c        end if


        fss(j,i)=1./(r3+r4)/dsqrt((r3-r4)**2+(z3-z4)**2)*
     1 (r2**2*(f23-f24)-r1**2*(f13-f14))
        fss(i,j)=1./(r1+r2)/dsqrt((r1-r2)**2+(z1-z2)**2)*
     1 (r2**2*(f23-f24)-r1**2*(f13-f14))
   2    continue
C        write(6,*)i,j,' fss=',fss(i,j),' fss=',fss(j,i)
	if(iprint.eq.1)then
         WRITE(UNIT,1005)I,J,FSS(I,J),FSS(J,I)     
	endif

       end do
      end do

c     Check view factors from surface

      do i=1,m
       sum=fsb(i)
       do j=1,n
       sum=sum+fsw(i,j)
       end do
       do j=1,m
       sum=sum+fss(i,j) 
       end do
	if(iprint.eq.1)then
       write(UNIT,*) i, ' SUM FOR SURFACE =',sum
	endif
      end do

c     Calculate total exchange factors

      do i=1,nt
        do j=1,nt
        if(i.eq.1)then
          if(j.eq.1)then
            a(i,j)=0.
            b(i,j)=0.
          end if
          if(j.gt.1.and.j.le.n+1)then
            a(i,j)=-fbw(j-1)*(1.-ew)
            b(i,j)=fbw(j-1)*ew
          end if
          if(j.gt.n+1)then
            a(i,j)=-fbs(j-n-1)*(1.-es) 
            b(i,j)=fbs(j-n-1)*es
          end if
        end if
        if(i.gt.1.and.i.le.n+1)then
          if(j.eq.1)then
            a(i,j)=-fwb(i-1)*(1.-eb)
            b(i,j)=fwb(i-1)*eb
          end if
          if(j.gt.1.and.j.le.n+1)then
            a(i,j)=-fww(i-1,j-1)*(1.-ew)
            b(i,j)=fww(i-1,j-1)*ew
          end if
          if(j.gt.n+1)then
            a(i,j)=-fws(i-1,j-n-1)*(1.-es)
            b(i,j)=fws(i-1,j-n-1)*es
          end if
        end if
        if(i.gt.n+1)then
          if(j.eq.1)then
            a(i,j)=-fsb(i-n-1)*(1.-eb)
            b(i,j)=fsb(i-n-1)*eb
          end if
          if(j.gt.1.and.j.le.n+1)then
            a(i,j)=-fsw(i-n-1,j-1)*(1.-ew)
            b(i,j)=fsw(i-n-1,j-1)*ew
          end if
          if(j.gt.n+1)then
            a(i,j)=-fss(i-n-1,j-n-1)*(1.-es)
            b(i,j)=fss(i-n-1,j-n-1)*es
          end if
        end if
        if(i.eq.j)a(i,j)=a(i,j)+1.
          end do
        end do

c       Evaluate inverse of matrix a

        call minv(a,nt)

        do i=1,nt
          do j=1,nt
          c(i,j)=0.
          do k=1,nt
          c(i,j)=c(i,j)+a(i,k)*b(k,j)
          end do
          end do
        end do
c
c       Convert c(i,j) to total exchange factors
c
        do i=1,nt
          do j=1,nt
          if(i.eq.1)then
           if(j.eq.1)tfbb=c(i,j)
           if(j.gt.1.and.j.le.n+1)tfbw(j-1)=c(i,j)
           if(j.gt.n+1)tfbs(j-n-1)=c(i,j)
          end if

          if(i.gt.1.and.i.le.n+1)then
           if(j.eq.1)tfwb(i-1)=c(i,j)
           if(j.gt.1.and.j.le.n+1)tfww(i-1,j-1)=c(i,j)
           if(j.gt.n+1)tfws(i-1,j-n-1)=c(i,j)
          end if

          if(i.gt.n+1)then
           if(j.eq.1)tfsb(i-n-1)=c(i,j)
           if(j.gt.1.and.j.le.n+1)tfsw(i-n-1,j-1)=c(i,j)
           if(j.gt.n+1)tfss(i-n-1,j-n-1)=c(i,j)
          end if
          end do
        end do

c      check view factors from the base

      sum=tfbb
      do i=1,n
        sum=sum+tfbw(i)
      end do
      do i=1,m
        sum=sum+tfbs(i)
      end do
	if(iprint.eq.1)then
      write(UNIT,*)'SUM OF TEF FROM BASE=',sum
	endif

c     Check view factors from the wall
c
      do i=1,n
        sum=tfwb(i)
       do j=1,n
        sum=sum+tfww(i,j)
       end do

       do j=1,m
        sum=sum+tfws(i,j)
       end do
	if(iprint.eq.1)then
        write(UNIT,*) i, ' SUM OF TEF FROM WALL=',sum
	endif
      end do

c     Check view factors from surface

      do i=1,m
       sum=tfsb(i)
       do j=1,n
       sum=sum+tfsw(i,j)
       end do
       do j=1,m
       sum=sum+tfss(i,j) 
       end do
	if(iprint.eq.1)then
        write(UNIT,*) i, ' SUM OF TEF FROM SURFACE=',sum
	endif
      end do

	end if


c     Calculate heat fluxes

      ab=pi*rc**2
      do i=1,n
       qw(i)=ab*eb*sigma*tb**4*tfbw(i)-aw(i)*ew*sigma*tw(i)**4
         do j=1,n
       qw(i)=qw(i)+aw(j)*ew*sigma*tw(j)**4*tfww(j,i)
         end do
         do j=1,m
       qw(i)=qw(i)+as(j)*es*sigma*ts(j)**4*tfsw(j,i)
         end do
C        write(6,*)i,' qw=', qw(i)
       end do

       do i=1,m
        qs(i)=ab*eb*sigma*tb**4*tfbs(i)-as(i)*es*sigma*ts(i)**4
          do j=1,n
        qs(i)=qs(i)+aw(j)*ew*sigma*tw(j)**4*tfws(j,i)
          end do
          do j=1,M
        qs(i)=qs(i)+as(j)*es*sigma*ts(j)**4*tfss(j,i)
          end do
C        write(6,*)i,' qs=', qs(i),' aw=',aw(i),' as=',as(i)
        end do

c ab: area, eb: emissivity, tb : temperature

	qb=ab*eb*sigma*tb**4*tfbb-ab*eb*sigma*tb**4
	do j=1,n
	qb=qb+aw(j)*ew*sigma*tw(j)**4*tfwb(j)
	enddo
	do j=1,m
	qb=qb+as(j)*es*sigma*ts(j)**4*tfsb(j)
	enddo

c	goto 100
c test
      do i=1,n
       qw(i)=aw(i)*ew*sigma*tb**4*tfwb(i)-aw(i)*ew*sigma*tw(i)**4
         do j=1,n
       qw(i)=qw(i)+aw(i)*ew*sigma*tw(j)**4*tfww(i,j)
         end do
         do j=1,m
       qw(i)=qw(i)+aw(i)*ew*sigma*ts(j)**4*tfws(i,j)
         end do
C        write(6,*)i,' qw=', qw(i)
       end do

       do i=1,m
        qs(i)=as(i)*es*sigma*tb**4*tfsb(i)-as(i)*es*sigma*ts(i)**4
          do j=1,n
        qs(i)=qs(i)+as(i)*es*sigma*tw(j)**4*tfsw(i,j)
          end do
          do j=1,M
        qs(i)=qs(i)+as(i)*es*sigma*ts(j)**4*tfss(i,j)
          end do
C        write(6,*)i,' qs=', qs(i),' aw=',aw(i),' as=',as(i)
        end do

c ab: area, eb: emissivity, tb : temperature

	qb=ab*eb*sigma*tb**4*tfbb-ab*eb*sigma*tb**4
	do j=1,n
	qb=qb+ab*eb*sigma*tw(j)**4*tfbw(j)
	enddo
	do j=1,m
	qb=qb+ab*eb*sigma*ts(j)**4*tfbs(j)
	enddo
100	continue

	sum=0

      do i=1,m
        Asurf(i)=As(m+1-i)
	sum=sum+qs(i)
        qsurf(i)=qs(m+1-i)/asURF(i)
	if(iprint.eq.1)then
        write(UNIT,*)i,' qsurf=', qsurf(i),' asURF=',asURF(i)
	endif
      end do

      do i=1,n
	sum=sum+qw(i)
        awmo(i)=aw(i)
        qw(i)=qw(i)/awmo(i)
	if(iprint.eq.1)then
        write(UNIT,*)i,' qw=', qw(i),' aw=',aw(i)
	endif
      end do
	sum=sum+qb

	qb=qb/ab
	if(iprint.eq.1)then
	write(unit,*)'qb=',qb,'sum(qb+qw+qs)=',sum
	endif

      return
      end


	subroutine view_2ring(m1,m2,m3,m,ra,za,rb,zb,rc,zc,rd,zd,
     1   r1m,z1m,ar,f)
	implicit real*8(a-h,o-z)
c f12 configuration factor from ring 1 to ring 2
c f21 configuration factor from ring 2 to ring 1
	dimension r1m(m),z1m(m),ar(m),f(201,201)
c	if(iprint.eq.1)then
c      WRITE(UNIT,1007)
c	endif
c 1007 FORMAT(' ',3X,'I',3X,'J',9X,'Fss(i,j)',10X,'Fss(j,i)')
	pai=3.14159265
      do i=1,m
       do j=1,m
	if(i.gt.j)goto 20
        if(i.eq.1)then
          r1=ra
          z1=za
	else if(i.eq.m1+1)then
	r1=rb
	z1=zb
	else if(i.eq.m1+m2+1)then
	r1=rc
	z1=zc
        else
          r1=(r1m(i)+r1m(i-1))/2.
          z1=(z1m(i)+z1m(i-1))/2.
        end if
	if(i.eq.m1)then
	r2=rb
	z2=zb
	else if(i.eq.m1+m2)then
	r2=rc
	z2=zc
	else if(i.eq.m1+m2+m3)then
	r3=rd
	z3=zd
	else
          r2=(r1m(i)+r1m(i+1))/2.
          z2=(z1m(i)+z1m(i+1))/2.
	endif

	if(j.eq.1)then
	r3=ra
	z3=za
	else if(j.eq.m1+1)then
	r3=rb
	z3=zb
	else if(j.eq.m1+m2+1)then
	r3=rc
	z3=zc
	else
          r3=(r1m(j)+r1m(j-1))/2.
          z3=(z1m(j)+z1m(j-1))/2.
	endif
	if(j.eq.m1)then
	r4=rb
	z4=zb
	else if(j.eq.m1+m2)then
	r4=rc
	z4=zc
	else if(j.eq.m1+m2+m3)then
	r4=rd
	z4=zd
	else
          r4=(r1m(j)+r1m(j+1))/2.
          z4=(z1m(j)+z1m(j+1))/2.
	end if
	if(i.eq.1)then
	ar(j)=pai*(r3+r4)*dsqrt((r3-r4)**2+(z3-z4)**2)
	endif
c	write(*,*)r1,z1,r2,z2,r3,z3,r4,z4
        z2z3=(z3-z2)
        z2z4=(z4-z2)
        z1z3=(z3-z1)
        z1z4=(z4-z1)
         if(dabs(z2z3).le.0.000001)then
        f23=1.
         else
        call disk1(r2,r3,z2z3,f23)
         end if
         if(dabs(z2z4).le.0.0000001)then
        f24=1.
         else
        call disk1(r2,r4,z2z4,f24)
         end if
         if(dabs(z1z3).le.0.0000001)then
        f13=1.
         else
        call disk1(r1,r3,z1z3,f13)
         end if
         if(dabs(z1z4).le.0.0000001)then
        f14=1.
         else
        call disk1(r1,r4,z1z4,f14)
         end if
          aa=f13-f14
          bb=f23-f24
C         write(6,*)aa,bb
c        if(i.eq.m-1.and.j.eq.m)then
c          fss(j,i)=1.-r1**2/r3**2*f13
c          fss(i,j)=r3**2/(r1+r3)/dsqrt((r1-r3)**2+(z1-z3)**2)
c     1   *fss(j,i)
c          go to 2
c        end if

	if(abs(r3+r4).le.0.0001.or.
     1     abs(r1+r2).le.0.0001.or.
     1     abs(r3-r4)+abs(z3-z4).le.0.0001.or.
     1     abs(r1-r2)+abs(z1-z2).le.0.0001)then
	f(j,i)=0
	f(i,j)=0
	else
        f(j,i)=1./ar(j)*pai*(r2**2*(f23-f24)-r1**2*(f13-f14))
	if(f(j,i).lt.0)f(j,i)=f(j,i)+1

c        f(i,j)=1./(r1+r2)/dsqrt((r1-r2)**2+(z1-z2)**2)*
c     1 (r2**2*(f23-f24)-r1**2*(f13-f14))
        f(i,j)=ar(j)/ar(i)*f(j,i)
	endif
c	write(*,*)f13,f23,f14,f24
c        write(*,*)i,j,' f(i,j)=',f(i,j),' f(j,i)=',f(j,i)
20	continue
c	if(iprint.eq.1)then
c         WRITE(UNIT,1005)I,J,FSS(I,J),FSS(J,I)     
c	endif

       end do
      end do
	do i=1,m
	sum=0
	do j=1,m
	sum=sum+f(i,j)
	enddo
c	write(*,*)sum
	enddo
	return
	end

      subroutine disk(r1,r2,h,f1)
      implicit real*8(a-h,o-z)

c     Subroutine for calculating view factors between two disks
      rr1=r1/h
      rr2=r2/h
      x=1.+(1.+rr2**2)/rr1**2
      f1=1./2.*(x-dsqrt(x**2-4.*(rr2/rr1)**2))
      return 
      end

      subroutine disk1(r1,r2,h,f1)
      implicit real*8(a-h,o-z)

c     Subroutine for calculating view factors between two disks
	if(r1.le.0.00001.or.r2.le.0.00001)then
	f1=0
	return
	endif
      rr1=r1/h
      rr2=r2/h
      x=1.+(1.+rr2**2)/rr1**2
      f1=1./2.*(x-dsqrt(x**2-4.*(rr2/rr1)**2))
      return 
      end


      subroutine minv(c,j3)

      implicit real*8(a-h,o-z)

      dimension c(201,201),j(1500)
      do 125 i=1,j3
  125 j(i+20)=i
      do 144 i=1,j3
      c0=0.
      j1=i
      do 135 k=i,j3
      if((dabs(c0)-dabs(c(i,k))).ge.0.)go to 135
  126 j1=k
      c0=c(i,k)
  135 continue
  127 if(i.eq.j1)go to 138
  128 k=j(j1+20)
      j(j1+20)=j(i+20)
      j(i+20)=k
      do 137 k=1,j3
      s0=c(k,i)
      c(k,i)=c(k,j1)
  137 c(k,j1)=s0
  138 c(i,i)=1
      do 139 j1=1,j3
  139 c(i,j1)=c(i,j1)/c0
      do 142 j1=1,j3
      if(i.eq.j1)go to 142
  129 c0=c(j1,i)
      if(c0.eq.0)go to 142
  130 c(j1,i)=0.0
      do 141 k=1,j3
  141 c(j1,k)=c(j1,k)-c0*c(i,k)
  142 continue
  144 continue
      do 143 i=1,j3
      if(j(i+20).eq.i)go to 143
  131 j1=i
  132 j1=j1+1
      if(j(j1+20).eq.i)go to 133
  136 if(j3.gt.j1)go to 132
  133 j(j1+20)=j(i+20)
      do 163 k=1,j3
      c0=c(i,k)
      c(i,k)=c(j1,k)
  163 c(j1,k)=c0
      j(i+20)=i
  143 continue
      return
      end


      subroutine change_grid()
      include 'sc.h'
	
        iun=1
	if(iun.eq.1)then
	
      open(unit=18, file='input')
      read(18,*) tinit,l1o,m1o
	else
      open(unit=18, file='input',form='unformatted')
      read(18) tinit,l1o,m1o
	endif
      write(*,*)'in change grid'
      write(*,*)'l1,m1,l1o,m1o=',l1,m1,l1o,m1o

	if(iun.eq.1)then
      read(18,*)((xo(i,j),i=1,l1o),j=1,m1o)
      read(18,*)((yo(i,j),i=1,l1o),j=1,m1o)
	else
      read(18)((xo(i,j),i=1,l1o),j=1,m1o)
      read(18)((yo(i,j),i=1,l1o),j=1,m1o)
	endif

	if(iun.eq.1)then
      read(18,*)((uch(i,j,1),i=1,l1o),j=1,m1o)
      read(18,*)((uch(i,j,2),i=1,l1o),j=1,m1o)
      read(18,*)((uch(i,j,3),i=1,l1o),j=1,m1o)
      read(18,*)((uch(i,j,5),i=1,l1o),j=1,m1o)
	else
      read(18)((uch(i,j,1),i=1,l1o),j=1,m1o)
      read(18)((uch(i,j,2),i=1,l1o),j=1,m1o)
      read(18)((uch(i,j,3),i=1,l1o),j=1,m1o)
      read(18)((uch(i,j,5),i=1,l1o),j=1,m1o)
	endif

	close(unit=18)
	do i=1,l1o
	do j=1,m1o
        uch(i,j,3)=uch(i,j,3)*(1+cdif)
	enddo
	enddo

      do loop=1,nfch
	write(*,*)'loop=',loop
	fmax=1.d-6
	do i=1,l1o
	do j=1,m1o
      fijk=dabs(uch(i,j,loop))
	fmax=dmax1(fmax,fijk)
	enddo
	enddo
	write(*,*)fmax
	enddo

      call interp(l1o,m1o)

      do loop=1,nfch
	fmax=1.d-6
	do i=1,l1
	do j=1,m1
      fijk=dabs(f(i,j,loop))
	fmax=dmax1(fmax,fijk)

	enddo
	enddo
	write(*,*)fmax
	enddo

        tpass=tinit
      call write_output

c get ruksi,rueta,ruzeta
	call initl
c get stream function
	call output

        write(*,*)'finish changing grid'
      stop
	end


      subroutine interp0(xxx,yyy,ni,nj,l1,m1,
     1   x,y,u,uuu)
	implicit real*8 (a-h,o-z)
        parameter (mss=400)
      dimension x(ni,nj),y(ni,nj)
      dimension u(ni,nj)
      common /interp_c/ssd,is,js
      dimension ssd(mss),is(mss),js(mss)
c      ss0=1d-8
      ss0=(1./ni)**2
10	continue
	iis=0
      do 100 i1=1,l1
      do 100 j1=1,m1
      xx1=x(i1,j1)
      yy1=y(i1,j1)
      ss=(xx1-xxx)**2+(yy1-yyy)**2
	if(ss.le.1d-8)then
      uuu=u(i1,j1)
	goto 500
	endif

	if(ss.le.ss0)then
	iis=iis+1
	ssd(iis)=ss

      if(iis.gt.mss)then
	write(*,*)'iis in interp0=',iis
      stop 'err in interp0'
      endif

	is(iis)=i1
	js(iis)=j1
	endif
100	continue

      if(iis.lt.1)then
      ss0=ss0*2
	goto 10
	endif

	uuu1=0
	sss=0
	do ii1=1,iis
	i1=is(ii1)
	j1=js(ii1)
	ss=1/sqrt(ssd(ii1))
      uuu1=uuu1+u(i1,j1)*ss
	sss=sss+ss
	enddo
	uuu=uuu1/sss
500	return
	end



      subroutine interp(l1o,m1o)
        parameter (mss=400)
      include 'sc.h'
      common /interp_c/ssd,is,js
      dimension ssd(mss),is(mss),js(mss)
      if(l1o.gt.ni.or.m1o.gt.nj)then
        write(*,*)'l1o',l1o,ni,m1o,nj
      stop 'err l1o>ni in interp'
	endif

	do 500 ii=1,l1
	write(*,*)ii
	do 500 jj=1,m1
c	write(*,*)ii,jj,kk
      xxx=x(ii,jj)
      yyy=y(ii,jj)
50	continue
c	write(*,*)id,delx,xxx,yyy,zzz
      ss_0=(dmin1(1./dfloat(l1o),aspect/dfloat(m1o)))**2
      ss_0=1d-6
      ss_min=1d-6
10	iis=0
	do 100 i1=1,l1o
	do 100 j1=1,m1o
      xx1=xo(i1,j1)
      yy1=yo(i1,j1)
        ss=(xx1-xxx)*(xx1-xxx)+(yy1-yyy)*(yy1-yyy)
        if(ss.le.ss_min)then
      do loop=1,nfch
      f(ii,jj,loop)=uch(i1,j1,loop)
	enddo
	goto 500
	endif

	if(ss.le.ss_0)then
	iis=iis+1
	ssd(iis)=ss
	is(iis)=i1
	js(iis)=j1
c     write(*,*)xo(i1,j1),yo(i1,j1)
	endif
        if(iis.gt.2)goto 101
	if(iis.gt.mss)then
	write(*,*)'iis=',iis,'ss_0=',ss_0
	stop 'err in interp'
	endif
100	continue

101	if(iis.lt.1)then
	ss_0=ss_0*1.5
	goto 10
	endif


c	write(*,*)x1,y1,z1
      do 400 loop=1,nfch
c	write(*,*)loop
	uuu=0
	sss=0
	do ii1=1,iis
	i1=is(ii1)
	j1=js(ii1)
c use distance
c        ss=1/dsqrt(ssd(ii1))
c use distance square
        ss=1/(ssd(ii1))

      uuu=uuu+uch(i1,j1,loop)*ss
	sss=sss+ss
	enddo
      f(ii,jj,loop)=uuu/sss
c	write(*,*)ii,jj,kk,xxx,yyy,zzz,uuu/sss
c	write(*,*)iis
400	continue
500	continue
	return
	end





      subroutine restart
	include 'sc.h'

      iter=0
        iun=1
        if(iun.eq.0)then
      open(unit=8, file='input', form='unformatted')
        else
        open(8,file='input')
        endif
	l1o=l1
	m1o=m1
        if(iun.eq.0)then
	read(8) tinit,l1o,m1o

        else
        read(8,*)tinit,l1o,m1o
        endif
	write(*,*)'tinit',tinit
        if(iun.eq.0)then
c      read(8)((x(i,j),i=1,l1o),j=1,m1o)
c      read(8)((y(i,j),i=1,l1o),j=1,m1o)
      read(8)((u(i,j),i=1,l1o),j=1,m1o)
      read(8)((u(i,j),i=1,l1o),j=1,m1o)
      read(8)((u(i,j),i=1,l1o),j=1,m1o)
      read(8)((v(i,j),i=1,l1o),j=1,m1o)
        else
      read(8,*)((u(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((u(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((u(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((v(i,j),i=1,l1o),j=1,m1o)
        endif
	goto 111
        if(iun.eq.0)then
      read(8)((w(i,j),i=1,l1o),j=1,m1o)
      read(8)((p(i,j),i=1,l1o),j=1,m1o)
      read(8)((t(i,j),i=1,l1o),j=1,m1o)
        else
      read(8,*)((w(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((p(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((t(i,j),i=1,l1o),j=1,m1o)


        endif
	goto 112
111	continue
        if(iun.eq.0)then
      read(8)((p(i,j),i=1,l1o),j=1,m1o)
      read(8)((t(i,j),i=1,l1o),j=1,m1o)
      read(8)((ccen(i,j),i=1,l1o),j=1,m1o)
        else
      read(8,*)((p(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((t(i,j),i=1,l1o),j=1,m1o)

c      read(8,*)((t(2*i,2*j),i=1,l1o),j=1,m1o)
c      do i=1,l1o
c      do j=1,m1o
c      t(2*i-1,2*j-1)=t(2*i,2*j)
c      t(2*i-1,2*j)=t(2*i,2*j)
c      t(2*i,2*j-1)=t(2*i,2*j)
c      enddo
c      enddo

      read(8,*)((ccen(i,j),i=1,l1o),j=1,m1o)
        endif
112	continue
        if(iun.eq.0)then
      read(8)((ruksi(i,j),i=1,l1o),j=1,m1o)
      read(8)((rueta(i,j),i=1,l1o),j=1,m1o)
      read(8)((f(i,j,10),i=1,l1o),j=1,m1o)
      read(8)((qth(i,j),i=1,l1o),j=1,m1o)
      read(8)((fc(i,j),i=1,l1o),j=1,m1o)
        else
      read(8,*)((ruksi(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((rueta(i,j),i=1,l1o),j=1,m1o)
      read(8,*)((f(i,j,10),i=1,l1o),j=1,m1o)
      read(8,*)((qth(i,j),i=1,l1o),j=1,m1o)

c (fc,f(1,1,8)
      read(8,*)((f(i,j,8),i=1,l1o),j=1,m1o)
      read(8,*)((f(i,j,9),i=1,l1o),j=1,m1o)

        endif
	close(8)
   95 format(e12.5)
   96 format(6e12.5)
	do i=1,l1
	do j=1,m1
c	rueta(i,j)=0
c	ruksi(i,j)=0
	do nf=1,4
c	f(i,j,nf)=0
	if(id_bound(i,j).eq.0)then
c	f(i,j,nf)=0
c	ccen(i,j)=0
	endif
	enddo
	enddo
	enddo


      if(l1.ne.l1o.or.m1.ne.m1o)then
      call change_grid
      endif

      return
      end

	subroutine cal_elemag
	include 'sc.h'
	do j=1,m1
	write(*,*)'conduc',y(l1/2,j),id_ele(l1/2,j),conduc(l1/2,j)
	enddo

        do 100 loop=1,30
c boudary condition
      do j=1,m1
      fc(1,j)=0
      fc(l1,j)=0
      enddo
      do i=1,l1
      fc(i,1)=0
      fc(i,m1)=0
      enddo
	call matr_ele
	write(*,*)'matr_ele',loop
	call solve_c

	fmax_r=0
	fmax_i=0
	fmax=0
	do i=1,l1
	do j=1,m1
	rsl=abs(dreal(fc(i,j))-dreal(fco(i,j)))
	if(rsl.ge.fmax_r)then
	fmax_r=rsl
	endif
	rsl=abs(dimag(fc(i,j))-dimag(fco(i,j)))
	if(rsl.ge.fmax_i)then
	fmax_i=rsl
	endif

	fco(i,j)=fc(i,j)
	fc0=abs(dreal(fc(i,j)))
	if(fmax.lt.fc0)fmax=fc0
	enddo
	enddo
	write(*,*)'fcmax=',fmax
	if(fmax_r.gt.1d-12.or.fmax_i.gt.1d-12)goto 100
	if(blc.gt.1d-12)goto 100
	goto 200
100	continue
200	write(*,*)fmax_r,fmax_i,loop
c	write(*,*)(a_mag(l1/2,j,1),j=1,m1)
c	write(*,*)(fc(l1/2,j),j=1,m1)
	write(*,*)omega
	fmax=0
	pow_c=0
	pow_c1=0
	do i=1,l1
	do j=1,m1
	qth0=0.5*conduc(i,j)*omega*omega
     1   *(cdabs(fc(i,j)))**2
	if(id_ele(i,j).ne.11)then
	qth(i,j)=qth0
	pow_c=pow_c+qth0*vol(i,j)
	if(fmax.lt.qth0)fmax=qth0
	endif
	pow_c1=pow_c1+qth0*vol(i,j)
	enddo
	enddo
	write(*,*)'qthmax=',fmax
	write(*,*)'power consuming=',pow_c
	write(*,*)'all power consuming=',pow_c1
	call zero(u,l1,m1)
	call zero(v,l1,m1)
	do i=2,l2
	do j=2,m2
	u(i,j)=(dreal(fc(i,j+1))-dreal(fc(i,j)))/(y(i,j+1)-y(i,j))
     1        +dreal(fc(i,j))/y(i,j)
	v(i,j)=-(dreal(fc(i+1,j))-dreal(fc(i,j)))/(x(i+1,j)-x(i,j))
	enddo
	enddo

c coil_xx0,coil_yy0,coil_dxx from radi_ini


	do i=1,l2
      if(abs(x(i,2)-(coil_xx0+coil_dxx*2)).lt.2d-2)then
	write(*,*)'x=',x(i,2),'for y, u(j+1)-u(j),u(j)'
	do j=1,m2
	write(*,'(3f10.6)')y(i,j),u(i,j+1)-u(i,j),u(i,j)
	enddo
	goto 111
	endif
	enddo
111	do j=1,m2
      if(abs(y(2,j)-coil_yy0).lt.0.5d-2)then
	write(*,*)'y=',y(2,j),'for x, u(i+1)-u(i),u(i)'
	do i=1,l2
	write(*,'(3f10.6)')x(i,j),u(i+1,j)-u(i,j),u(i,j)
	enddo
	goto 112
	endif
	enddo
112	open(11,file='mag.plt')
	write(11,82)
	write(11,83)l1,m1
	write(11,88)((y(i,j),i=1,l1),j=1,m1)
	write(11,88)((x(i,j),i=1,l1),j=1,m1)
	write(11,88)((dreal(fc(i,j)),i=1,l1),j=1,m1)
	write(11,88)((dimag(fc(i,j)),i=1,l1),j=1,m1)
	write(11,88)((qth(i,j),i=1,l1),j=1,m1)
	write(11,88)((conduc(i,j),i=1,l1),j=1,m1)
	write(11,88)((u(i,j),i=1,l1),j=1,m1)
	write(11,88)((v(i,j),i=1,l1),j=1,m1)
	close(11)
	call zero(u,l1,m1)
	call zero(v,l1,m1)
   82 format(1x,'variables ="x","y","mag1","mag2","qth","con","u","v"')
   83 format(1x,'zone t= "zone 1" , i=',i3,' , j=',i3,' , f= block')
   88 format(6(1pe12.5,1x))

	open(11,file='mag.dat')
	write(11,*)l1,m1
	write(11,88)((qth(i,j),i=1,l1),j=1,m1)
	write(11,88)((fc(i,j),i=1,l1),j=1,m1)
	close(11)

	xi0=0.
	xi1=0.05
	xi2=0.1
	xi3=0.15
	xi4=0.2
	i0=interid(xi0,x,l1)
	i1=interid(xi1,x,l1)
	i2=interid(xi2,x,l1)
	i3=interid(xi3,x,l1)
	i4=interid(xi4,x,l1)
	open(11,file='magy.dat')
	do j=1,m1
	write(11,'(6e10.3)')y(i0,j),qth(i0,j),qth(i1,j),qth(i2,j),
     1  qth(i3,j)
     1       ,qth(i4,j)
	enddo
	close(11)
	return
	end

	integer function insd(x,y,uu1,uu2,vv1,vv2)
	implicit real*8 (a-h,o-z)
	insd=0
	x1=uu1
	x2=uu2
	y1=vv1
	y2=vv2
	if(x1.gt.x2)then
	x2=uu1
	x1=uu2
	y2=vv1
	y1=vv2
	endif
	if(x.ge.x1.and.x.le.x2.and.y.ge.y1.and.y.le.y2)insd=1
	return
	end

	integer function insdh(x,y,uu1,uu2,vv1,vv2)
	implicit real*8 (a-h,o-z)

	insdh=0
	x1=uu1
	x2=uu2
	y1=vv1
	y2=vv2
	if(x1.gt.x2)then
	x2=uu1
	x1=uu2
	y2=vv1
	y1=vv2
	endif

	yy1=dmin1(y1,y2)
	yy2=y1+(x-x1)/(x2-x1)*(y2-y1)
	if(x.ge.x1.and.x.le.x2.and.y.ge.yy1.and.y.le.yy2)insdh=1
	return
	end

	integer function interid(x0,x,m)
	implicit real*8 (a-h,o-z)
	real*8 x(m)
	xx=1d9
	do i=1,m
	if(abs(x0-x(i)).lt.xx)then
	i0=i
	xx=abs(x0-x(i))
	endif
	enddo
	interid=i0
	return
	end

	real*8 function rinterp(x,y,n,x0)
	implicit real*8 (a-h,o-z)
	real*8 x(n),y(n),x0
	if(x(1).gt.x(n))stop 'err in rinterp'

	if(x0.lt.x(1))ii=1
	if(x0.gt.x(n))ii=n-1
	do i=1,n-1
	if(x0.ge.x(i).and.x0.le.x(i+1))then
	ii=i
	goto 10
	endif
	enddo
10	rinterp=y(ii)+(y(ii+1)-y(ii))*(x0-x(ii))/(x(ii+1)-x(ii))
	return
	end

	real*8 function viscosity(i,j)
	include 'sc.h'
	real*8 visc(8),targ(8)
c viscosity of argon, micro Pa s, versus celsus degree
	data visc/54.65,58.36,61.94,65.39,81.33,95.84,109.45,122.39/
	data targ/700,800,900,1000,1500,2000,2500,3000/
	id=id_ele(i,j)
	viscosity=1d20
	tt=t(i,j)-273
	if(id.eq.3.or.id.eq.4)then
	viscosity=rinterp(targ,visc,8,tt)*1d-6
c	viscosity=1
	endif
	return
	end

	real*8 function cal_diffu(t)
	implicit real*8 (a-h,o-z)
	dimension cn_heat(14),cn_heatt(14)
        dimension temd(14),heatc(14)
c graphite density : 2200 kg/m3
	data cn_heat/1.0333,.8867,.785,.6167,.50,.43,.3667,.3333,
     1    .3017,.2850,.2667,.2633,.2533,.2533/
	data cn_heatt/0,125,250,500,750,1000,1250,1500,
     1    1750,2000,2250,2500,2750,3000/
        data temd/298,300,400,500,600,700,800,
     1  900,1000,1100,1200,1300,1400,1500/
        data heatc/710,716,997,1210,1382,1524,1640,
     1     1734,1810,1863,1915,1949,1973,1991/
	tcel=t-273
	th_conduc=rinterp(cn_heatt,cn_heat,14,tcel)
	th_conduc=107.9975*th_conduc
	capac=rinterp(temd,heatc,14,t)
	alfa=th_conduc/2200/capac
c	write(*,*)'con,rho,cp,alfa',th_conduc,2200,capac,alfa
	write(*,*)t,th_conduc,alfa
	cal_diffu=alfa
	return
	end

	real*8 function th_conduc(i,j)
c for sterling
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: foam 
c 7:felt  8:graphite stick 9: crucible 10: air with condution
c 11: copper 12: bottom graphite 13: crystal
c 14: gas with apparent conductivity 50 W/m/K

	include 'sc.h'
	dimension gam_d(12)
	dimension cnabh(11)
	dimension cnabs(13),cnt_ad(13),cnch(11),cncs(11),grh(11)
	dimension cn_felt(8),cn_feltt(8)
	dimension cn_foamt(10),
     1    cn_foamstdc_vac(8),cn_foamd_vac(8),
     1   cn_foamstdb(10),cn_foamstdc(10),
     1   cn_foamd(10),cn_foamhd(10)
	dimension cn_heat(14),cn_heatt(14),cn_gra(14),cn_grat(14)

      dimension tem_arg(10),con_arg(10)
      dimension tem_sic(11),con_sic(11)
c thermal conductivity of sic  in celsius
      data tem_sic/0,250,500,750,1000,1250,1500,
     1     1750,2000,2250,2500/
      data con_sic/650,220,120,110,100,90,85,75,65,60,50/


c thermal conductivity of copper 4.3744x10^^6 /ohm/m
	data gam_d/20,2d0,5d0,0.05d0,2d0,0.1,0,50,50,50,0,7/
c in celsius
	data tem_arg/100,200,400,600,
     1   800,1000,1500,2000,2500,3000/
	data con_arg/0.02133,0.02564,0.03318,0.03976,
     1  0.04572,0.05123,0.06373,0.07508,0.08573,0.09586/

c hard felt ab
	data cnabh/0.4,0.5,0.8,0.9,1.1,1.8,2.0,3.0,3.5,5.0,6.5/
	data cnch/0.18,0.20,0.3,0.4,0.5,0.7,0.9,1.0,1.7,2.,3./

	data cnt_ad/0,250,500,750,1000,1250,1500,1750,2000,2250,
     1    2500,2750,3000/
c soft felt ab
	data cnabs/0.09,0.105,0.19,0.22,0.32,0.5,0.7,1.,1.2,2.,
     1    3., 5., 7./
	data cncs/0.04,0.05,0.08,0.1,0.18,0.2,0.3,0.45,0.65,0.95,1.1/
c graphite
	data grh/60.,60.,40.,35.,30.,25.,22.,20.,20.,20.,20./
c felt for Sterling, graphite (argon) (in celsius degree)
	data cn_felt/0.22,0.2676,0.3190,0.4190,0.5613,0.8442,1.297,2.02/
	data cn_feltt/538,   816,  1093,  1371,  1649,  1927, 2204,2482/
c foam conductivity for Sterling, (in celsius degree)
	data cn_foamt/    400, 800,1000,1200,1400,1600,1800,2000,
     1    2200,2400/
	data cn_foamstdc_vac/0.12,0.16,0.22,0.32,0.40,0.51,0.65,0.76/
	data cn_foamd_vac/   0.35,0.53,0.58,0.66,0.78,0.90,1.05,1.24/
	data cn_foamstdb/0.25,0.35,0.42,0.48,0.56,0.68,0.84,0.98,
     1     1.31,1.5/
	data cn_foamstdc/0.15,0.24,0.34,0.4,0.48,0.62,0.78,0.92,
     1    1.06,1.2/
	data cn_foamd/   0.55,0.68,0.75,0.85,0.95,1.10,1.28,1.45,
     1   1.61,1.8/
	data cn_foamhd/2.05,2.1,2.15,2.18,2.21,2.22,2.3,2.32,
     1     2.35,2.4/
c graphite conductivity for Sterling (in celsius)
	data cn_gra/1.,0.747,0.651,0.521,0.427,0.36,0.307,0.273,
     1      0.24,0.224,0.213,0.2,0.2,0.211/
	data cn_grat/0,125,250,500,750,1000,1250,1500,
     1     1750,2000,2250,2500,2750,3000/
c heater conductivity for Sterling (in celsius)
	data cn_heat/1.0333,.8867,.785,.6167,.50,.43,.3667,.3333,
     1    .3017,.2850,.2667,.2633,.2533,.2533/
	data cn_heatt/0,125,250,500,750,1000,1250,1500,
     1    1750,2000,2250,2500,2750,3000/

	th_conduc=0
	id=id_ele(i,j)

c default thermal conductivity
	th_conduc=0.

c for  thermal conductivity
	tcel=t(i,j)-273

	if(id.eq.8.or.id.eq.9.or.id.eq.12)then
	th_conduc=gam_d(id)
	endif

	if(id.eq.0.or.id.eq.4)then
	th_conduc=rinterp(tem_arg,con_arg,10,tcel)
	endif
c increase thermal cond. for argon gas around thermal pack
c for th_conduct=1, t_top=1775, t_bottom=2290
        if(id.eq.4)th_conduc=2

c from Pons, J. Electrochem. Soc., Vol. 143, p. 3727
        if(id.eq.14) th_conduc=50

c for air
	if(id.eq.10)then
	th_conduc=5
	tref=t(i,j)
      th_conduc=0.9*5.67d-8*tref*tref*tref*4*air_len
	endif


	if(id.eq.3)then
	void_s=0.4
	emiss=0.3
	diam_p=0.0016
      con_sicp=5
	con_gas=0.05

	ste=5.67051d-8
      th_conduc=(1-void_s)*con_sicp+void_s*(con_gas
     1   +8./3.*emiss*4*ste*t(i,j)**3*diam_p)

c from Pons, J. Electrochem. Soc., Vol. 143, p. 3727
c        th_conduc=5
c	write(*,*)th_conduc
c	pause
	endif

	if(id.eq.13)then
c SiC
      th_conduc=50
      th_conduc=rinterp(tem_sic,con_sic,11,tcel)
	endif

	ee=(t(i,j)-273.)/250.+1
	iee=int(ee)
	if(iee.gt.10)iee=10
	if(iee.lt.1)iee=1
	if(iee.lt.1.or.iee.gt.11)then
	
	write(*,*)"iee,t(i,j),i,j=",iee,t(i,j),i,j
	stop 'err in th_conduc'
	endif

c for graphite heater
	if(id.eq.1)then
c	th_conduc=grh(iee)*(iee+1-ee)+grh(iee+1)*(ee-iee)
c	th_conduc=107.9975*(1./(0.115*t(i,j)**0.1)-3.95+0.0001*t(i,j))
	th_conduc=rinterp(cn_heatt,cn_heat,14,tcel)
	th_conduc=107.9975*th_conduc
	endif

c for soft felt of ATMI
	if(id.eq.2)then
c	th_conduc=cncs(iee)*(iee+1-ee)+cncs(iee+1)*(ee-iee)
c	th_conduc=cnabs(iee)*(iee+1-ee)+cnabs(iee+1)*(ee-iee)
	th_conduc=rinterp(cnt_ad,cnabs,13,tcel)
	endif

cfor hard felt of ATMI
	if(id.eq.5)then
c	th_conduc=cnch(iee)*(iee+1-ee)+cnch(iee+1)*(ee-iee)
	th_conduc=cnabh(iee)*(iee+1-ee)+cnabh(iee+1)*(ee-iee)
	endif

c for foam
	if(id.eq.6)then
c	th_conduc=0.15+(0.36-0.15)*(t(i,j)-273-200)/(1650-200)
c	th_conduc=0.17*1d-6*t(i,j)*t(i,j)+0.08
	th_conduc=rinterp(cn_foamt,cn_foamstdc,10,tcel)
c	th_conduc=rinterp(cn_foamt,cn_foamstdb,10,tcel)
c	th_conduc=rinterp(cn_foamt,cn_foamd,10,tcel)
c	th_conduc=rinterp(cn_foamt,cn_foamhd,10,tcel)
	endif

c for felt
	if(id.eq.7)then

c	th_conduc=0.07+(0.12-0.07)*(t(i,j)-273-50)/(1000-50)

c	th_conduc=rinterp(cn_feltt,cn_felt,8,tcel)
c use ATMI data for felt
	th_conduc=rinterp(cnt_ad,cnabs,13,tcel)
	endif

c for graphite stick in Sterling
	if(id.eq.8)then
	th_conduc=rinterp(cn_heatt,cn_heat,14,tcel)
	th_conduc=107.9975*th_conduc
	endif
c for crucible in Sterling
	if(id.eq.9)then
c	th_conduc=107.9975*(1./(0.11*t(i,j)**0.2)+1.2d-4*t(i,j)-1.99)
	th_conduc=rinterp(cn_grat,cn_gra,14,tcel)
      th_conduc=107.9975*th_conduc
	endif

	if(id.eq.12)then
	tref=t(i,j)
	th_conduc=5.67d-8*tref*tref*tref*4*0.002
	endif
	return
	end

	real*8 function rmu_c(i,j)
      implicit real*8 (a-h,o-z)
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: foam 
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11: copper 12: bottom graphite
	rmu_c=1
	id=id_ele(i,j)
	if(id.eq.2)then
	rmu_c=1
	endif
	if(id.eq.5)then
	rmu_c=1
	endif
	return
	end

	real*8 function heat_cap(i,j)
c for sterling
	include "sc.h"
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: griphe foam 
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11: copper 12: bottom graphite 13: crystal
c 14: gas with apparent conductivity 50 W/m/K

c temperature and heat capacity for graphite
        dimension temd(14),heatc(14),temd_sic(3),heatc_sic(3)

        data temd/298,300,400,500,600,700,800,
     1  900,1000,1100,1200,1300,1400,1500/
        data heatc/710,716,997,1210,1382,1524,1640,
     1     1734,1810,1863,1915,1949,1973,1991/
c heat capacity for sic powder (in Kelvin)
        data temd_sic/293,1273,1823/
        data heatc_sic/586,1297,1465/

	tkel=t(i,j)

	capac=0
	id=id_ele(i,j)
        if(id.eq.0.or.id.eq.4.or.id.eq.10.or.id.eq.14)then
c	capac=(20.8+51.7d-6*t(i,j)*t(i,j))/0.0399
	capac=520.326
	endif

	if(id.eq.1.or.id.eq.2.or.id.eq.5.or.id.eq.6.or.
     1    id.eq.7.or.id.eq.8.or.id.eq.9.or.id.eq.12)then
	capac=rinterp(temd,heatc,14,tkel)
c        do i1=1,13
c        if(tkel.ge.temd(i1).and.tkel.le.temd(i1+1))then
c        capac=heatc(i1)+(heatc(i1+1)-heatc(i1))/(temd(i1+1)-temd(i1))
c     1    *(tkel-temd(i1))
c        goto 10
c        endif
c        enddo
c        if(tkel.le.temd(1))capac=heatc(1)
c        if(tkel.ge.temd(14))capac=heatc(14)
	endif


        if(id.eq.3)then
c        if(tkel.le.temd_sic(1))capac=heatc_sic(1)
c        if(tkel.ge.temd_sic(3))capac=heatc_sic(3)
c        do i1=1,2
c        if(tkel.ge.temd_sic(i1).and.tkel.le.temd_sic(i1+1))then
c        capac=heatc_sic(i1)+(heatc_sic(i1+1)-heatc_sic(i1))
c     1   /(temd_sic(i1+1)-temd_sic(i1))
c     1    *(tkel-temd_sic(i1))
c        goto 10
c        endif
c        enddo
c10	continue
	capac=rinterp(temd_sic,heatc_sic,3,tkel)
        endif

	if(id.eq.13)then
	capac=rinterp(temd_sic,heatc_sic,3,tkel)
	endif

	heat_cap=capac
	return
	end

	real*8 function conduc(i,j)
c for sterling
	include "sc.h"
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: griphe foam 
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11: copper 12: bottom graphite 13: crystal
c 14: gas with apparent conductivity 50 W/m/K


c only consider conduction of SiC and graphite heater
c conductivy: /ohm/m
	dimension conduc_d(12)
c resistivity of hard felt (c): 220x10^-5 ohm.m from ATMI
c resistivity of soft felt : 200x10^-5 ohm.m from ATMI
        data conduc_d/1d5,0.5d3,0.2d4,0,454,0,500,0,0,0,4.3744d6,0/


	dimension con_sict(7),con_sic(7)
	data con_sic/1.,10.,0.5d2,0.2d3,0.4d3,1d3,0.17d4/
	data con_sict/1000,1250,1500,1750,2000,2250,2500/
	dimension con_gra(10),con_grat(10)
	data con_gra/11,9,8.3,8.3,8.5,8.6,9.1,9.6,9.9,10.3/
c in celsius degree
	data con_grat/0,250,500,750,1000,1250,1500,1750,2000,2250/
	conduc=0
	id=id_ele(i,j)
	if(id.eq.2.or.(id.ge.5.and.id.le.12))then
	conduc=conduc_d(id)
	endif

	tcel=t(i,j)-273
c	if(id_ele(i,j).eq.11)conduc=0

	ee=(t(i,j)-273.)/250.+1
	iee=int(ee)

	if(iee.gt.9)iee=9
	if(iee.lt.1)iee=1

	if(id.eq.1.or.id.eq.8)then
	conduc=rinterp(con_grat,con_gra,10,tcel)

c	if(t(i,j).ge.2250+273)then
c	conduc=10.3
c	else
c	conduc=con_gra(iee)*(iee+1-ee)+con_gra(iee+1)*(ee-iee)
c	endif

	conduc=conduc*1d4
	endif

	if(id.eq.3)then
	temp=t(i,j)-273
	if(temp.le.1000)then
	conduc=con_sic(1)
	else if(temp.ge.2500)then
	conduc=con_sic(7)
	else
	ee=(temp-1000)/250.+1
	iee=int(ee)
	conduc=con_sic(iee)*(1+iee-ee)+con_sic(iee+1)*(ee-iee)
	endif
	endif

	if(id.eq.13)then
	conduc=rinterp(con_sict,con_sic,7,tcel)
	endif
	return
	end

	integer function id_bound_f(i1,j1)
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: griphe foam 
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11: copper
	include 'sc.h'
	id=id_ele(i1,j1)

	if(nf.ne.5)then
	if(id.eq.3.or.id.eq.4)then
	id_bound_f=1
	else
	id_bound_f=0
	endif
	endif

	if(nf.eq.5)then
	if(id.eq.0.or.id.eq.11)then
	id_bound_f=0
	else
	id_bound_f=1
	endif
	endif
	return
	end



	integer function id_ele(i1,j1)
c for sterling
	include 'sc.h'
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: graphite foam 
c 7:felt  8:temperature constant 9: crucible 10: air with condution
c 11: copper, 12: bottom graphite under crucible, 13: crystal
c 14: gas with apparent conductivity 50 W/m/K


c 2+0.75+1.5
c 0.4375

      if(.not.(idesign.eq.3))goto 100

c  design of 1.2'' crystal and 2'' crucible

c for design III
	id_el=0
	xx=x(i1,j1)
	yy=y(i1,j1)
        do i9=1,id_total
        if(insd(xx,yy,xx_0(i9),xx_1(i9),yy_0(i9),yy_1(i9)).eq.1)then
        id_el=id_el0(i9)
        if(id_el.eq.10)air_len=air_len0(i9)
        endif
        enddo
	id_ele=id_el
	return



	if(insd(xx,yy,-0.080d0,0.280d0,0.d0,0.182d0).eq.1)id_el=4
c felt
	if(insd(xx,yy,-0.055d0,0.241d0,0.105d0,0.125d0).eq.1)id_el=7
c foam
	if(insd(xx,yy,0d0,0.226d0,0.084d0,0.105d0).eq.1)id_el=6
	if(insd(xx,yy,0.211d0,0.226d0,0.008d0,0.083d0).eq.1)id_el=6
	if(insd(xx,yy,0.226d0,0.241d0,0.008d0,0.105d0).eq.1)id_el=6
	if(insd(xx,yy,-0.055d0,0d0,0.0075d0,0.1225d0).eq.1)id_el=6
c air
      if(insd(xx,yy,0.0d0,0.226d0,0.083d0,0.084d0).eq.1)then
      id_el=10
      air_len=0.007
      endif
      if(insd(xx,yy,0.0d0,0.205d0,0.d0,0.035d0).eq.1)then
      id_el=10
      air_len=0.007
      endif
c felt
	if(insd(xx,yy,0.205d0,0.211d0,0.008d0,0.084d0).eq.1)id_el=7
	if(insd(xx,yy,0.200d0,0.205d0,0.050d0,0.084d0).eq.1)id_el=7
c felt near outer wall of heater
	if(insd(xx,yy,0.d0,0.200d0,0.050d0,0.083d0).eq.1)id_el=7
c graphite heater
	if(insd(xx,yy,0d0,0.200d0,0.035d0,0.050d0).eq.1)id_el=1
c lower heater
	if(insd(xx,yy,-0.005d0,0.005d0,0.d0,0.050d0).eq.1)id_el=1
c upper heater
c make sure heater is contact with stick
      if(insd(xx,yy,0.197d0,0.205d0,0.008d0,0.0085d0).eq.1)then
      id_el=10
      air_len=0.007
      endif
	if(insd(xx,yy,0.197d0,0.205d0,0.0085d0,0.050d0).eq.1)id_el=1
c felt above crucible
	if(insd(xx,yy,0.184d0,0.197d0,0.015d0,0.035d0).eq.1)id_el=6
c for temperature stick
	if(insd(xx,yy,0.184d0,0.224d0,0.0d0,0.008d0).eq.1)id_el=1
	if(insd(xx,yy,0.184d0,0.193d0,0.0d0,0.015d0).eq.1)id_el=1
	if(insdh(xx,yy,0.193d0,0.199d0,0.015d0,0.008d0).eq.1)id_el=1

c  bottom graphite
c	if(insd(xx,yy,0.005d0,.047d0,.0d0,.032d0).eq.1)id_el=7
	if(insd(xx,yy,0.005d0,.047d0,.0d0,.010d0).eq.1)id_el=9
	if(insd(xx,yy,0.041d0,.047d0,.0d0,.032d0).eq.1)id_el=9

	if(insd(xx,yy,0.005d0,.047d0,.010d0,.032d0).eq.1)id_el=12
c for crucible

	if(insd(xx,yy,.045d0,0.175d0,0.0d0,.028d0).eq.1)id_el=9
	if(insd(xx,yy,.164d0,0.184d0,0.0d0,.0305d0).eq.1)id_el=9

c for graphite near seed
	if(insd(xx,yy,0.175d0,0.181d0,0.0d0,0.028d0).eq.1)id_el=1
c for felt near seed
	if(insd(xx,yy,0.176d0,0.181d0,0.015d0,0.028d0).eq.1)id_el=7
c crystal 
      if(insd(xx,yy,0.175d0-crys_len,0.175d0,.0d0,.025d0).eq.1)
     1   id_el=13
c argon
	if(insd(xx,yy,0.128d0,.175d0-crys_len,.0d0,.025d0).eq.1)id_el=4
c for powder
	if(insd(xx,yy,.048d0,0.128d0,0.0d0,.025d0).eq.1)id_el=3
	

c for part under foam
c	if(insd(xx,yy,.00976d0,0.04589d0,0.0d0,.00814d0).eq.1)id_el=6
c	if(insd(xx,yy,.00976d0,0.02034d0,0.0d0,.01464d0).eq.1)id_el=6
	
C For air upon graphite stick
	if(insd(xx,yy,.224d0,0.241d0,0.0d0,0.008d0).eq.1)id_el=4
c for air below heater
	if(insd(xx,yy,-.055d0,-0.005d0,0.0d0,.0075d0).eq.1)id_el=4

c	if(insd(xx,yy,-.055d0,0.0d0,0.0d0,.007973d0).eq.1)id_el=0
c	if(insd(xx,yy,.0d0,0.04589d0,0.0d0,.00651d0).eq.1)id_el=0

c	if(insd(xx,yy,0.2107d0,0.240d0,0.0d0,.00960d0).eq.1)id_el=0
c	if(insd(xx,yy,0.2146d0,0.240d0,0.0d0,.03254d0).eq.1)id_el=0
c for copper


	do ii=1,5
      xxx=coil_xx0+coil_dxx*(ii-1)
c	if(insd(xx,yy,xxx-0.01905,xxx+0.01905,0.196d0,.205525d0)
c     1  .eq.1)id_el=11
c	if(insd(xx,yy,xxx-0.01905+0.001219,xxx+0.01905-0.001219,
c     1  0.196d0+0.001219,.205525d0-0.001219).eq.1)id_el=0
	enddo

	id_ele=id_el
	return
100	continue

      if(.not.(idesign.eq.2))goto 200
c This for design of 2.5'' crystal and 3'' crucible
c for design II
	id_el=0
	xx=x(i1,j1)
	yy=y(i1,j1)
      if(insd(xx,yy,-0.080d0,0.280d0+felt_incr,
     1   0.d0,0.182d0).eq.1)id_el=4
c foam
      if(insd(xx,yy,0.240d0+felt_incr,0.270d0+felt_incr,
     1    0.009d0,0.105d0).eq.1)id_el=6
      if(insd(xx,yy,0d0,0.255d0+felt_incr,0.085d0,0.105d0).eq.1)
     1  id_el=6
	if(insd(xx,yy,-0.055d0,0d0,0.0075d0,0.1225d0).eq.1)id_el=6
c felt
      if(insd(xx,yy,0.d0,0.255d0+felt_incr,0.105d0,0.125d0).eq.1)
     1 id_el=7
      if(insd(xx,yy,0.207d0,0.232d0+felt_incr,0.009d0,.085d0)
     1    .eq.1)id_el=7
      if(insd(xx,yy,0.d0,.207d0,0.075d0,0.085d0).eq.1)id_el=7

c air
      if(insd(xx,yy,0.005d0,0.203d0,0.d0,0.060d0).eq.1)then
      id_el=10
      air_len=0.009
      endif
        if(insd(xx,yy,0.232d0+felt_incr,0.240d0+felt_incr,
     1   0.0d0,0.085d0).eq.1)then
      id_el=10
      air_len=0.008
      endif

c graphite susceptor
      if(insd(xx,yy,0d0,0.200d0,0.060d0,0.075d0).eq.1)id_el=1
c lower heater
      if(insd(xx,yy,-0.007d0,0.005d0,0.d0,0.075d0).eq.1)id_el=1
c upper heater
      if(insd(xx,yy,0.195d0,0.207d0,0.023d0,0.075d0).eq.1)id_el=1
c no bottom felt
c	if(insd(xx,yy,0.04589d0,.05533d0,.0d0,.03336d0).eq.1)id_el=7
c upper felt

c  bottom graphite
	if(insd(xx,yy,0.005d0,.050d0,.0d0,.010d0).eq.1)id_el=9
	if(insd(xx,yy,0.046d0,.050d0,.0d0,.056d0).eq.1)id_el=9

	if(insd(xx,yy,0.005d0,.050d0,.010d0,.056d0).eq.1)id_el=12

c for crucible

	if(insd(xx,yy,.050d0,0.193d0,0.0d0,.051d0).eq.1)id_el=9
	if(insd(xx,yy,.153d0,0.193d0,0.0d0,.055d0).eq.1)id_el=9
c argon
      if(insd(xx,yy,0.140d0,.178d0,.0d0,rc).eq.1)id_el=4

c for argon with apparent conductivity from Pons 1996
c        if(insd(xx,yy,0.140d0,.178d0,.0d0,rc).eq.1)id_el=14
c for powder
      if(insd(xx,yy,0.058d0,0.140d0,0.0d0,rc).eq.1)id_el=3
c for felt near seed
      if(insd(xx,yy,0.1795d0,0.190d0,0.020d0,0.051d0).eq.1)id_el=7
c      if(insd(xx,yy,0.165d0,0.190d0,0.040d0,0.0435d0).eq.1)id_el=7

c for temperature constant
	if(insd(xx,yy,.193d0,0.203d0,0.0d0,.020d0).eq.1)id_el=8
      if(insd(xx,yy,0.203d0,zb_c_t,0.0d0,.008d0).eq.1)id_el=8
      if(insdh(xx,yy,0.203d0,0.213d0,0.020d0,0.008d0).eq.1)id_el=8


c for part under foam
c	if(insd(xx,yy,.00976d0,0.04589d0,0.0d0,.00814d0).eq.1)id_el=6
c	if(insd(xx,yy,.00976d0,0.02034d0,0.0d0,.01464d0).eq.1)id_el=6
	
C For air upon temperature
      if(insd(xx,yy,zb_c_t,zt_c_t,0.0d0,rr_c_t).eq.1)id_el=4
	if(insd(xx,yy,-.055d0,-0.007d0,0.0d0,0.0075d0).eq.1)id_el=4

c	if(insd(xx,yy,.0d0,0.04589d0,0.0d0,.00651d0).eq.1)id_el=0
c	if(insd(xx,yy,0.2107d0,0.240d0,0.0d0,.00960d0).eq.1)id_el=0
c	if(insd(xx,yy,0.2146d0,0.240d0,0.0d0,.03254d0).eq.1)id_el=0
c for copper


	do ii=1,5
      xxx=coil_xx0+coil_dxx*(ii-1)
c	if(insd(xx,yy,xxx-0.01905,xxx+0.01905,0.196d0,.205525d0)
c     1  .eq.1)id_el=11
c	if(insd(xx,yy,xxx-0.01905+0.001219,xxx+0.01905-0.001219,
c     1  0.196d0+0.001219,.205525d0-0.001219).eq.1)id_el=0
	enddo

	id_ele=id_el
	return

200	continue
      if(.not.(idesign.eq.1))goto 300
c design of 
c for design 1
c 2+0.75+1.5
c 0.4375
	id_el=0
        xx=x(i1,j1)
        yy=y(i1,j1)
c foam
	if(insd(xx,yy,-0.055d0,0.240d0,0.00d0,0.1269d0).eq.1)id_el=6
	if(insd(xx,yy,0.00976d0,0.2146d0,0.00d0,0.03588d0)
     1   .eq.1)then
      id_el=10
      air_len=0.007
      endif
c graphite heater
	if(insd(xx,yy,0d0,0.20340d0,0.03588d0,0.05134d0).eq.1)id_el=1
c felt
	if(insd(xx,yy,0.d0,0.240d0,0.1087d0,0.1269d0).eq.1)id_el=7
c bottom felt
	if(insd(xx,yy,0.04589d0,.05533d0,.0d0,.03336d0).eq.1)id_el=7
c upper felt
	if(insd(xx,yy,0.1904d0,.2001d0,0.01139d0,.03588d0)
     1    .eq.1)id_el=7
	if(insd(xx,yy,0.2034d0,.2146d0,.05134d0,.08771d0).eq.1)id_el=7
	if(insd(xx,yy,0.2083d0,.2146d0,.0d0,.05134d0).eq.1)id_el=7

c for crucible
	if(insd(xx,yy,.1709d0,0.1904d0,0.0d0,.03157d0).eq.1)id_el=9
	if(insd(xx,yy,.05533d0,0.1904d0,0.0d0,.02880d0).eq.1)id_el=9
c for powder
	if(insd(xx,yy,.05858d0,0.1455d0,0.0d0,.02262d0).eq.1)id_el=3
c argon
	if(insd(xx,yy,0.1455d0,.1839d0,.0d0,.02262d0).eq.1)id_el=4

c for temperature constant
	if(insd(xx,yy,.1904d0,0.2107d0,0.0d0,.00822d0).eq.1)id_el=8
	if(insd(xx,yy,.1904d0,0.2001d0,0.0d0,.01139d0).eq.1)id_el=8

c for part under felt 
	if(insd(xx,yy,.00976d0,0.04589d0,0.0d0,.00814d0).eq.1)id_el=6
	if(insd(xx,yy,.00976d0,0.02034d0,0.0d0,.01464d0).eq.1)id_el=6
	
c for air
	if(insd(xx,yy,-.055d0,0.0d0,0.0d0,.007973d0).eq.1)id_el=0
	if(insd(xx,yy,.0d0,0.04589d0,0.0d0,.00651d0).eq.1)id_el=0

	if(insd(xx,yy,0.2107d0,0.240d0,0.0d0,.00960d0).eq.1)id_el=0
	if(insd(xx,yy,0.2146d0,0.240d0,0.0d0,.03254d0).eq.1)id_el=0
c for copper
c coil_xx0,coil_dxx from radi_ini

	do ii=1,5

      xxx=coil_xx0+coil_dxx*(ii-1)
      yyy=coil_yy0

c        xx_r=0.01905
c        yy_r=0.0047625

      if(insd(xx,yy,xxx-xx_r,xxx+xx_r,yyy-yy_r,yyy+yy_r)
     1  .eq.1)id_el=11
c        xx_r=0.01905-width_coil
c        yy_r=0.0047625-width_coil

      if(insd(xx,yy,xxx-xx_r+width_coil,xxx+xx_r-width_coil,
     1  yyy-yy_r+width_coil,yyy+yy_r-width_coil).eq.1)id_el=0
	enddo

	id_ele=id_el

300	continue


        if(.not.(idesign.eq.4))goto 400
c This is for design of a 3'' crystal and 3.5'' crucible
c for design II
	id_el=0
	xx=x(i1,j1)
	yy=y(i1,j1)
        if(insd(xx,yy,-0.080d0,0.280d0+z_incr,0.d0,0.182d0).eq.1)
     1  id_el=4
c upper foam
        if(insd(xx,yy,z_incr+0.240d0,0.270d0+z_incr,tos_a,0.105d0)
     1  .eq.1)id_el=6
c side foam
        if(insd(xx,yy,0d0,0.255d0+z_incr,0.085d0,0.105d0).eq.1)
     1   id_el=6
	if(insd(xx,yy,-0.055d0,0d0,0.0075d0,0.1225d0).eq.1)id_el=6
c felt
        if(insd(xx,yy,0.d0,0.255d0+z_incr,0.105d0,0.125d0).eq.1)
     1   id_el=7
c upper felt
      if(insd(xx,yy,z_incr+0.207d0,.232d0+z_incr,0.00d0,.085d0)
     1    .eq.1)id_el=7
c felt outside heater
        if(insd(xx,yy,0.d0,.207d0+z_incr,.075d0,.085d0).eq.1)id_el=7

c air
        if(insd(xx,yy,0.0d0,0.207d0+z_incr,0.d0,0.060d0).eq.1)then
        id_el=10
        air_len=0.007
        endif
        if(insd(xx,yy,z_incr+0.232d0,0.240d0+z_incr,0.00d0,0.085d0)
     1    .eq.1)then
      id_el=10
      air_len=0.007
      endif

c graphite heater
        if(insd(xx,yy,0d0,0.200d0+z_incr,0.060d0,0.075d0).eq.1)
     1  id_el=1
c lower heater
        if(insd(xx,yy,-0.007d0,0.005d0,0.d0,0.075d0).eq.1)id_el=1
c upper heater
        if(insd(xx,yy,z_incr+0.195d0,z_incr+0.207d0,
     1   tos_d+0.002d0,0.075d0).eq.1) id_el=1

c  bottom graphite
        if(insd(xx,yy,0.005d0,.039d0,.0d0,.010d0).eq.1)id_el=9
        if(insd(xx,yy,0.035d0,.041d0,.0d0,.056d0).eq.1)id_el=9

        if(insd(xx,yy,0.005d0,.041d0,.010d0,.056d0).eq.1)id_el=12

c for crucible

        if(insd(xx,yy,.039d0,0.193d0+z_incr,0.0d0,.051d0).eq.1)id_el=9
        if(insd(xx,yy,z_incr+0.153d0,0.193d0+z_incr,0.0d0,.055d0)
     1  .eq.1) id_el=9
c argon
        if(insd(xx,yy,0.140d0,.178d0+z_incr,.0d0,.0455d0).eq.1)id_el=4

c for argon with apparent conductivity from Pons 1996
c        if(insd(xx,yy,0.140d0,.178d0+z_incr,.0d0,.0385d0).eq.1)id_el=14
c for powder
        if(insd(xx,yy,.044d0,0.140d0,0.0d0,.0455d0).eq.1)id_el=3
c for felt near seed
        if(insd(xx,yy,0.1795d0+z_incr,0.190d0+z_incr,tos_d,0.051d0)
     1  .eq.1)id_el=7

c for temperature constant
      if(insd(xx,yy,.193d0+z_incr,z_incr+0.193d0+tos_h2,0.0d0,tos_d)
     1   .eq.1)id_el=8
      if(insd(xx,yy,z_incr+0.193d0+tos_h2,
     1   z_incr+0.193d0+tos_h,0.0d0,tos_a).eq.1) id_el=8
      if(insdh(xx,yy,z_incr+0.193d0+tos_h2,
     1   z_incr+0.193d0+tos_h-tos_h1,tos_d,tos_a).eq.1)id_el=8
c to check the effect of the shape of TOS
     



c for part under foam
c	if(insd(xx,yy,.00976d0,0.04589d0,0.0d0,.00814d0).eq.1)id_el=6
c	if(insd(xx,yy,.00976d0,0.02034d0,0.0d0,.01464d0).eq.1)id_el=6
	
C For air upon temperature
        if(insd(xx,yy,z_incr+.232d0,z_incr+0.270d0,0.0d0,tos_a)
     1   .eq.1)id_el=4
        if(insd(xx,yy,z_incr+0.193d0+tos_h,
     1   z_incr+0.232d0,0.0d0,tos_a).eq.1)id_el=4
	if(insd(xx,yy,-.055d0,-0.007d0,0.0d0,0.0075d0).eq.1)id_el=4

c	if(insd(xx,yy,.0d0,0.04589d0,0.0d0,.00651d0).eq.1)id_el=0
c       if(insd(xx,yy,z_incr+0.2107d0,z_incr+0.240d0,0.0d0,.00960d0)
c     1 .eq.1)id_el=0
c       if(insd(xx,yy,z_incr+0.2146d0,z_incr+0.240d0,0.0d0,.03254d0)
c      1  .eq.1)id_el=0
c for copper

	do ii=1,5
      xxx=coil_xx0+coil_dxx*(ii-1)
c	if(insd(xx,yy,xxx-0.01905,xxx+0.01905,0.196d0,.205525d0)
c     1  .eq.1)id_el=11
c	if(insd(xx,yy,xxx-0.01905+0.001219,xxx+0.01905-0.001219,
c     1  0.196d0+0.001219,.205525d0-0.001219).eq.1)id_el=0
	enddo

	id_ele=id_el
	return

400     continue
      if(idesign.ne.5)goto 500
c this is for design of a 4'' crystal and 4.5'' crucible

	id_el=0
	xx=x(i1,j1)
	yy=y(i1,j1)
      if(insd(xx,yy,-0.080d0,z_incr+0.280d0,0.d0,0.222d0).eq.1)id_el=4
c top foam
      if(insd(xx,yy,z_incr+0.240d0,z_incr+0.270d0,0.009d0,0.145d0)
     1  .eq.1)id_el=6
c side foam
      if(insd(xx,yy,0d0,z_incr+0.255d0,0.125d0,0.145d0).eq.1)id_el=6
c bottom foam
      if(insd(xx,yy,-0.055d0,0d0,0.0075d0,0.1625d0).eq.1)id_el=6
c felt
      if(insd(xx,yy,0.d0,z_incr+0.255d0,0.145d0,0.165d0).eq.1)id_el=7
c upper felt
      if(insd(xx,yy,z_incr+0.207d0,z_incr+0.232d0,0.00d0,0.125d0)
     1    .eq.1)id_el=7
c felt outside heater
        if(insd(xx,yy,0.d0,z_incr+0.207d0,.095d0,.125d0).eq.1)id_el=7

c air
        if(insd(xx,yy,z_incr+0.232d0,z_incr+0.240d0,0.00d0,0.125d0)
     1  .eq.1)then
      id_el=10
      air_len=0.007
      endif
        if(insd(xx,yy,0.0d0,z_incr+0.207d0,0.d0,0.075d0).eq.1)then
        id_el=10
      air_len=0.007
        endif

c graphite heater
        if(insd(xx,yy,0d0,z_incr+0.200d0,0.075d0,0.095d0).eq.1)id_el=1
c lower heater
        if(insd(xx,yy,-0.007d0,0.005d0,0.d0,0.095d0).eq.1)id_el=1
c upper heater
        if(insd(xx,yy,z_incr+0.195d0,z_incr+0.207d0,
     1   tos_d+0.002d0,0.095d0).eq.1) id_el=1

c  bottom graphite
        if(insd(xx,yy,0.005d0,.039d0,.0d0,.010d0).eq.1)id_el=9
        if(insd(xx,yy,0.035d0,.041d0,.0d0,.068d0).eq.1)id_el=9

        if(insd(xx,yy,0.005d0,.041d0,.010d0,.068d0).eq.1)id_el=12

c for crucible

        if(insd(xx,yy,.040d0,z_incr+0.193d0,0.0d0,.064d0).eq.1)id_el=9
      if(insd(xx,yy,z_incr+0.153d0,z_incr+0.193d0,0.0d0,.068d0)
     1   .eq.1)id_el=9
c argon
        if(insd(xx,yy,0.140d0,z_incr+0.178d0,.0d0,.060d0).eq.1)id_el=4

c for argon with apparent conductivity from Pons 1996
c        if(insd(xx,yy,0.140d0,z_incr+0.178d0,.0d0,.0385d0)
c     1.eq.1)id_el=14
c for powder
        if(insd(xx,yy,.044d0,0.140d0,0.0d0,.060d0).eq.1)id_el=3
c for felt near seed
        if(insd(xx,yy,z_incr+0.1795d0,z_incr+0.190d0,tos_d,0.064d0)
     1   .eq.1)id_el=7

c for temperature constant
      if(insd(xx,yy,z_incr+0.193d0,z_incr+0.193d0+tos_h2,0.0d0,tos_d)
     1   .eq.1)id_el=8
      if(insd(xx,yy,z_incr+0.193d0+tos_h2,z_incr+0.193d0+tos_h,
     1   0.0d0,tos_a).eq.1) id_el=8
      if(insdh(xx,yy,z_incr+0.193d0+tos_h2,z_incr+0.193d0+tos_h-tos_h1,
     1  tos_d,tos_a).eq.1)id_el=8
c to check the effect of the shape of TOS
     



	
C For air upon temperature
        if(insd(xx,yy,z_incr+0.232d0,z_incr+0.270d0,0.0d0,0.009d0)
     1   .eq.1)id_el=4
        if(insd(xx,yy,z_incr+0.193d0+tos_h,z_incr+0.232d0,0.0d0,tos_a)
     1  .eq.1)id_el=4
	if(insd(xx,yy,-.055d0,-0.007d0,0.0d0,0.0075d0).eq.1)id_el=4


	id_ele=id_el
	return



500   continue

        if(.not.(idesign.eq.6))goto 600
c This is for design of a 3'' system
c for design of 02 Furnace 3'' crucible and heat pack as of 16 Feb 2000
	id_el=0
	xx=x(i1,j1)
	yy=y(i1,j1)
c air 
        if(insd(xx,yy,-0.080d0,0.280d0+z_incr,0.d0,0.182d0).eq.1)
     1  id_el=4

c air with effective conductivity
        if(insd(xx,yy,0.0d0,0.220d0+z_incr,0.d0,0.060d0).eq.1)then
        id_el=10
      air_len=0.007
        endif
c upper foam
        if(insd(xx,yy,z_incr+0.240d0,0.270d0+z_incr,0.0225d0,0.105d0)
     1  .eq.1)id_el=6
c side foam
        if(insd(xx,yy,0d0,0.255d0+z_incr,0.085d0,0.105d0).eq.1)id_el=6
	if(insd(xx,yy,-0.055d0,0d0,0.0075d0,0.1225d0).eq.1)id_el=6
c felt
        if(insd(xx,yy,0.d0,0.255d0+z_incr,0.105d0,0.1225d0)
     1   .eq.1)id_el=7
c upper felt
      if(insd(xx,yy,z_incr+0.207d0,0.240d0+z_incr,0.015d0,.085d0)
     1    .eq.1)id_el=7
c felt outside heater
        if(insd(xx,yy,0.d0,.207d0+z_incr,.075d0,.085d0).eq.1)id_el=7


c graphite heater
        if(insd(xx,yy,0d0,0.200d0+z_incr,0.060d0,0.075d0).eq.1)id_el=1
c lower heater
        if(insd(xx,yy,-0.007d0,0.005d0,0.d0,0.075d0).eq.1)id_el=1
c upper heater
        if(insd(xx,yy,z_incr+0.195d0,z_incr+0.207d0,
     1   tos_d+0.002d0,0.075d0).eq.1) id_el=1

c  bottom graphite
        if(insd(xx,yy,0.005d0,.0365d0,.0d0,.020d0).eq.1)id_el=9
        if(insd(xx,yy,0.0325d0,.0365d0,.0d0,.055d0).eq.1)id_el=9

        if(insd(xx,yy,0.005d0,.0365d0,.020d0,.055d0).eq.1)id_el=12

c for crucible

        if(insd(xx,yy,.0365d0,0.0365d0+0.1495d0+z_incr,0.0d0,.051d0)
     1  .eq.1)id_el=9
        if(insd(xx,yy,0.0365d0+0.1095d0+z_incr,
     1    0.0365d0+0.1495d0+z_incr,0.0d0,.055d0).eq.1) id_el=9
c argon
        if(insd(xx,yy,0.0365d0+0.095d0,0.0365d0+0.140d0+z_incr,
     1    .0d0,.046d0).eq.1)id_el=4

c for argon with apparent conductivity from Pons 1996
c        if(insd(xx,yy,0.1315d0,.1765d0+z_incr,.0d0,.046d0).eq.1)id_el=14
c for powder
        if(insd(xx,yy,0.0365d0+0.005d0,0.0365d0+0.095d0,
     1   0.0d0,.046d0).eq.1)id_el=3
c for felt near seed
        if(insd(xx,yy,0.0365d0+0.142d0+z_incr,
     1   0.0365d0+0.1465d0+z_incr,tos_d,0.051d0)
     1  .eq.1)id_el=7

c for temperature constant
      if(insd(xx,yy,0.186d0+z_incr,z_incr+0.186d0+tos_h2,0.0d0,tos_d)
     1   .eq.1)id_el=8
      if(insd(xx,yy,z_incr+0.186d0+tos_h2,
     1   z_incr+0.186d0+tos_h,0.0d0,tos_a).eq.1) id_el=8
      if(insdh(xx,yy,z_incr+0.186d0+tos_h2,
     1  z_incr+0.186d0+tos_h-tos_h1,
     1  tos_d,tos_a).eq.1)id_el=8
c to check the effect of the shape of TOS
     
	
C For air upon temperature
        if(insd(xx,yy,z_incr+0.240d0,z_incr+0.270d0,0.0d0,0.0225d0)
     1   .eq.1)id_el=4
        if(insd(xx,yy,z_incr+0.186d0+tos_h,z_incr+0.240d0,0.0d0,tos_a)
     1  .eq.1)id_el=4
      if(insd(xx,yy,-.055d0,-0.007d0,0.0d0,0.005d0).eq.1)id_el=4

c	if(insd(xx,yy,.0d0,0.04589d0,0.0d0,.00651d0).eq.1)id_el=0
c       if(insd(xx,yy,z_incr+0.2107d0,z_incr+0.240d0,0.0d0,.00960d0)
c     1 .eq.1)id_el=0
c       if(insd(xx,yy,z_incr+0.2146d0,z_incr+0.240d0,0.0d0,.03254d0)
c      1  .eq.1)id_el=0
c for copper

	do ii=1,5
      xxx=coil_xx0+coil_dxx*(ii-1)
c	if(insd(xx,yy,xxx-0.01905,xxx+0.01905,0.196d0,.205525d0)
c     1  .eq.1)id_el=11
c	if(insd(xx,yy,xxx-0.01905+0.001219,xxx+0.01905-0.001219,
c     1  0.196d0+0.001219,.205525d0-0.001219).eq.1)id_el=0
	enddo

	id_ele=id_el
	return

600     continue

        return
	end

	subroutine zero(u,ni,nj)
	implicit real*8 (a-h,o-z)
	dimension u(ni,nj)
	do i=1,ni
	do j=1,nj
	u(i,j)=0
	enddo
	enddo
	return
	end

	subroutine us_elemag
	include 'sc.h'

	complex*16 fi1,fil1,fj1,fjm1
	entry matr_ele

	call reset1
	do 110 j=2,m2
	diff=1/(0.5*hksi(2,j))
	aim(2,j)=diff*ak1(2,j)
	diff=1/(0.5*hksi(l2,j))
	aip(l2,j)=diff*ak1(l1,j)
	do 110 i=2,l3
	diff=1/(0.5*hksi(i,j)+0.5*hksi(i+1,j))
	aip(i,j)=diff*ak1(i+1,j)
	aim(i+1,j)=aip(i,j)
110	continue

	do 120 i=2,l2
	diff=1/(0.5*heta(i,2))
	ajm(i,2)=diff*ae1(i,2)
	diff=1/(0.5*heta(i,m2))
	ajp(i,m2)=diff*ae1(i,m1)
	do 120 j=2,m3
	diff=1/(0.5*heta(i,j)+0.5*heta(i,j+1))
	ajp(i,j)=diff*ae1(i,j+1)
	ajm(i,j+1)=ajp(i,j)
120	continue

c boundary
	do j=2,m2
	goto 10
	fi1=
     1   (fc(2,j)-fc(3,j))
     1   /(x(2,j)-x(3,j))*0.5*hksi(2,j)
c	fi1=0
	con(2,j)=con(2,j)+aim(2,j)*fi1
	apc(2,j)=apc(2,j)-aim(2,j)
	apc(2,j)=apc(2,j)+aim(2,j)
	aim(2,j)=0
	fc(1,j)=fc(2,j)+fi1


	fil1=
     1   (fc(l3,j)-fc(l2,j))
     1   /(x(l3,j)-x(l2,j))*0.5*hksi(l2,j)
c	fil1=0
	con(l2,j)=con(l2,j)-aip(l2,j)*fil1
	apc(l2,j)=apc(l2,j)-aip(l2,j)
	apc(l2,j)=apc(l2,j)+aip(l2,j)
	aip(l2,j)=0
	fc(l1,j)=fc(l2,j)-fil1
10	continue
	enddo

	do i=2,l2
	goto 20
c	fj1=(fc(i,2)-fc(i,3))
c     1  /(y(i,2)-y(i,3))*0.5*heta(i,2)
	fj1=0

	con(i,2)=con(i,2)+ajm(i,2)*fj1
	apc(i,2)=apc(i,2)-ajm(i,2)
	apc(i,2)=apc(i,2)+ajm(i,2)
	ajm(i,2)=0
	fc(i,1)=fc(i,2)+fj1

	fjm1=
     1   (fc(i,m3)-fc(i,m2))
     1   /(y(i,m3)-y(i,m2))*0.5*heta(i,m2)
c	fjm1=0
	con(i,m2)=con(i,m2)-ajp(i,m2)*fjm1
	ap(i,m2)=ap(i,m2)-ajp(i,m2)
	ap(i,m2)=ap(i,m2)+ajp(i,m2)
	ajp(i,m2)=0
	fc(i,l1)=fc(i,l2)-fjm1
20	continue
	enddo


	xmiu=4*3.14159265*1d-7
	epsilon=8.8542d-12

	do 600 i=2,l2
	do 600 j=2,m2
	rmu=rmu_c(i,j)
	conduc0=conduc(i,j)
	if(id_ele(i,j).eq.11)conduc0=0
	apc(i,j)=
     1  -xmiu*epsilon*omega*omega*vol(i,j)
     1    +aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
     1   +1/r(i,j)/r(i,j)*vol(i,j)
     1  +(0d0,1d0)*rmu*xmiu*conduc0*omega*vol(i,j)
600	continue
c	write(*,*)(apc(2,j),j=2,m2)


      if(idesign.eq.5)then
      do 700 id=1,10
      xx0=coil_xx0+coil_dxx*(id-1)
      yy0=coil_yy0
      call add_con_circle(xx0,yy0)
700	continue

      endif
      if(idesign.eq.1.or.idesign.eq.2.or.idesign.eq.3.or.
     1  idesign.eq.4.or.idesign.eq.6)then
      do 800 id=1,5
      xx0=coil_xx0+coil_dxx*(id-1)
      yy0=coil_yy0
      call add_con(xx0,yy0)
800   continue
      endif

	return
	end

	subroutine add_con(x0,y0)
c for sterling
	include "sc.h"

c	write(*,*)'coil x',x0-0.75*0.0254,x0+0.75*0.0254,
c     1         x0-(0.75-0.048)*0.0254,x0+(0.75-0.048)*0.0254
c	write(*,*)'coil y ',y0-0.1875*0.0254,y0+0.1875*0.0254,
c     1         y0-(0.1875-0.048)*0.0254,y0+(0.1875-0.048)*0.0254

	id=0
      area=0
	do i=2,l2
	do j=2,m2
c	rr2=(x(i,j)-x0)**2+(y(i,j)-y0)**2
c	if(rr2.le.(0.25*0.0254)**2)then

	xxo=0.75*0.0254
	xxi=(0.75-0.048)*0.0254
	yyo=0.1875*0.0254
	yyi=(0.1875-0.048)*0.0254
	if(dabs(x(i,j)-x0).le.xxo.and.
     1     dabs(y(i,j)-y0).le.yyo.and.
     1     (dabs(x(i,j)-x0).ge.xxi.or.
     1     dabs(y(i,j)-y0).ge.yyi))then
	id=id+1

      dxdeta=xc(i,j+1)-xc(i,j)
      dydeta=yc(i,j+1)-yc(i,j)
      dxdksi=x(i,j)-x(i-1,j)
      dydksi=y(i,j)-y(i-1,j)
      if(i.eq.2.or.i.eq.l1)dxdksi=dxdksi*2.0
      if(i.eq.2.or.i.eq.l1)dydksi=dydksi*2.0
      xjacb=dxdksi*dydeta-dxdeta*dydksi

      area=area+xjacb
      if(id.gt.6000)then
      write(*,*)'need iid with elements of',id
      stop 'err in add_con, check iid'
      endif
	iid(id)=i
	jjd(id)=j
	endif
	enddo
	enddo
	write(*,*)'cells in coil',id

	do i=1,id
	i1=iid(i)
	j1=jjd(i)
c I/area*vol=I*2*pai*r
c area=dr*dz
c vol=2*pai*r*dr*dz
      con(i1,j1)=con(i1,j1)+xmiu*cden/area*vol(i1,j1)
	enddo
	return
	end

      subroutine add_con_circle(x0,y0)
c this subroutine is the same as in sc_advan.f
	include "sc.h"


	id=0
      area=0
	do i=2,l2
	do j=2,m2
	rr2=(x(i,j)-x0)**2+(y(i,j)-y0)**2
	if(rr2.ge.((0.25-0.048)*0.0254)**2
     1     .and.rr2.le.(0.25*0.0254)**2)then
	id=id+1

      dxdeta=xc(i,j+1)-xc(i,j)
      dydeta=yc(i,j+1)-yc(i,j)
      dxdksi=x(i,j)-x(i-1,j)
      dydksi=y(i,j)-y(i-1,j)
      if(i.eq.2.or.i.eq.l1)dxdksi=dxdksi*2.0
      if(i.eq.2.or.i.eq.l1)dydksi=dydksi*2.0
      xjacb=dxdksi*dydeta-dxdeta*dydksi

      area=area+xjacb


      if(id.gt.4000)stop 'err in add_con'
	iid(id)=i
	jjd(id)=j
	endif
	enddo
	enddo
	write(*,*)'cells in coil',id  

	do i=1,id
	i1=iid(i)
	j1=jjd(i)
      con(i1,j1)=con(i1,j1)+xmiu*cden/area*vol(i1,j1)
	enddo
	return
	end


	subroutine test_cal_elemag
	include 'sc.h'
	complex*16 tes

	call test_grid_mag
c for xc,yc,hksi,heta,ak1,ae1...
	call setup1
	call start

c for large grid, iteration has to be more

	do 100 loop=1,100
	write(*,*)'matr_ele'
	call test_matr_ele

c	call gamsor
c	call flux_s
c	call extra2

	do loop1=1,100
	call solve_tdma_c
	enddo

	fmax_r=0
	fmax_i=0
	fmax=0
	do i=1,l1
	do j=1,m1
	rsl=abs(dreal(fc(i,j))-dreal(fco(i,j)))
	if(rsl.ge.fmax_r)then
	fmax_r=rsl
	endif
	rsl=abs(dimag(fc(i,j))-dimag(fco(i,j)))
	if(rsl.ge.fmax_i)then
	fmax_i=rsl
	endif

	fco(i,j)=fc(i,j)
	fc0=abs(dreal(fc(i,j)))
	if(fmax.lt.fc0)fmax=fc0
	enddo
	enddo
	write(*,*)'fcmax=',fmax
	fsmall=1d-15
c blc is the err of the solved equation
c fmax_r and fmax_i are not good to determine convergence
	if(fmax_r.gt.fsmall.or.fmax_i.gt.fsmall)goto 100
	if(blc.gt.1d-15)goto 100
	goto 200
100	continue
200	write(*,*)fmax_r,fmax_i,loop
c	write(*,*)(a_mag(l1/2,j,1),j=1,m1)
c	write(*,*)(fc(l1/2,j),j=1,m1)
	write(*,*)omega
	fmax=0
	pow_c=0
	do i=1,l1
	do j=1,m1
	qth0=0.5*conduc(i,j)*omega*omega
     1   *(cdabs(fc(i,j)))**2
	qth(i,j)=qth0
c vol: dz*dr*2*pai*r
	pow_c=pow_c+qth0*vol(i,j)
	if(fmax.lt.qth0)fmax=qth0
	enddo
	enddo
	write(*,*)'qthmax=',fmax
	write(*,*)'power consuming=',pow_c
	call zero(u,l1,m1)
	call zero(v,l1,m1)
	do i=2,l2
        u(i,1)=2*(dreal(fc(i,2))-dreal(fc(i,1)))/(y(i,2)-y(i,1))
	do j=2,m2
	u(i,j)=
c     1   (fc(i,j+1)*r(i,j+1)-fc(i,j)*r(i,j))/(y(i,j+1)-y(i,j))
c     1  /((r(i,j)+r(i,j+1))/2)
     1   (dreal(fc(i,j+1))-dreal(fc(i,j-1)))/(y(i,j+1)-y(i,j-1))
     1  +dreal(fc(i,j))/r(i,j)
	v(i,j)=-(dreal(fc(i+1,j))-dreal(fc(i,j)))/(x(i+1,j)-x(i,j))
	tes=
     1  (fc(i+1,j)-fc(i,j))/(x(i+1,j)-x(i,j))*ak1(i+1,j)
     1 -(fc(i,j)-fc(i-1,j))/(x(i,j)-x(i-1,j))*ak1(i,j)
     1 +(fc(i,j+1)-fc(i,j))/(y(i,j+1)-y(i,j))*ae1(i,j+1)
     1 -(fc(i,j)-fc(i,j-1))/(y(i,j)-y(i,j-1))*ae1(i,j)
	write(*,*)i,j,tes,fc(i,j)
	enddo
	enddo
	write(*,*)'u=',(u(i,2),i=2,l2)
	i=interid(0d0,x,l1)
	do j=1,m2
	write(*,*)'r,z,u,du,real(fc)=',y(i,j),x(i,j),u(i,j),
     1    u(i,j+1)-u(i,j),dreal(fc(i,j))
	enddo
	open(11,file='mag.plt')
	write(11,82)
	write(11,83)l1,m1
	write(11,88)((y(i,j),i=1,l1),j=1,m1)
	write(11,88)((x(i,j),i=1,l1),j=1,m1)
	write(11,88)((dreal(fc(i,j)),i=1,l1),j=1,m1)
	write(11,88)((dimag(fc(i,j)),i=1,l1),j=1,m1)
	write(11,88)((qth(i,j),i=1,l1),j=1,m1)
	write(11,88)((conduc(i,j),i=1,l1),j=1,m1)
	write(11,88)((u(i,j),i=1,l1),j=1,m1)
	write(11,88)((v(i,j),i=1,l1),j=1,m1)
	close(11)
	call zero(u,l1,m1)
	call zero(v,l1,m1)
   82 format(1x,'variables ="x","y","mag1","mag2","qth","con","u","v"')
   83 format(1x,'zone t= "zone 1" , i=',i3,' , j=',i3,' , f= block')
   88 format(6(1pe12.5,1x))
	open(11,file='mag.dat')
	write(11,*)l1,m1
	write(11,88)((qth(i,j),i=1,l1),j=1,m1)
	close(11)
	stop
	end

	subroutine test_us_elemag
	include 'sc.h'

	complex*16 fi1,fil1,fj1,fjm1

	entry test_matr_ele
c set omega cden
	omega=0
	cden=10

	call reset1
	do 110 j=2,m2
	diff=1/(0.5*hksi(2,j))
	aim(2,j)=diff*ak1(2,j)
	diff=1/(0.5*hksi(l2,j))
	aip(l2,j)=diff*ak1(l1,j)
	do 110 i=2,l3
	diff=1/(0.5*hksi(i,j)+0.5*hksi(i+1,j))
	aip(i,j)=diff*ak1(i+1,j)
	aim(i+1,j)=aip(i,j)
110	continue
c	write(*,*)'ak1=',(ak1(i,2),i=2,l1)
c	write(*,*)'ak2=',(ak2(i,2),i=2,l1)
c	write(*,*)'ae1=',(ae1(2,j),j=2,m1)
c	write(*,*)'ae2=',(ae2(2,j),j=2,m1)
	do 120 i=2,l2
	diff=1/(0.5*heta(i,2))
	ajm(i,2)=diff*ae1(i,2)

	diff=1/(0.5*heta(i,m2))
	ajp(i,m2)=diff*ae1(i,m1)
	do 115 j=2,m3
	diff=1/(0.5*heta(i,j)+0.5*heta(i,j+1))
	ajp(i,j)=diff*ae1(i,j+1)
	ajm(i,j+1)=ajp(i,j)
115	continue

c	diff=1/(heta(i,2)+0.5*heta(i,3))
c	ajp(i,2)=diff*ae1(i,3)
c	ajm(i,2)=-ajp(i,2)

120	continue

c boundary
	do j=2,m2
	goto 10
	fi1=
     1   (fc(2,j)-fc(3,j))
     1   /(x(2,j)-x(3,j))*0.5*hksi(2,j)
c	fi1=0
	con(2,j)=con(2,j)+aim(2,j)*fi1
	apc(2,j)=apc(2,j)-aim(2,j)
	apc(2,j)=apc(2,j)+aim(2,j)
	aim(2,j)=0
	fc(1,j)=fc(2,j)+fi1


	fil1=
     1   (fc(l3,j)-fc(l2,j))
     1   /(x(l3,j)-x(l2,j))*0.5*hksi(l2,j)
c	fil1=0
	con(l2,j)=con(l2,j)-aip(l2,j)*fil1
	apc(l2,j)=apc(l2,j)-aip(l2,j)
	apc(l2,j)=apc(l2,j)+aip(l2,j)
	aip(l2,j)=0
	fc(l1,j)=fc(l2,j)-fil1
10	continue
	enddo

	do i=2,l2
	goto 20
c	fj1=(fc(i,2)-fc(i,3))
c     1  /(y(i,2)-y(i,3))*0.5*heta(i,2)
	fj1=0

	con(i,2)=con(i,2)+ajm(i,2)*fj1
	apc(i,2)=apc(i,2)-ajm(i,2)
	apc(i,2)=apc(i,2)+ajm(i,2)
	ajm(i,2)=0
	fc(i,1)=fc(i,2)+fj1

	fjm1=
     1   (fc(i,m3)-fc(i,m2))
     1   /(y(i,m3)-y(i,m2))*0.5*heta(i,m2)
c	fjm1=0
	con(i,m2)=con(i,m2)-ajp(i,m2)*fjm1
	ap(i,m2)=ap(i,m2)-ajp(i,m2)
	ap(i,m2)=ap(i,m2)+ajp(i,m2)
	ajp(i,m2)=0
	fc(i,l1)=fc(i,l2)-fjm1
20	continue
	enddo


	xmiu=4*3.14159265*1d-7
	epsilon=8.8542d-12
	omega=0

	do 600 i=2,l2
	do 600 j=2,m2
	rmu=rmu_c(i,j)
	apc(i,j)=
     1  -xmiu*epsilon*omega*omega*vol(i,j)
     1    +aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
     1   +1/r(i,j)/r(i,j)*vol(i,j)
     1  +(0d0,1d0)*rmu*xmiu*conduc(i,j)*omega*vol(i,j)
600	continue
c	write(*,*)(apc(2,j),j=2,m2)

c      coil_xx0=0.0
c      coil_yy0=0.02
c      coil_dxx=0.0

	do 700 id=1,1
      xx0=coil_xx0+coil_dxx*(id-1)
      yy0=coil_yy0
	call test_add_con(xx0,yy0)
700	continue

	return
	end

	subroutine test_add_con(x0,y0)
	include "sc.h"
	id=0
      area=0
	do i=2,l2
	do j=2,m2
	rr2=sqrt((x(i,j)-x0)**2+(y(i,j)-y0)**2)
	if(rr2.le.1d-5)then
c	if(dabs(x(i,j)-x0).le.1d-4.and.
c     1     dabs(y(i,j)-y0).le.1d-4)then
	id=id+1

      dxdeta=xc(i,j+1)-xc(i,j)
      dydeta=yc(i,j+1)-yc(i,j)
      dxdksi=x(i,j)-x(i-1,j)
      dydksi=y(i,j)-y(i-1,j)
      if(i.eq.2.or.i.eq.l1)dxdksi=dxdksi*2.0
      if(i.eq.2.or.i.eq.l1)dydksi=dydksi*2.0
      xjacb=dxdksi*dydeta-dxdeta*dydksi

      area=area+xjacb

	if(id.gt.1000)stop 'err in add_con'
	iid(id)=i
	jjd(id)=j
	endif
	enddo
	enddo
c	write(*,*)(x(i,m1/2),i=1,l1)
c	write(*,*)(y(l1/2,j),j=1,m1)
	write(*,*)'test id,x0,y0',id,x0,y0
	xmiu=4*3.14159265*1d-7
	do i=1,id
	i1=iid(i)
	j1=jjd(i)
c I/area*ak1=I*2*pai*r
c area=dr*dz
c ak1=2*pai*r*dr*dz
      con(i1,j1)=con(i1,j1)+xmiu*cden/area*vol(i1,j1)
	write(*,*)'id=,cden,r',id,cden,r(i1,j1),vol(i1,j1)/r(i1,j1)
	write(*,*)ak1(i1,j1)/r(i1,j1),ae1(i1,j1)/r(i1,j1)
	enddo
	return
	end

	subroutine test_grid_mag
	include 'sc.h'
c points  u0         du_max
c  83x 83  0.000309   0.0015
c  83x103  0.0003098  0.0016
c 103x103  0.0003087  0.0018

c	data xxd/-0.1,0.1/
c	data yyd/0,0.1/
c	data iddx/0,81/
c	data iddy/0,81/
c points in z direction    B at r=0 and z=0
c    43                3.257e-4 T
c    83                3.184e-4
c   103                3.176e-4

c	data xxd/-0.1,0.1/
c	data yyd/0,0.0195,0.0205,0.1/
c	data iddx/0,81/
c	data iddy/0,20,11,20/

c points in z direction    B at r=0 and z=0
c         93x63           3.145e-4 T
c        103x63           3.143e-4

c	data xxd/-0.1,-0.0005,0.0005,0.1/
c	data yyd/0,0.0195,0.0205,0.1/
c	data iddx/0,40,21,40/
c	data iddy/0,40,21,40/
c B0=3.241e-4

c	data xxd/-0.1,-0.01,-0.00005,0.00005,0.01,0.1/
c	data yyd/0,0.01,0.01995,0.02005,0.03,0.1/
c	data iddx/0,20,10,21,10,20/
c	data iddy/0,20,10,21,10,20/

c	data xxd/-0.1,-0.01,-0.001,-0.00005,0.00005,0.001,0.01,0.1/
c	data yyd/0,0.01,0.019,0.01995,0.02005,0.021,0.03,0.1/
c	data iddx/0,20,10,10,21,10,10,20/
c	data iddy/0,20,10,10,21,10,10,20/


c   B0=0.0003176   fcm=1.094d-5
c        parameter (nxx=8,nyy=8)
c        dimension xxd(nxx),yyd(nyy),iddx(nxx),iddy(nyy)
c        data xxd/-1,-0.1,-0.02,-0.005,0.005,0.02,0.1,1/
c        data yyd/0,0.015,0.019,0.0195,0.0205,0.025,0.1,1/
c        data iddx/0,5,20,20,21,20,20,5/
c        data iddy/0,10,10,10,31,20,20,5/


	xxl=xxd(nxx)-xxd(1)
	l1=2
	do i=2,nxx
	l1=l1+iddx(i)
	enddo
	yyl=yyd(nyy)-yyd(1)
	m1=2
	do i=2,nyy
	m1=m1+iddy(i)
	enddo

	if(l1.gt.ni-1.or.m1.gt.nj-1)then
	stop 'err in test_grid_mag'
	endif

	i=2
	do j=2,m1
	xc(2,j)=xxd(1)
	enddo

	do id=2,nxx
c	idd=xxd(id)/xxl*(l1-2)+0.5
	idd=iddx(id)
	if(idd.lt.1)idd=1
	write(*,*)idd
	if(id.eq.nxx)idd=l1-i
	do idx=1,idd
	do j=2,m1
	xc(i+1,j)=xxd(id-1)+(xxd(id)-xxd(id-1))*idx/idd
	enddo
	i=i+1
	enddo
	enddo
	write(*,*)l1,'=?',i,xxl,'=?',xc(l1,m1/2)-xc(2,m1/2)
c	write(*,*)(xc(i,m1/2),i=2,l1)

	j=2
	do i=2,l1
	yc(i,2)=yyd(1)
	enddo

	do id=2,nyy
c	idd=yyd(id)/yyl*(m1-2)+0.5
	idd=iddy(id)
	write(*,*)idd
	if(idd.lt.1)idd=1
	if(id.eq.nyy)idd=m1-j
	do idx=1,idd
	do i=2,l1
	yc(i,j+1)=yyd(id-1)+(yyd(id)-yyd(id-1))*idx/idd
	enddo
	j=j+1
	enddo
	enddo
	write(*,*)m1,'=?',j,yyl,'=?',yc(l1/2,m1)-yc(l1/2,2)
c	write(*,*)(yc(l1/2,j),j=2,m1)
	return
	end

	subroutine grid_mag_design3
	include 'sc.h'


	xxl=xxd(nxx)-xxd(1)
	l1=2
	do i=2,nxx
	if(xxd(i).le.xxd(i-1))stop 'err in grid:x(i)>x(i+1)'
	l1=l1+iddx(i)
	enddo
	yyl=yyd(nyy)-yyd(1)
	m1=2
	do i=2,nyy
	if(yyd(i).le.yyd(i-1))stop 'err in grid:y(i)>y(i+1)'
	m1=m1+iddy(i)
	enddo

	i=2
	do j=2,m1
	xc(2,j)=xxd(1)
	enddo

	do id=2,nxx
c	idd=xxd(id)/xxl*(l1-2)+0.5
	idd=iddx(id)
	if(idd.lt.1)idd=1
	write(*,*)idd
	if(id.eq.nxx)idd=l1-i
	do idx=1,idd
	do j=2,m1
	xc(i+1,j)=xxd(id-1)+(xxd(id)-xxd(id-1))*idx/idd
	enddo
	i=i+1
	enddo
	enddo
	write(*,*)l1,'=?',i,xxl,'=?',xc(l1,m1/2)-xc(2,m1/2)
c	write(*,*)(xc(i,m1/2),i=2,l1)

	j=2
	do i=2,l1
	yc(i,2)=yyd(1)
	enddo

	do id=2,nyy
c	idd=yyd(id)/yyl*(m1-2)+0.5
	idd=iddy(id)
	write(*,*)idd
	if(idd.lt.1)idd=1
	if(id.eq.nyy)idd=m1-j
	do idx=1,idd
	do i=2,l1
	yc(i,j+1)=yyd(id-1)+(yyd(id)-yyd(id-1))*idx/idd
	enddo
	j=j+1
	enddo
	enddo
	write(*,*)m1,'=?',j,yyl,'=?',yc(l1/2,m1)-yc(l1/2,2)
c	write(*,*)(yc(l1/2,j),j=2,m1)
	return
	end

	subroutine grid_mag_design2
	include 'sc.h'

c	parameter (nxx=33,nyy=24)
c	data xxd/-1.2,-0.080,-0.055,-0.03,-0.007,0,
c     1  0.005,0.01,0.02,0.03,0.04,0.046,0.050,0.058,
c     1  0.07,0.08,0.10,0.11,0.122,
c     1   0.130,0.140,0.178,0.190,
c     1   0.193,0.195,0.200,0.207,
c     1   0.213,0.235,0.255,0.270,0.280,1.2/
c	data yyd/0,0.008,0.020,0.023,0.0385,0.0405,0.0435,0.051,0.055,
c     1   0.060,0.080,0.085,0.105,
c     1    0.1225,0.125,0.19,
c     1  0.196,0.197219,0.204306,0.205525,0.213,0.28,0.4,1.2/
c	data iddx/0,8,2,5,8,5,
c     1   6,3,7,5,4,6,3,3,
c     1  4,5,10,10,6,
c     1   16,8,10,8,
c     1    6,3,7,7,
c     1    4,10,4,3,2,8/
c	data iddy/0,4,6,4,6,2,2,3,3,
c     1   10,16,3,5,
c     1   8,2,8,4,
c     1   15,4,15,5,5,3,12/

	xxl=xxd(nxx)-xxd(1)
	l1=2
	do i=2,nxx
c      iddx(i)=iddx(i)*2
	if(xxd(i).le.xxd(i-1))stop 'err in grid:x(i)>x(i+1)'
	l1=l1+iddx(i)
	enddo
	yyl=yyd(nyy)-yyd(1)
	m1=2
	do i=2,nyy
c      iddy(i)=iddy(i)*2
	if(yyd(i).le.yyd(i-1))stop 'err in grid:y(i)>y(i+1)'
	m1=m1+iddy(i)
	enddo

	i=2
	do j=2,m1
	xc(2,j)=xxd(1)
	enddo

	do id=2,nxx
c	idd=xxd(id)/xxl*(l1-2)+0.5
	idd=iddx(id)
	if(idd.lt.1)idd=1
	write(*,*)idd
	if(id.eq.nxx)idd=l1-i
	do idx=1,idd
	do j=2,m1
	xc(i+1,j)=xxd(id-1)+(xxd(id)-xxd(id-1))*idx/idd
	enddo
	i=i+1
	enddo
	enddo
	write(*,*)l1,'=?',i,xxl,'=?',xc(l1,m1/2)-xc(2,m1/2)
c	write(*,*)(xc(i,m1/2),i=2,l1)

	j=2
	do i=2,l1
	yc(i,2)=yyd(1)
	enddo

	do id=2,nyy
c	idd=yyd(id)/yyl*(m1-2)+0.5
	idd=iddy(id)
	write(*,*)idd
	if(idd.lt.1)idd=1
	if(id.eq.nyy)idd=m1-j
	do idx=1,idd
	do i=2,l1
	yc(i,j+1)=yyd(id-1)+(yyd(id)-yyd(id-1))*idx/idd
	enddo
	j=j+1
	enddo
	enddo
	write(*,*)m1,'=?',j,yyl,'=?',yc(l1/2,m1)-yc(l1/2,2)
c	write(*,*)(yc(l1/2,j),j=2,m1)
	return
	end

        subroutine grid_mag_design4
	include 'sc.h'
c	parameter (nxx=33,nyy=24)

c	data xxd/-1.2,-0.080,-0.055,-0.03,-0.007,0,
c     1  0.005,0.01,0.02,0.03,0.04,0.046,0.050,0.058,
c     1  0.07,0.08,0.10,0.11,0.122,
c     1   0.130,0.140,0.178,0.190,
c     1   0.193,0.195,0.200,0.207,
c     1   0.213,0.235,0.255,0.270,0.280,1.2/
c	data yyd/0,0.008,0.020,0.023,0.0385,0.0405,0.0435,0.051,0.055,
c     1   0.060,0.080,0.085,0.105,
c     1    0.1225,0.125,0.19,
c     1  0.196,0.197219,0.204306,0.205525,0.213,0.28,0.4,1.2/
c	data iddx/0,8,2,5,8,5,
c     1   6,3,7,5,4,6,3,3,
c     1  4,5,10,10,6,
c     1   16,8,10,8,
c     1    6,3,7,7,
c     1    4,10,4,3,2,8/
c	data iddy/0,4,6,4,6,2,2,3,3,
c     1   10,16,3,5,
c     1   8,2,8,4,
c     1   15,4,15,5,5,3,12/

	xxl=xxd(nxx)-xxd(1)
	l1=2
	do i=2,nxx
	if(xxd(i).le.xxd(i-1))stop 'err in grid:x(i)>x(i+1)'
	l1=l1+iddx(i)
	enddo
	yyl=yyd(nyy)-yyd(1)
	m1=2
	do i=2,nyy
	if(yyd(i).le.yyd(i-1))stop 'err in grid:y(i)>y(i+1)'
	m1=m1+iddy(i)
	enddo

	i=2
	do j=2,m1
	xc(2,j)=xxd(1)
	enddo

	do id=2,nxx
c	idd=xxd(id)/xxl*(l1-2)+0.5
	idd=iddx(id)
	if(idd.lt.1)idd=1
	write(*,*)idd
	if(id.eq.nxx)idd=l1-i
	do idx=1,idd
	do j=2,m1
	xc(i+1,j)=xxd(id-1)+(xxd(id)-xxd(id-1))*idx/idd
	enddo
	i=i+1
	enddo
	enddo
	write(*,*)l1,'=?',i,xxl,'=?',xc(l1,m1/2)-xc(2,m1/2)
c	write(*,*)(xc(i,m1/2),i=2,l1)

	j=2
	do i=2,l1
	yc(i,2)=yyd(1)
	enddo

	do id=2,nyy
c	idd=yyd(id)/yyl*(m1-2)+0.5
	idd=iddy(id)
	write(*,*)idd
	if(idd.lt.1)idd=1
	if(id.eq.nyy)idd=m1-j
	do idx=1,idd
	do i=2,l1
	yc(i,j+1)=yyd(id-1)+(yyd(id)-yyd(id-1))*idx/idd
	enddo
	j=j+1
	enddo
	enddo
	write(*,*)m1,'=?',j,yyl,'=?',yc(l1/2,m1)-yc(l1/2,2)
c	write(*,*)(yc(l1/2,j),j=2,m1)
	return
	end

        subroutine grid_mag_design5
	include 'sc.h'
c	parameter (nxx=33,nyy=24)

c	data xxd/-1.2,-0.080,-0.055,-0.03,-0.007,0,
c     1  0.005,0.01,0.02,0.03,0.04,0.046,0.050,0.058,
c     1  0.07,0.08,0.10,0.11,0.122,
c     1   0.130,0.140,0.178,0.190,
c     1   0.193,0.195,0.200,0.207,
c     1   0.213,0.235,0.255,0.270,0.280,1.2/
c      data yyd/0,0.008,0.020,0.023,0.0385,0.0405,0.0435,0.064,0.068,
c     1   0.075,0.100,0.125,0.145,
c     1    0.1625,0.165,0.23,
c     1  0.233,0.235,0.256,0.258,0.26,0.3,0.4,1.2/
c	data iddx/0,8,2,5,8,5,
c     1   6,3,7,5,4,6,3,3,
c     1  4,5,10,10,6,
c     1   16,8,10,8,
c     1    6,3,7,7,
c     1    4,10,4,3,2,8/
c	data iddy/0,4,6,4,6,2,2,3,3,
c     1   10,16,3,5,
c     1   8,2,8,4,
c     1   15,4,15,5,5,3,12/

	xxl=xxd(nxx)-xxd(1)
	l1=2
	do i=2,nxx
	if(xxd(i).le.xxd(i-1))stop 'err in grid:x(i)>x(i+1)'
	l1=l1+iddx(i)
	enddo
	yyl=yyd(nyy)-yyd(1)
	m1=2
	do i=2,nyy
	if(yyd(i).le.yyd(i-1))stop 'err in grid:y(i)>y(i+1)'
	m1=m1+iddy(i)
	enddo

	i=2
	do j=2,m1
	xc(2,j)=xxd(1)
	enddo

	do id=2,nxx
	idd=iddx(id)
	if(idd.lt.1)idd=1
	write(*,*)idd
	if(id.eq.nxx)idd=l1-i
	do idx=1,idd
	do j=2,m1
	xc(i+1,j)=xxd(id-1)+(xxd(id)-xxd(id-1))*idx/idd
	enddo
	i=i+1
	enddo
	enddo
	write(*,*)l1,'=?',i,xxl,'=?',xc(l1,m1/2)-xc(2,m1/2)
c	write(*,*)(xc(i,m1/2),i=2,l1)

	j=2
	do i=2,l1
	yc(i,2)=yyd(1)
	enddo

	do id=2,nyy
	idd=iddy(id)
	write(*,*)idd
	if(idd.lt.1)idd=1
	if(id.eq.nyy)idd=m1-j
	do idx=1,idd
	do i=2,l1
	yc(i,j+1)=yyd(id-1)+(yyd(id)-yyd(id-1))*idx/idd
	enddo
	j=j+1
	enddo
	enddo
	write(*,*)m1,'=?',j,yyl,'=?',yc(l1/2,m1)-yc(l1/2,2)
c	write(*,*)(yc(l1/2,j),j=2,m1)
	return
	end




      subroutine imgrid
	include 'sc.h'
	dimension xb(ni),zb1(ni),high(ni),flow1(ni),pres(ni)


      entry initgrid
      ll1=lst3/2+1
      do 22 i=2,lst3
      do 22 j=2,m1
      power=1.5
      if(i.le.ll1)t1=0.5*(dfloat(i-2)/dfloat(ll1-2))**power
      if(i.ge.ll1)t1=1.0-0.5*(dfloat(lst3-i)/dfloat(lst3-ll1))**power
      xc(i,j)=t1*xlen*0.95
   22 continue
      do 23 j=2,m1
      do 23 i=lst3,l1
      power=1.0
      t1=1.0*(dfloat(i-lst3)/dfloat(l1-lst3))**power
      xc(i,j)=xc(lst3,j)+t1*xlen*0.05
   23 continue
      llst1=m1/2+1
      do 24 i=2,l1
      do 24 j=2,m1
      power=1.5
      if(j.le.llst1)t1=0.5*(dfloat(j-2)/dfloat(llst1-2))**power
      if(j.ge.llst1)
     1  t1=1.0-0.5*(dfloat(m1-j)/dfloat(m1-llst1))**power
      yc(i,j)=t1*ylen
   24 continue
      return

c --------------------------------------------------------------------

      entry newstep

      iter=0
      l2=l1-1
      m2=m1-1
      do 80 i=1,l1
      do 80 j=1,m1
      fo(i,j,1)=f(i,j,1)
      fo(i,j,2)=f(i,j,2)
      fo(i,j,5)=f(i,j,5)
      fo(i,j,6)=f(i,j,6)
	rho_o(i,j)=rho(i,j)
   80 continue
      return

c --------------------------------------------------------------------

      entry movegrid

      do 110 j=2,m1
  110 zb1(j)=xc(lst3,j)
c
      do 201 j=2,m2
  201 b1jbl(j)=ak1(lst3,j)/heta(lst3,j)*(t(lst3-1,j)-t(lst3,j))
     1 /(0.5*hksi(lst3-1,j)+0.5*hksi(lst3,j))
      b1jbl(1)=b1jbl(2)
      b1jbl(m1)=b1jbl(m2)
c
      do 202 j=2,m2
  202 b2jbl(j)=ak1(lst3,j)/heta(lst3,j)*(t(lst3,j)-t(lst3+1,j))
     1 /(0.5*hksi(lst3,j)+0.5*hksi(lst3+1,j))
      b2jbl(1)=b2jbl(2)
      b2jbl(m1)=b2jbl(m2)
c
      do 203 j=2,m2
      dxm1=xc(lst3,j+1)-xc(lst3,j)
      dym1=yc(lst3,j+1)-yc(lst3,j)
      dlen=sqrt(dxm1*dxm1+dym1*dym1)
      as=dxm1/dlen
      ac=dym1/dlen
      dx1=dtm*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel*ac
      dy1=dtm*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel*as
      flow1(j)=dx1
  203 continue
      flowmax=0.0
      flowmin=0.0
      do 204 j=2,m2
      flowmax=dmax1(flowmax,flow1(j))
      flowmin=dmin1(flowmin,flow1(j))
  204 continue
      do 207 j=2,m2
      xb(j)=xc(lst3,j)+flow1(j)
  207 xc(lst3,j)=xb(j)
      xc(lst3,2)=xc(lst3,3)
      xc(lst3,m1)=xc(lst3,m2)
c
      print 227,flowmin,flowmax,dtm
      print 228,flow1(21),flow1(19),flow1(17),flow1(15),flow1(13)
     1  ,flow1(11),flow1(10)
      print 228,flow1(9),flow1(8),flow1(7),flow1(6),flow1(5),flow1(4),
     1   flow1(3)
	print *,'     solid-liquid interface movement'
  227 format(/2x,'movemin = ',f9.4,'movemax = ',f9.4,'dtm =',e9.4/)
  228 format(2x,7f9.4)
c
      do 261 j=2,m1
      coeff=(xc(2,j)-xc(lst3,j))/(xc(2,j)-zb1(j))
      coff=(xlen-xc(lst3,j))/(xlen-zb1(j))
      do 262 i=2,lst3-1
  262 xc(i,j)=xc(2,j)-(xc(2,j)-xc(i,j))*coeff
      do 263 i=lst3+1,l1
      xc(i,j)=xlen-(xlen-xc(i,j))*coff
  263 continue
  261 continue
      call adaptive
c
      if(rdtm.eq.0.0) return
      if(rdtm.ne.0.0) return
c
      do 270 i=2,l2
      do 270 j=2,m2
      delx=0.25*(xc(i,j)+xc(i+1,j)+xc(i,j+1)+xc(i+1,j+1))-x(i,j)
      dely=0.25*(yc(i,j)+yc(i+1,j)+yc(i,j+1)+yc(i+1,j+1))-y(i,j)
      dxdksi=x(i+1,j)-x(i-1,j)
      dydksi=y(i+1,j)-y(i-1,j)
      dxdeta=x(i,j+1)-x(i,j-1)
      dydeta=y(i,j+1)-y(i,j-1)
      xyjac=dxdksi*dydeta-dydksi*dxdeta
      if(abs(xyjac).le.1.e-9) go to 272
      do 280 k=1,nfmax
      if(k.eq.3) go to 281
      if(k.eq.4) go to 281
      dfdksi=f(i+1,j,k)-f(i-1,j,k)
      dfdeta=f(i,j+1,k)-f(i,j-1,k)
      fo(i,j,k)=fo(i,j,k)+((dfdksi*dydeta-dfdeta*dydksi)*delx
     1     -(dfdksi*dxdeta-dfdeta*dxdksi)*dely)/xyjac
  281 continue
  280 continue
  272 continue
  270 continue
c
      return
      end

	function diff_i(i,j)
	include "sc.h"
	if(gam(i,j).eq.0.or.gam(i+1,j).eq.0)then
	diff_i=0
	return
	endif
	diff_i=gam(i,j)*gam(i+1,j)/(0.5*hksi(i,j)*gam(i+1,j)
     1    +0.5*hksi(i+1,j)*gam(i,j))
	return
	end

	function diff_j(i,j)
	include "sc.h"
	if(gam(i,j).eq.0.or.gam(i,j+1).eq.0)then
	diff_j=0
	return
	endif

	diff_j=gam(i,j)*gam(i,j+1)/(0.5*heta(i,j)*gam(i,j+1)
     1    +0.5*heta(i,j+1)*gam(i,j))
	return
	end

      subroutine setup
	include "sc.h"

      data fl00/0./

c -------------------------------------------------------------------

      entry initl

	write(*,*)'initl'
	write(*,*)'rho(2,2)',rho(2,2)
c-----find rhou and pressure at the interface for the initial guess---
c --- calculate p_xi, rho*u_xi, p_eta, rho*u_eta for inital guess
c --- boundary condition is given by u and v. then change to xi, eta 
c --- direction here, then remain unchange in the program. 

      if(iter.ne.0) return
      do 50 i=2,l1
      do 50 j=1,m1
      ucs=(u(i-1,j)*hksi(i,j)+u(i,j)*hksi(i-1,j))/(hksi(i,j)
     1        +hksi(i-1,j))
      vcs=(v(i-1,j)*hksi(i,j)+v(i,j)*hksi(i-1,j))/(hksi(i,j)
     1        +hksi(i-1,j))
      pksi(i,j)=(p(i-1,j)*hksi(i,j)+p(i,j)*hksi(i-1,j))/(hksi(i,j)
     1          +hksi(i-1,j))
      xcomp=x(i,j)-x(i-1,j)
      ycomp=y(i,j)-y(i-1,j)
      if(abs(xcomp).lt.1.e-6) xcomp=0.0
      if(abs(ycomp).lt.1.e-6) ycomp=0.0
      ulen=sqrt(xcomp**2+ycomp**2)
      uksi=(ucs*xcomp+vcs*ycomp)/ulen
      rhov=(rho(i-1,j)*hksi(i,j)+rho(i,j)*hksi(i-1,j))/(hksi(i,j)
     1          +hksi(i-1,j))
      ruksi(i,j)=rhov*uksi
   50 continue

      do 60 i=1,l1
      do 60 j=2,m1
      ucs=(u(i,j-1)*heta(i,j)+u(i,j)*heta(i,j-1))/(heta(i,j)
     1        +heta(i,j-1))
      vcs=(v(i,j-1)*heta(i,j)+v(i,j)*heta(i,j-1))/(heta(i,j)
     1        +heta(i,j-1))
      peta(i,j)=(p(i,j-1)*heta(i,j)+p(i,j)*heta(i,j-1))/(heta(i,j)
     1          +heta(i,j-1))
      xcomp=x(i,j)-x(i,j-1)
      ycomp=y(i,j)-y(i,j-1)
      if(abs(xcomp).lt.1.e-6) xcomp=0.0
      if(abs(ycomp).lt.1.e-6) ycomp=0.0
      ulen=sqrt(xcomp**2+ycomp**2)
      ueta=(ucs*xcomp+vcs*ycomp)/ulen
      rhov=(rho(i,j-1)*heta(i,j)+rho(i,j)*heta(i,j-1))/(heta(i,j)
     1          +heta(i,j-1))
      rueta(i,j)=rhov*ueta
   60 continue
      return


      entry setup1

      l2=l1-1
      l3=l2-1
      m2=m1-1
      m3=m2-1

c-----calculate and output the position of the primary nodes
c-----and boundary nodes------
c --- calculate x(i,j), y(i,j) from inital xc(i,j), yc(i,j)----

      do 10 i=2,l2
      do 10 j=2,m2
      x(i,j)=0.25*(xc(i,j)+xc(i+1,j)+xc(i,j+1)+xc(i+1,j+1))
      y(i,j)=0.25*(yc(i,j)+yc(i+1,j)+yc(i,j+1)+yc(i+1,j+1))
   10 continue
      do 11 i=2,l2
      x(i,1)=0.5*(xc(i,2)+xc(i+1,2))
      y(i,1)=0.5*(yc(i,2)+yc(i+1,2))
      x(i,m1)=0.5*(xc(i,m1)+xc(i+1,m1))
      y(i,m1)=0.5*(yc(i,m1)+yc(i+1,m1))
   11 continue
      do 12 j=2,m2
      x(1,j)=0.5*(xc(2,j)+xc(2,j+1))
      y(1,j)=0.5*(yc(2,j)+yc(2,j+1))
      x(l1,j)=0.5*(xc(l1,j)+xc(l1,j+1))
      y(l1,j)=0.5*(yc(l1,j)+yc(l1,j+1))
   12 continue
      x(1,1)=xc(2,2)
      y(1,1)=yc(2,2)
      x(l1,m1)=xc(l1,m1)
      y(l1,m1)=yc(l1,m1)
      x(1,m1)=xc(2,m1)
      y(1,m1)=yc(2,m1)
      x(l1,1)=xc(l1,2)
      y(l1,1)=yc(l1,2)
c
c-----calculate local scale factors and volumes-----
c --- calculate h_xi, h_eta, Ja. -----------------
c
      do 20 i=2,l2
      do 20 j=2,m2
      xp=0.5*(xc(i+1,j+1)+xc(i+1,j))
      xm=0.5*(xc(i,j+1)+xc(i,j))
      yp=0.5*(yc(i+1,j+1)+yc(i+1,j))
      ym=0.5*(yc(i,j+1)+yc(i,j))
      xu=0.5*(xc(i,j+1)+xc(i+1,j+1))
      xd=0.5*(xc(i,j)+xc(i+1,j))
      yu=0.5*(yc(i,j+1)+yc(i+1,j+1))
      yd=0.5*(yc(i,j)+yc(i+1,j))
      dxdksi=xp-xm
      dxdeta=xu-xd
      dydksi=yp-ym
      dydeta=yu-yd
      hksi(i,j)=sqrt(dxdksi**2+dydksi**2)
      heta(i,j)=sqrt(dxdeta**2+dydeta**2)
      vol(i,j)=dxdksi*dydeta-dxdeta*dydksi
   20 continue
      do 21 i=1,l1
      heta(i,1)=0.0
      heta(i,m1)=0.0
      vol(i,1)=0.0
      vol(i,m1)=0.0
      if(i.eq.1.or.i.eq.l1) go to 22
      hksi(i,1)=sqrt((xc(i+1,2)-xc(i,2))**2+(yc(i+1,2)-yc(i,2))**2)
      hksi(i,m1)=
     1   sqrt((xc(i+1,m1)-xc(i,m1))**2+(yc(i+1,m1)-yc(i,m1))**2)
      go to 21
   22 continue
      hksi(i,1)=0.0
      hksi(i,m1)=0.0
   21 continue
      do 25 j=1,m1
      hksi(1,j)=0.0
      hksi(l1,j)=0.0
      vol(1,j)=0.0
      vol(l1,j)=0.0
      if(j.eq.1.or.j.eq.m1) go to 26
      heta(1,j)=sqrt((xc(2,j)-xc(2,j+1))**2+(yc(2,j)-yc(2,j+1))**2)
      heta(l1,j)
     1   =sqrt((xc(l1,j)-xc(l1,j+1))**2+(yc(l1,j)-yc(l1,j+1))**2)
      go to 25
   26 continue
      heta(1,j)=0.0
      heta(l1,j)=0.0
   25 continue
c
c-----calculate areas on the control-volume faces-----
c --- calculate alpha_xi, beta_xi, alpha_eta, beta_eta ----
c
      do 30 j=2,m2
      do 30 i=2,l1
      dxdeta=xc(i,j+1)-xc(i,j)
      dydeta=yc(i,j+1)-yc(i,j)
      dxdksi=x(i,j)-x(i-1,j)
      dydksi=y(i,j)-y(i-1,j)
      if(i.eq.2.or.i.eq.l1)dxdksi=dxdksi*2.0
      if(i.eq.2.or.i.eq.l1)dydksi=dydksi*2.0
      t1=dxdeta**2+dydeta**2
      t2=dxdksi**2+dydksi**2
      xjacb=dxdksi*dydeta-dxdeta*dydksi
      t3=dxdksi*dxdeta+dydeta*dydksi
      ak1(i,j)=sqrt(t2)*t1/xjacb
      ak2(i,j)=sqrt(t1)*t3/xjacb
   30 continue
c
      do 40 i=2,l2
      do 40 j=2,m1
      dxdeta=x(i,j)-x(i,j-1)
      dydeta=y(i,j)-y(i,j-1)
      dxdksi=xc(i+1,j)-xc(i,j)
      dydksi=yc(i+1,j)-yc(i,j)
      if(j.eq.2.or.j.eq.m1)dxdeta=dxdeta*2.0
      if(j.eq.2.or.j.eq.m1)dydeta=dydeta*2.0
      t1=dxdeta**2+dydeta**2
      t2=dxdksi**2+dydksi**2
      xjacb=dxdksi*dydeta-dxdeta*dydksi
      t3=dxdksi*dxdeta+dydeta*dydksi
      ae1(i,j)=sqrt(t1)*t2/xjacb
      ae2(i,j)=sqrt(t2)*t3/xjacb
   40 continue

      do 41 i=1,l1
      do 41 j=1,m1
      if(mode.eq.0) r(i,j)=1.0
      if(mode.eq.1) r(i,j)=y(i,j)+1d-8
	if(i.ge.2)then
      ak1(i,j)=ak1(i,j)*r(i,j)*2*pai
      ak2(i,j)=ak2(i,j)*r(i,j)*2*pai
	endif
	if(j.ge.2)then
      ae1(i,j)=ae1(i,j)*yc(i,j)*2*pai
      ae2(i,j)=ae2(i,j)*yc(i,j)*2*pai
	endif
      vol(i,j)=vol(i,j)*r(i,j)*2*pai
   41 continue

      return
	end

	subroutine setup2
	include 'sc.h'

c-----coefficients for all equations------
	write(*,*)'v=',v(i_rb,j_rr-1),v(i_rb,j_rr),v(i_rb,j_rr+1)
      do 100 nvar=1,nfmax
c	write(*,*)'nf,rueta(12,69)',nf,rueta(12,69)
	write(*,*)'var=',nvar
      nf=nvar







      if(.not.lsolve(nf))go to 100
	if(nf.le.4)then
	do i=1,4
	if(.not.lcon(i))goto 101
	enddo
	goto 100
	endif

101	if(nf.ge.5.and.lcon(nf))goto 100
      if(nf.eq.3) go to 400
      if(nf.eq.4) go to 500
c for 1:u 2:v 5:t
c reset ap,con
	call reset
	call gamsor

	if(nf.eq.6)then
	write(*,*)'ccen',ap(i_rb,2),aip(i_rb,2),aim(i_rb,2),
     1     ajp(i_rb,2),ajm(i_rb,2),con(i_rb,2),ccen(i_rb,2),
     1   fo(i_rb,2,6)
	endif
      rel=1.0-relax(nf)

c --- (3.14) find a_e, a_w, a_n, a_s --- 

      do 110 j=2,m2
      diff=gam(1,j)/(0.5*hksi(2,j))
c     if(nf.eq.1.or.nf.eq.2) diff=gam(2,j)/(0.5*hksi(2,j))
      flow=-ruksi(2,j)
	vv=void(2,j)
      call diflow(vv,2,j)
      aim(2,j)=acof*ak1(2,j)
  111 continue
      diff=gam(l1,j)/(0.5*hksi(l2,j))
c     if(nf.eq.1.or.nf.eq.2) diff=gam(l2,j)/(0.5*hksi(l2,j))
      flow=ruksi(l1,j)
	vv=void(l1,j)
      call diflow(vv,l1,j)
      aip(l2,j)=acof*ak1(l1,j)
  112 continue
      do 110 i=2,l3

      diff=diff_i(i,j)

      flow=ruksi(i+1,j)
	vv=dmax1(void(i,j),void(i+1,j))
      call diflow(vv,i+1,j)
      aip(i,j)=acof*ak1(i+1,j)
      aim(i+1,j)=aip(i,j)+flow*ak1(i+1,j)
  110 continue

      do 120 i=2,l2
      diff=gam(i,2)/(0.5*heta(i,2))
c     if(nf.eq.1.or.nf.eq.2) diff=gam(i,2)/(0.5*heta(i,2))
      flow=-rueta(i,2)
	vv=void(i,2)
      call diflow(vv,i,2)
      ajm(i,2)=acof*ae1(i,2)
  121 continue
      diff=gam(i,m1)/(0.5*heta(i,m2))
c     if(nf.eq.1.or.nf.eq.2) diff=gam(i,m2)/(0.5*heta(i,m2))
      flow=rueta(i,m1)
	vv=void(i,m1)
      call diflow(vv,i,m1)
      ajp(i,m2)=acof*ae1(i,m1)
  122 continue
      do 120 j=2,m3
      diff=diff_j(i,j)
      flow=rueta(i,j+1)
	vv=dmax1(void(i,j),void(i,j+1))
      call diflow(vv,i,j+1)
      ajp(i,j)=acof*ae1(i,j+1)
      ajm(i,j+1)=ajp(i,j)+flow*ae1(i,j+1)
	if(i.eq.2.and.j.eq.2)then
	write(*,*)'ajp',diff,flow,ajp(2,2),gam(2,2),gam(2,3),gam(2,1)
	endif
  120 continue

c	write(*,*)'ruksi',ruksi(i_rb-1,1),ruksi(i_rb-1,2),
c     1    ruksi(i_rb,1),ruksi(i_rb,2)
	write(*,*)i_rb,i_rt
	do i=i_rb,i_rt
	write(*,*)"ruksi=",ruksi(i,2)
	enddo
	if(nf.eq.1)then
      call awrite(1)
	endif
	if(nf.eq.2)then
      call awrite(2)
	endif
	write(*,*)'ajp',ajp(2,2)

c --- (3-14) finish ---
c --- (3-25) momentum interpolation ---

      if(nf.ne.1) go to 200 
c for 1:u

	t_ref=t(i_rt,2)
      do 130 i=2,l2
      do 130 j=2,m2
      apt=rdtm*rho(i,j)*vol(i,j)/void(i,j)

	sour=0
c density of argon
     1     +rho(i,j)*9.81*(1-t_ref/t(i,j))*vol(i,j)
      con(i,j)=con(i,j)+fo(i,j,1)*apt
     1    +sour
c     ap(i,j)=apt+ex1(i,j)/void(i,j)/void(i,j)
      ap(i,j)=ap(i,j)+apt
      if(id_ele(i,j).eq.3)then
	ap(i,j)=ap(i,j)
c +1/darcy*vol(i,j)
	endif
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)

130	continue
	call bound
	write(*,*)'sour=',con(i0-4,2),i0


	write(*,*)'sour=',con(124,2)

      call flux(0)
      do 215 i=2,l2
      do 215 j=2,m2
      hutop(i,j)=(aip(i,j)*u(i+1,j)+aim(i,j)*u(i-1,j)
     1  +ajp(i,j)*u(i,j+1)+ajm(i,j)*u(i,j-1)+con(i,j)
c     1    -rdtm*rho(i,j)*vol(i,j)/void(i,j)*fo(i,j,1)
     1 )/vol(i,j)
  215 continue
	write(*,*)'hutop=',hutop(124,2)
      call reset
      go to 100
c
c --- (3-36) momentum interpolation ---
c
  200 if(nf.ne.2) go to 600

      do 210 i=2,l2
      do 210 j=2,m2
      apt=rdtm*rho(i,j)*vol(i,j)/void(i,j)

      barbeta=1.0
      con(i,j)=con(i,j)+fo(i,j,2)*apt
c     1    +barbeta*rho(i,j)*grash*pr*pr*(f(i,j,5)-1.0)*vol(i,j)
c     ap(i,j)=apt+ex1(i,j)/void(i,j)/void(i,j)
      ap(i,j)=ap(i,j)+apt

	if(id_ele(i,j).eq.3)then
	ap(i,j)=ap(i,j)
c  +1/darcy*vol(i,j)
	endif
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)

  210 continue
	call bound



c for 2:v

      call flux(0)
      do 315 i=2,l2
      do 315 j=2,m2
      hvtop(i,j)=(aip(i,j)*v(i+1,j)+aim(i,j)*v(i-1,j)
     1   +ajp(i,j)*v(i,j+1)
     1           +ajm(i,j)*v(i,j-1)+con(i,j)
c     1    -rdtm*rho(i,j)*vol(i,j)/void(i,j)*fo(i,j,2)
     1   )/vol(i,j)
      dxdksi=0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
      dydksi=0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
      hksip(i,j)=0.5*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
      dxdeta=0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
      dydeta=0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
      hetap(i,j)=0.5*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
  315 continue
c
c --- (3-36) finish ---

      call reset
      go to 100

c for pressure correction pc
500	continue
	goto 400
      call aread(1)
      do 510 i=2,l2
      do 510 j=2,m2
      apt=rdtm*rho(i,j)*vol(i,j)/void(i,j)

510   con(i,j)=(con(i,j)+apt*u(i,j))
      nf=1
      call flux(0)
      nf=npc
      do 515 i=2,l2
      do 515 j=2,m2
      hutop(i,j)=(aip(i,j)*u(i+1,j)+aim(i,j)*u(i-1,j)
     1   +ajp(i,j)*u(i,j+1)
     1          +ajm(i,j)*u(i,j-1)+con(i,j))/vol(i,j)
515	continue
      call reset
      call aread(2)
      do 520 i=2,l2
      do 520 j=2,m2
      apt=rdtm*rho(i,j)*vol(i,j)/void(i,j)

520   con(i,j)=(con(i,j)+apt*v(i,j))
      nf=2
      call flux(0)
      nf=npc
      do 525 i=2,l2
      do 525 j=2,m2
      hvtop(i,j)=(aip(i,j)*v(i+1,j)+aim(i,j)*v(i-1,j)
     1   +ajp(i,j)*v(i,j+1)
     1           +ajm(i,j)*v(i,j-1)+con(i,j))/vol(i,j)
      dxdksi=0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
      dydksi=0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
      hksip(i,j)=0.5*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
      dxdeta=0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
      dydeta=0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
      hetap(i,j)=0.5*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
525	continue

c --- (3-36) finish ---

      call reset


c-----coefficients for p and p prime eqns--------
  400 continue

c --- hksip, hetap are based on momentum interpolation.
c --- and t3 is based on the two-dimenisonal correction of MIS  
c --- and t4 is based on the secondary correction of 2-d of MIS  
c --- mis = 0 using all the correction. mis = 1 using one correction
c --- mis = 2 no correction. 
	call aread(1)
      mis=0
      do 410 i=2,l2
      do 410 j=2,m2
	con(i,j)=0
c	con(i,j)=ex1(i,j)
c-(rho(i,j)-rho_o(i,j))*rdtm*vol(i,j)
c see p35, Ap
      ap(i,j)=(aip(i,j)+aim(i,j)+ajm(i,j)+ajp(i,j))/vol(i,j)/rho(i,j)
c	if(i.eq.12.and.j.eq.70)write(*,*)'ap',rho(i,j),vol(i,j),
c     1   ap(i,j),
c     1    aip(i,j),aim(i,j),ajp(i,j),ajm(i,j)
	tmp(i,j)=ap(i,j)
	if(j.eq.2.and.i.eq.2)then
	write(*,*)aip(i,j),aim(i,j),ajm(i,j),ajp(i,j)
	endif
410	continue
      if(.not.lortho) call extra1
      do 422 j=2,m2
      do 420 i=2,l3
      a1=0.5*(ak1(i,j)+ak1(i+1,j))
      a2=0.5*(ak1(i+1,j)+ak1(i+2,j))
c see p36 equ.31
	x1=ap(i,j)*hksi(i,j)/a1
	x2=ap(i+1,j)*hksi(i+1,j)/a2
	denom=0.5*(x1+x2)
c	if(id_bound(i,j).ne.1.or.id_bound(i+1,j).ne.1)then
c	denom=2*x1*x2/(x1+x2)
c	endif
	if(nf.eq.npc)goto 419
      apr=2.0*denom/(hksi(i,j)+hksi(i+1,j))
      dxdksi=0.5*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
      dydksi=0.5*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
      rukij=(dxdksi*u(i,j)+dydksi*v(i,j))*rho(i,j)/hksi(i,j)*a1
      dxdksi=0.5*(xc(i+2,j+1)+xc(i+2,j)-xc(i+1,j+1)-xc(i+1,j))
      dydksi=0.5*(yc(i+2,j+1)+yc(i+2,j)-yc(i+1,j+1)-yc(i+1,j))
      ruki1j=(dxdksi*u(i+1,j)+dydksi*v(i+1,j))*rho(i+1,j)
     1   /hksi(i+1,j)*a2
      t1=rukij*(apr*hksi(i+1,j)-ap(i,j)*hksi(i,j)/a1)
      t2=ruki1j*(apr*hksi(i,j)-ap(i+1,j)*hksi(i+1,j)/a2)
      t3=0.5*(t1+t2)
      rukt=ruksi(2,j)
      if(i.ne.2) rukt=temp
      temp=ruksi(i+1,j)
      frc=hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
      t4=0.5*(frc*(ruksi(i+1,j)*ak1(i+1,j)-rukt*ak1(i,j))
     1  +(1.0-frc)*(ruksi(i+1,j)*ak1(i+1,j)-ruksi(i+2,j)*ak1(i+2,j)))
      if(mis.eq.2)t3=0
      if(mis.eq.1)t4=0
      if(mis.eq.2)t4=0
	t4=0
	if(id_bound(i,j).eq.1.or.id_bound(i+1,j).eq.1)then
      ruksi(i+1,j)=((hksip(i,j)+hksip(i+1,j)+t3)/denom+t4)/ak1(i+1,j)
	else
	ruksi(i+1,j)=(rukij*x1+ruki1j*x2)/(x1+x2)
	endif
      if(nf.eq.npc) ruksi(i+1,j)=ruksi(i+1,j)+(p(i,j)-p(i+1,j))/denom
     1             /ak1(i+1,j)
419	continue
      aip(i,j)=1.0/denom
  420 aim(i+1,j)=aip(i,j)
      aip(l2,j)=0.0
  422 aim(2,j)=0.0
	write(*,*)'ruksi=',ruksi(i_rt,2),ruksi(i_rt+1,2),
     1   u(i_rt,2),u(i_rt+1,2)

      do 432 i=2,l2
      do 430 j=2,m3
      a1=0.5*(ae1(i,j)+ae1(i,j+1))
      a2=0.5*(ae1(i,j+1)+ae1(i,j+2))
	x1=ap(i,j)*heta(i,j)/a1
	x2=ap(i,j+1)*heta(i,j+1)/a2
	denom=0.5*(x1+x2)
c	if(id_bound(i,j).ne.1.or.id_bound(i,j+1).ne.1)then
c	denom=2*x1*x2/(x1+x2)
c	endif
	if(nf.eq.npc)goto 429
      apr=2.0*denom/(heta(i,j)+heta(i,j+1))
      dxdeta=0.5*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
      dydeta=0.5*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
      rueij=(dxdeta*u(i,j)+dydeta*v(i,j))*rho(i,j)/heta(i,j)*a1
      dxdeta=0.5*(xc(i+1,j+2)+xc(i,j+2)-xc(i+1,j+1)-xc(i,j+1))
      dydeta=0.5*(yc(i+1,j+2)+yc(i,j+2)-yc(i+1,j+1)-yc(i,j+1))
      rueij1=(dxdeta*u(i,j+1)+dydeta*v(i,j+1))*rho(i,j+1)
     1   /heta(i,j+1)*a2
      t1=rueij*(apr*heta(i,j+1)-ap(i,j)*heta(i,j)/a1)
      t2=rueij1*(apr*heta(i,j)-ap(i,j+1)*heta(i,j+1)/a2)
      t3=0.5*(t1+t2)
      ruet=rueta(i,2)
      if(j.ne.2) ruet=temp
      temp=rueta(i,j+1)
      frc=heta(i,j+1)/(heta(i,j)+heta(i,j+1))
      t4=0.5*(frc*(rueta(i,j+1)*ae1(i,j+1)-ruet*ae1(i,j))
     1  +(1.0-frc)*(rueta(i,j+1)*ae1(i,j+1)-rueta(i,j+2)*ae1(i,j+2)))
      if(mis.eq.2)t3=0
      if(mis.eq.1)t4=0
      if(mis.eq.2)t4=0
	t4=0
	if(id_bound(i,j).eq.1.or.id_bound(i,j+1).eq.1)then
      rueta(i,j+1)=((hetap(i,j)+hetap(i,j+1)+t3)/denom+t4)/ae1(i,j+1)
	else
	rueta(i,j+1)=(rueij*x1+rueij1*x2)/(x1+x2)
	endif
      if(nf.eq.npc) rueta(i,j+1)=rueta(i,j+1)+(p(i,j)-p(i,j+1))/denom
     1             /ae1(i,j+1)
c	if(i.eq.12.and.j.eq.69)write(*,*)'denom',denom,hetap(i,j),
c     1   hetap(i,j+1),t3,t4,ae1(i,j+1),ap(i,j),ap(i,j+1)
429	continue
      ajp(i,j)=1.0/denom
  430 ajm(i,j+1)=ajp(i,j)
      ajp(i,m2)=0.0
  432 ajm(i,2)=0.0
c	write(*,*)'after p,nf,rueta(12,69)',nf,rueta(12,69),
c     1   rueta(12,70)

      do 440 i=2,l2
      do 440 j=2,m2
      con(i,j)=con(i,j)+ak1(i,j)*ruksi(i,j)-ak1(i+1,j)*ruksi(i+1,j)
     1          +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
  440 ap(i,j)=aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
	call bound

	if(nf.eq.np)then
	i=i_rb
	j=j_rr
	write(*,*)'test p',ap(i,j),aip(i,j),aim(i,j),ajp(i,j),ajm(i,j),
     1    smass(i,j),rueta(i,j),rueta(i,j+1)
	endif
c only for pressure and pressure correction

	write(*,*)'p4=',pc(i_rb-10,1),pc(i_rb-10,2)
c      fmax1=0
c      fmax2=0
c      do 445 i=1,l1
c      do 445 j=1,m1
c      fij1=abs(ruksi(i,j))
c      fij2=abs(rueta(i,j))
c	fmax1=dmax1(fmax1,fij1)
c445	fmax2=dmax1(fmax2,fij2)
c	write(*,*)'nf,ruksi,rueta=',nf,fmax1,fmax2
c --- finish pressure and pressure correction equation coeff. and 
c --- source term calculation ---------

      if(nf.ne.np) go to 704
      smax=0.0
      ssum=0.0
      ismax=0
      jsmax=0
      do 442 i=2,l2
      do 442 j=2,m2
      soor=ap(i,j)*p(i,j)-aip(i,j)*p(i+1,j)-aim(i,j)*p(i-1,j)-
     1      ajp(i,j)*p(i,j+1)-ajm(i,j)*p(i,j-1)-con(i,j)
      ssum=ssum+soor/rho(i,j)
      smax=dmax1(abs(soor),smax)
      if(smax.eq.abs(soor)) then
      ismax=i
      jsmax=j
      endif
  442 continue

      do 450 i=2,l2
      do 450 j=2,m2
      ap(i,j)=ap(i,j)/relax(np)
  450 con(i,j)=con(i,j)+(1.0-relax(np))*ap(i,j)*p(i,j)

  704 continue

c for pressure and pc
      call solve

	write(*,*)'p5=',p(1,1),p(1,2),p(1,3)

      fmax=0
      do 455 i=1,l1
      do 455 j=1,m1
      fij=abs(f(i,j,nf))
455	fmax=dmax1(fmax,fij)
	write(*,*)'in setup2 : nf,fmax=',nf,fmax,f(1,1,nf)

c --- correct ruksi and rueta---------
c --- (3.34) and (3.35)

      do 460 j=2,m2
      do 460 i=2,l3
	if(id_bound(i,j).eq.1.and.id_bound(i+1,j).eq.1)then
	ruksi(i+1,j)=ruksi(i+1,j)+aip(i,j)/ak1(i+1,j)*(f(i,j,nf)-
     1             f(i+1,j,nf))
	endif
460	continue
      do 470 i=2,l2
      do 470 j=2,m3
	if(id_bound(i,j).eq.1.and.id_bound(i,j+1).eq.1)then
	rueta(i,j+1)=rueta(i,j+1)+ajp(i,j)/ae1(i,j+1)*(f(i,j,nf)-
     1             f(i,j+1,nf))
	endif
470	continue

	do i=2,l2
	do j=2,m2
	smass(i,j)=ak1(i,j)*ruksi(i,j)-ak1(i+1,j)*ruksi(i+1,j)
     1    +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
c     1    +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)+ex1(i,j)
	enddo
	enddo
	fmax=1d-6
	do i=2,l2
	do j=2,m2
	fmax1=fmax
	fij=dabs(smass(i,j))
	fmax=dmax1(fmax,fij)
	if(fij.gt.fmax1)then
	i1=i
	j1=j
	endif
	enddo
	enddo
	ssum=0
	do i=2,l2
	do j=2,m2
	ssum=ssum+dabs(smass(i,j))
	enddo
	enddo
	smax=smass(i1,j1)
	write(*,*)'i1,j1,smass,ssum',i1,j1,smass(i1,j1),ssum
c-----interpolate to get pressure on the interface-------

	call boun_b

	write(*,*)'ps=',p(1,1),p(1,2),p(1,3)
      do 518 j=2,m2
      do 508 i=2,l3
c	write(*,*)i,j
	if(id_bound(i,j).ne.1.or.id_bound(i+1,j).ne.1)then
	em=0.5*hksi(i,j)
	ep=0.5*hksi(i+1,j)
      de=em+ep
c	write(*,*)em,ep,de
	pksi(i+1,j)=(f(i,j,nf)*ep+f(i+1,j,nf)*em)/de
	else
      bpi=p(i,j)
      bpi1=p(i+1,j)-(hksip(i,j)+hksip(i+1,j))
c see equ 31
      a1=0.5*(ak1(i,j)+ak1(i+1,j))
      a2=0.5*(ak1(i+1,j)+ak1(i+2,j))
      em=0.5*hksi(i,j)*tmp(i,j)/a1
      ep=0.5*hksi(i+1,j)*tmp(i+1,j)/a2
      de=em+ep
      if(nf.eq.np)then
c in case of infinite tmp(i,j=2)
      t1=0.5*em*ep*(ruksi(i+2,j)*ak1(i+2,j)-ruksi(i,j)*ak1(i,j))/de
      bpm=(bpi*ep+bpi1*em)/de+t1
      pksi(i+1,j)=bpm+hksip(i,j)
	else if(nf.eq.npc)then
	pksi(i+1,j)=(f(i,j,nf)*ep+f(i+1,j,nf)*em)/de
	endif
	endif
	if(i.eq.2.and.j.eq.2)then
	write(*,*)tmp(i,j),tmp(i+1,j),hksi(i,j),hksi(i+1,j),em,ep
	endif
508	continue
      pksi(l1,j)=2.0*f(l2,j,nf)-pksi(l2,j)
      pksi(2,j)=2.0*f(2,j,nf)-pksi(3,j)
      f(1,j,nf)=pksi(2,j)
518	f(l1,j,nf)=pksi(l1,j)

	write(*,*)ruksi(2,2),ruksi(3,2),pksi(3,2),p(2,2)
	write(*,*)'ps=',p(1,1),p(1,2),p(1,3)
      do 534 i=2,l2
      do 530 j=2,m3
	if(id_bound(i,j).ne.1.or.id_bound(i,j+1).ne.1)then
	em=0.5*heta(i,j)
	ep=0.5*heta(i,j+1)
      de=em+ep
	peta(i,j+1)=(f(i,j,nf)*ep+f(i,j+1,nf)*em)/de
	else
      bpj=p(i,j)
      bpj1=p(i,j+1)-(hetap(i,j)+hetap(i,j+1))
      a1=0.5*(ae1(i,j)+ae1(i,j+1))
      a2=0.5*(ae1(i,j+1)+ae1(i,j+2))
      em=0.5*heta(i,j)*tmp(i,j)/a1
      ep=0.5*heta(i,j+1)*tmp(i,j+1)/a2
      de=em+ep
      if(nf.eq.np)then
      t1=0.5*em*ep*(rueta(i,j+2)*ae1(i,j+2)-rueta(i,j)*ae1(i,j))/de
      bpm=(bpj*ep+bpj1*em)/de+t1
      peta(i,j+1)=bpm+hetap(i,j)
	else if(nf.eq.npc)then
	peta(i,j+1)=(f(i,j,nf)*ep+f(i,j+1,nf)*em)/de
	endif
	endif
530	continue
      peta(i,m1)=2.0*f(i,m2,nf)-peta(i,m2)
      peta(i,2)=2.0*f(i,2,nf)-peta(i,3)
      f(i,1,nf)=peta(i,2)
534	f(i,m1,nf)=peta(i,m1)
	write(*,*)'p6=',p(1,1),p(1,2),p(1,3)

c --- finish pressure calculation and interpolation ---
c-----solve u----------------------------------

      if(nf.eq.npc) go to 705
      call aread(1)
      do 550 i=2,l2
      if(gam(i,1).eq.0) ajm(i,2)=0.0
  550 if(gam(i,m1).eq.0) ajp(i,m2)=0.0
      do 555 j=2,m2
      if(gam(1,j).eq.0) aim(2,j)=0.0
  555 if(gam(l1,j).eq.0) aip(l2,j)=0.0
  705 continue
      do 560 i=2,l2
      do 560 j=2,m2
      xr=0.5*(xc(i+1,j+1)+xc(i+1,j))
      yr=0.5*(yc(i+1,j+1)+yc(i+1,j))
      xl=0.5*(xc(i,j+1)+xc(i,j))
      yl=0.5*(yc(i,j+1)+yc(i,j))
      xu=0.5*(xc(i+1,j+1)+xc(i,j+1))
      yu=0.5*(yc(i+1,j+1)+yc(i,j+1))
      xd=0.5*(xc(i+1,j)+xc(i,j))
      yd=0.5*(yc(i+1,j)+yc(i,j))
      dxdksi=xr-xl
      dydksi=yr-yl
      dxdeta=xu-xd
      dydeta=yu-yd
      dpksi=pksi(i+1,j)-pksi(i,j)
      dpeta=peta(i,j+1)-peta(i,j)
      t1=-dydeta*dpksi+dydksi*dpeta
      t2=dxdeta*dpksi-dxdksi*dpeta
      if(mode.eq.1) t1=t1*r(i,j)
      if(mode.eq.1) t2=t2*r(i,j)
      if(nf.ne.npc) go to 561

c --- correct velocity or step 6 ---
	if(id_bound(i,j).eq.1)then
      u(i,j)=u(i,j)+t1/tmp(i,j)/rho(i,j)/vol(i,j)
      v(i,j)=v(i,j)+t2/tmp(i,j)/rho(i,j)/vol(i,j)
	endif
      go to 560
  561 continue
      tmp(i,j)=t2
	if(id_bound(i,j).eq.1)then
      con(i,j)=con(i,j)+t1
	endif
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
  560 continue
	write(*,*)'p7=',p(1,1),p(1,2),p(1,3)
      if(nf.eq.npc) go to 100
      nf=1
      call flux(1)

c-----solve v----------------------------------
      call aread(2)
      nf=2
      do 570 i=2,l2
      if(gam(i,1).eq.0) ajm(i,2)=0.0
  570 if(gam(i,m1).eq.0) ajp(i,m2)=0.0
      do 575 j=2,m2
      if(gam(1,j).eq.0) aim(2,j)=0.0
  575 if(gam(l1,j).eq.0) aip(l2,j)=0.0
      do 580 i=2,l2
      do 580 j=2,m2
	if(id_bound(i,j).eq.1)then
      con(i,j)=con(i,j)+tmp(i,j)
	endif
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
  580 continue
      call flux(1)
      go to 100

c-------------------------------------------------

  600 continue
c for nf=5:t
      if(nf.ne.5)goto 1000
c for temperature
      do 310 i=1,l1
      do 310 j=1,m1
c rho_s here is different from rho
	cp(i,j)=heat_cap(i,j)
      apt=rdtm*rho_s(i,j)*cp(i,j)*vol(i,j)

      con(i,j)=con(i,j)+fo(i,j,nf)*apt
c     ap(i,j)=apt+ex1(i,j)
c It is difficult to consider the source term
      ap(i,j)=ap(i,j)+apt
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
310	continue
	call bound



c	if(nf.eq.6)then
c	write(*,*)'ccen',ap(i_rb,2),aip(i_rb,2),aim(i_rb,2),
c     1     ajp(i_rb,2),ajm(i_rb,2),con(i_rb,2),
c     1    ccen(i_rb+1,2),ccen(i_rb,2),ccen(i_rb-1,2)
c	endif
      call flux(1)
c	if(nf.eq.6)then
c	write(*,*)'ccen',ap(i_rb,2),aip(i_rb,2),aim(i_rb,2),
c     1     ajp(i_rb,2),ajm(i_rb,2),con(i_rb,2),
c     1   ccen(i_rb+1,2),ccen(i_rb,2),ccen(i_rb-1,2)
c	endif

1000  continue
      if(nf.ne.6)goto 100
c concentration
	do i=1,l1
	do j=1,m1
	apt=rdtm*0.0521*void(i,j)*vol(i,j)
c molecular weight of SiC2, Mt=0.0521 kg/mol

c	ap(i,j)=0
	con(i,j)=con(i,j)+fo(i,j,nf)*apt
c	ap(i,j)=ap(i,j)+ex1(i,j)
	ap(i,j)=ap(i,j)+apt
	if(id_ele(i,j).eq.3)then
	con(i,j)=con(i,j)+scargon(i,j)*0.0521
c P (SiC2)=R*T*C(SiC2)
	sour1=6/0.0016*0.606/sqrt(t(i,j))*
     1    8.31*t(i,j)*0.0521*vol(i,j)
	ap(i,j)=ap(i,j)+sour1
	endif
	enddo
	enddo
c	write(*,*)'ccen',ap(i_rb,2),aip(i_rb,2),aim(i_rb,2),
c     1     ajp(i_rb,2),ajm(i_rb,2),con(i_rb,2),ccen(i_rb,2),
c     1   fo(i_rb,2,6),rho(i_rb,2)

      do i=1,l1
      do j=1,m1
      ap(i,j)=ap(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)

      enddo
      enddo
	call bound
	call flux(0)
  100 continue
      return
      end

c ------------------------------------------------------------------
      subroutine flux(isol)

	include "sc.h"
      data fl00/0./

c	write(*,*)'nf,rueta(12,69)',nf,rueta(12,69)
      do 10 i=2,l2
      do 10 j=2,m2
      con0(i,j)=con(i,j)
   10 ap0(i,j)=ap(i,j)

	goto 30
      if(iter.ne.0) go to 30
      if(isol.eq.0) return

   30 continue

      do 999 nt_f=1,ntimes(nf)
      if(nt_f.ne.1) go to 998

c--- find flux,needed for the first iteration in the ntimes loop ---

      do 300 j=1,m1
      diff=gam(2,j)/(0.5*hksi(2,j))
      flow=-ruksi(2,j)
	vv=void(2,j)
      call diflow(vv,2,j)
      qt(2)=-flow*f(2,j,nf)+acof*(f(1,j,nf)-f(2,j,nf))
      do 320 i=3,l2
      jj=j
      if(j.eq.1) jj=2
      if(j.eq.m1) jj=m2
      diff=diff_i(i-1,jj)
      flow=-ruksi(i,j)
	vv=dmax1(void(i-1,j),void(i,j))
      call diflow(vv,i,j)
      qt(i)=-flow*f(i,j,nf)+acof*(f(i-1,j,nf)-f(i,j,nf))
  320 continue
      diff=gam(l2,j)/(0.5*hksi(l2,j))
      flow=ruksi(l1,j)
	vv=void(l1,j)
      call diflow(vv,l1,j)
      qt(l1)=flow*f(l2,j,nf)+acof*(f(l2,j,nf)-f(l1,j,nf))
      do 330 i=2,l1
      fjksi(i,j)=qt(i)
  330 continue
  300 continue

c---- along j ------

      do 200 i=1,l1
      diff=gam(i,2)/(0.5*heta(i,2))
      flow=-rueta(i,2)
	vv=void(i,2)
      call diflow(vv,i,2)
      qt(2)=-flow*f(i,2,nf)+acof*(f(i,1,nf)-f(i,2,nf))
      do 220 j=3,m2
      ii=i
      if(i.eq.1) ii=2
      if(i.eq.l1) ii=l2
      diff=diff_j(ii,j-1)
      flow=-rueta(i,j)
	vv=dmax1(void(i,j-1),void(i,j))
      call diflow(vv,i,j)
      qt(j)=-flow*f(i,j,nf)+acof*(f(i,j-1,nf)-f(i,j,nf))
  220 continue
      diff=gam(i,m2)/(0.5*heta(i,m2))
      flow=rueta(i,m1)
	vv=void(i,m1)
      call diflow(vv,i,m1)
      qt(m1)=flow*f(i,m2,nf)+acof*(f(i,m2,nf)-f(i,m1,nf))
      do 230 j=2,m1
      fjeta(i,j)=qt(j)
  230 continue
  200 continue

	goto 997
c --- finish jksi and jeta calculation.
c------ contstract boundary flux for j type b.c. -------
c --- boundary condition con(1,j),*,*,* go into jksi and jeta
c --- ak2/ak1,* is project on xi - direction from eta - direction
c --- variable, in program let ak2=0.0, because sym. orth.

      do 62 j=2,m2
      if(gam(1,j).ne.0) go to 61
	if(ak1(2,j).eq.0)then
	fjksi(2,j)=0
	else
      fjksi(2,j)=(con(1,j)*heta(1,j)+0.5*ak2(2,j)*(fjeta(1,j)+
     &          fjeta(1,j+1)))/ak1(2,j)
	endif
   61 continue
      if(gam(l1,j).ne.0) go to 62
      fjksi(l1,j)=(con(l1,j)*heta(l1,j)+0.5*ak2(l1,j)*(fjeta(l1,j)+
     &          fjeta(l1,j+1)))/ak1(l1,j)
   62 continue

      do 64 i=2,l2
      if(gam(i,1).ne.0) go to 63
	if(ae1(i,2).eq.0)then
	fjeta(i,2)=0
	else
      fjeta(i,2)=(con(i,1)*hksi(i,1)+0.5*ae2(i,2)*(fjksi(i,1)+
     &          fjksi(i+1,1)))/ae1(i,2)
	endif
   63 continue
      if(gam(i,m1).ne.0) go to 64
      fjeta(i,m1)=(con(i,m1)*hksi(i,m1)+0.5*ae2(i,m1)*(fjksi(i,m1)+
     &          fjksi(i+1,m1)))/ae1(i,m1)
   64 continue
c	goto 700
c	goto 997   
c --- (2-31) and (2-42) find out vector J cdot vector n or J_n in bound.

      do 881 i=2,l2
      t1=0.5*(fjksi(i,1)+fjksi(i+1,1))
      t2=0.5*(fjksi(i,m1)+fjksi(i+1,m1))
      fjbj1(i,nf)=(fjeta(i,2)*ae1(i,2)-t1*ae2(i,2))/hksi(i,1)
      fjbm1(i,nf)=(fjeta(i,m1)*ae1(i,m1)-t2*ae2(i,m1))/hksi(i,m1)
	if(mode.eq.1)then
	if(r(i,2).eq.0)then
	fjbj1(i,nf)=0
	else
	fjbj1(i,nf)=fjbj1(i,nf)/r(i,2)
	endif
	fjbm1(i,nf)=fjbm1(i,nf)/r(i,m1)
	endif
  881 continue
      do 882 j=2,m2
      t1=0.5*(fjeta(1,j)+fjeta(1,j+1))
      t2=0.5*(fjeta(l1,j)+fjeta(l1,j+1))
      fjbi1(j,nf)=(fjksi(2,j)*ak1(2,j)-t1*ak2(2,j))/heta(1,j)
      fjbl1(j,nf)=(fjksi(l1,j)*ak1(l1,j)-t2*ak2(l1,j))/heta(l1,j)
	if(mode.eq.1)then
	if(r(2,j).eq.0)then
	fjbi1(j,nf)=0
	else
	fjbi1(j,nf)=fjbi1(j,nf)/r(2,j)
	endif
	fjbl1(j,nf)=fjbl1(j,nf)/r(l1,j)
	endif
  882 continue

      go to 997

  998 continue
	goto 997
c--- correct j from known phi and j
c------ correct jksi ------

      do 50  j=2,m2
      if(gam(1,j).eq.0) go to 51
      t1=aim(2,j)/ak1(2,j)*(f(1,j,nf)-f(2,j,nf))+ruksi(2,j)*f(2,j,nf)
      fjksi(2,j)=t1
   51 continue
      do 50 i=3,l1
      if(gam(l1,j).eq.0.and.i.eq.l1) go to 50
      t1=aip(i-1,j)/ak1(i,j)*(f(i-1,j,nf)-f(i,j,nf))+ruksi(i,j)*
     &   f(i-1,j,nf)
      fjksi(i,j)=t1
   50 continue

c--- correct jeta  ------

      do 60 i=2,l2
      if(gam(i,1).eq.0) go to 67
      t1=ajm(i,2)/ae1(i,2)*(f(i,1,nf)-f(i,2,nf))+rueta(i,2)*f(i,2,nf)
      fjeta(i,2)=t1
   67 continue
      do 60 j=3,m1
      if(gam(i,m1).eq.0.and.j.eq.m1) go to 60
      t1=ajp(i,j-1)/ae1(i,j)*(f(i,j-1,nf)-f(i,j,nf))+rueta(i,j)*
     &   f(i,j-1,nf)
      fjeta(i,j)=t1
   60 continue

  997 continue

c--- reconstruct con and ap

      do 65 i=2,l2
      do 65 j=2,m2
      con(i,j)=con0(i,j)
   65 ap(i,j)=ap0(i,j)
      if(.not.lortho) call extra2
c	goto 700

c---- calculate jhat ----
c --- step 5 use equation (2-97) and (2-102), (2-103), (2-104) ---
c --- calculate boundary data from J and other purpose for b_sp. ---
c --- from here to end. pc = jksi hat or jeta hat ----
c -----jksi hat-----
c
      do 70 j=2,m2
      if(gam(1,j).eq.0)then
      diff=gam(2,j)/(0.5*hksi(2,j))
      flow=-ruksi(2,j)
	vv=void(2,j)
      call diflow(vv,2,j)
      t1=fjksi(2,j)+flow*f(2,j,nf)
      t2=acof+1.e-10
      f(1,j,nf)=f(2,j,nf)+t1/t2
	else
	fjksi(2,j)=0.0
	endif
c      pc(2,j)=fjksi(2,j)
   71 continue
      if(gam(l1,j).eq.0)then
      diff=gam(l2,j)/(0.5*hksi(l2,j))
      flow=ruksi(l1,j)
	vv=void(l1,j)
      call diflow(vv,l1,j)
      t1=fjksi(l1,j)-flow*f(l2,j,nf)
      t2=acof+1.e-10
	if(t2.eq.0)then
	t12=0
	else
	t12=t1/t2
	endif
      f(l1,j,nf)=f(l2,j,nf)-t12
c      pc(l1,j)=fjksi(l1,j)
	else
	fjksi(l1,j)=0.0
	endif
   72 continue
      do 70 i=3,l2
      fjksi(i,j)=0.0
   70 continue

c      do 80 j=2,m2
c      do 80 i=2,l1
c      fjksi(i,j)=pc(i,j)
c   80 continue

c---jeta hat ------

      do 90 i=2,l2
      if(gam(i,1).eq.0)then
      diff=gam(i,2)/(0.5*heta(i,2))
      flow=-rueta(i,2)
	vv=void(i,2)
      call diflow(vv,i,2)
      t1=fjeta(i,2)+flow*f(i,2,nf)
      t2=acof+1.e-10
      f(i,1,nf)=f(i,2,nf)+t1/t2
c      pc(i,2)=fjeta(i,2)
	else
	fjeta(i,2)=0.0
	endif
   91 continue
      if(gam(i,m1).eq.0)then
      diff=gam(i,m2)/(0.5*heta(i,m2))
      flow=rueta(i,m1)
	vv=void(i,m1)
      call diflow(vv,i,m1)
      t1=fjeta(i,m1)-flow*f(i,m2,nf)
      t2=acof+1.e-10
	if(t2.eq.0)then
	t12=0
	else
	t12=t1/t2
	endif
      f(i,m1,nf)=f(i,m2,nf)-t12
c      pc(i,m1)=fjeta(i,m1)
	else
	fjeta(i,m1)=0.0
	endif
   92 continue
      do 90 j=3,m2
	fjeta(i,j)=0.0
   90 continue

c      do 100 i=2,l2
c      do 100 j=2,m1
c      fjeta(i,j)=pc(i,j)
c  100 continue

c--- calculate and solve phy equation
c --- (2-106) and (2-107) b = b_s + b_no + b_sp, t1=b_sp ---

      rel=1.0-relax(nf)
      do 120 i=2,l2
      do 120 j=2,m2
      t1=fjksi(i,j)*ak1(i,j)-fjksi(i+1,j)*ak1(i+1,j)+fjeta(i,j)*
     &    ae1(i,j)-fjeta(i,j+1)*ae1(i,j+1)
      con(i,j)=con(i,j)+t1
c	write(*,*)con(i,j)
  120 continue

700	if(isol.eq.0) return

      do 125 i=2,l2
      do 125 j=2,m2
      ap(i,j)=ap(i,j)/relax(nf)
      con(i,j)=con(i,j)+rel*ap(i,j)*f(i,j,nf)
  125 continue

      call solve

  999 continue

	write(*,*)nf,i_rt,j_rr,ccen(i_rt,j_rr-1)

      res(nf)=0.0
      fmax=0
      do 400 i=1,l1
      do 400 j=1,m1
      fij=abs(f(i,j,nf))
  400 fmax=dmax1(fmax,fij)
	write(*,*)'nf,fmax=',nf,fmax
  405 continue
      do 410 i=1,l1
      do 410 j=1,m1
      foij=fn(i,j,nf)
      fij=f(i,j,nf)
      dffo=abs(fij-foij)
      error=dffo/dmax1(fmax,1d-5)
  410 res(nf)=dmax1(res(nf),error)
	epsil=0.01
	if(nf.eq.1.or.nf.eq.2)epsil=0.04
	write(*,*)'res(nf)=',res(nf)
	if(res(nf).gt.epsil.or.iter.le.6)then
	lconv=.false.
	else
	lcon(nf)=.true.
	lconv=.true.
	endif
      do 420 i=1,l1
      do 420 j=1,m1
  420 fn(i,j,nf)=f(i,j,nf)
	write(*,*)ccen(i_rt,j_rr-1)

      return
      end

	subroutine flux_s
	include "sc.h"
      do 300 j=1,m1
      diff=gam(2,j)/(0.5*hksi(2,j))
      flow=-ruksi(2,j)
	vv=void(2,j)
      call diflow(vv,2,j)
      qt(2)=-flow*f(2,j,nf)+acof*(f(1,j,nf)-f(2,j,nf))
      do 320 i=3,l2
      jj=j
      if(j.eq.1) jj=2
      if(j.eq.m1) jj=m2
      diff=diff_i(i-1,jj)
      flow=-ruksi(i,j)
	vv=dmax1(void(i-1,j),void(i,j))
      call diflow(vv,i,j)
      qt(i)=-flow*f(i,j,nf)+acof*(f(i-1,j,nf)-f(i,j,nf))
  320 continue
      diff=gam(l2,j)/(0.5*hksi(l2,j))
      flow=ruksi(l1,j)
	vv=void(l1,j)
      call diflow(vv,l1,j)
      qt(l1)=flow*f(l2,j,nf)+acof*(f(l2,j,nf)-f(l1,j,nf))
      do 330 i=2,l1
      fjksi(i,j)=qt(i)
  330 continue
  300 continue

c---- along j ------

      do 200 i=1,l1
      diff=gam(i,2)/(0.5*heta(i,2))
      flow=-rueta(i,2)
	vv=void(i,2)
      call diflow(vv,i,2)
      qt(2)=-flow*f(i,2,nf)+acof*(f(i,1,nf)-f(i,2,nf))
      do 220 j=3,m2
      ii=i
      if(i.eq.1) ii=2
      if(i.eq.l1) ii=l2
      diff=diff_j(ii,j-1)
      flow=-rueta(i,j)
	vv=dmax1(void(i,j-1),void(i,j))
      call diflow(vv,i,j)
      qt(j)=-flow*f(i,j,nf)+acof*(f(i,j-1,nf)-f(i,j,nf))
  220 continue
      diff=gam(i,m2)/(0.5*heta(i,m2))
      flow=rueta(i,m1)
	vv=void(i,m1)
      call diflow(vv,i,m1)
      qt(m1)=flow*f(i,m2,nf)+acof*(f(i,m2,nf)-f(i,m1,nf))
      do 230 j=2,m1
      fjeta(i,j)=qt(j)
  230 continue
  200 continue
	return
	end

c -------------------------------------------------------------------
      subroutine extra
c -------------------------------------------------------------------
	include "sc.h"

      entry extra1
c --- construct extra con for non-orthogonal coordinate system ---
c --- calculate b_no for pressure and pressure correction equation.
c --- b_no = 0.0 for ortho. gird. so control by lortho = .false.

      do 100 j=2,m2
      do 110 i=1,l1
      qt(i)=0.5*(rueta(i,j)+rueta(i,j+1))
  110 continue
      do 100 i=2,l2
      tem1=hksi(i+1,j)+hksi(i,j)
      tem2=hksi(i-1,j)+hksi(i,j)
      t1=hksi(i,j)/tem1
      t2=hksi(i+1,j)/tem1
      t3=hksi(i,j)/tem2
      t4=hksi(i-1,j)/tem2
      con(i,j)=con(i,j)+(t1*qt(i+1)+t2*qt(i))*ak2(i+1,j)
     1        -(t4*qt(i)+t3*qt(i-1))*ak2(i,j)
  100 continue

      do 200 i=2,l2
      do 210 j=1,m1
      qt(j)=0.5*(ruksi(i,j)+ruksi(i+1,j))
  210 continue
      do 200 j=2,m2
      tem1=heta(i,j+1)+heta(i,j)
      tem2=heta(i,j-1)+heta(i,j)
      t1=heta(i,j)/tem1
      t2=heta(i,j+1)/tem1
      t3=heta(i,j)/tem2
      t4=heta(i,j-1)/tem2
      con(i,j)=con(i,j)+(t1*qt(j+1)+t2*qt(j))*ae2(i,j+1)
     1        -(t4*qt(j)+t3*qt(j-1))*ae2(i,j)
  200 continue
      return

      entry extra2
c --- construct extra con for non-orthogonal coordinate system ---
c --- calculate b_no for velocity and transport correction equation.
     
      do 300 j=2,m2
      do 310 i=1,l1
      pt(i)=0.5*(fjeta(i,j)+fjeta(i,j+1))
      qt(i)=0.5*(rueta(i,j)+rueta(i,j+1))
  310 continue
      do 300 i=2,l2
      tem1=hksi(i+1,j)+hksi(i,j)
      tem2=hksi(i-1,j)+hksi(i,j)
      t1=hksi(i,j)/tem1
      t2=hksi(i+1,j)/tem1
      t3=hksi(i,j)/tem2
      t4=hksi(i-1,j)/tem2
      xje=t2*pt(i)+t1*pt(i+1)
      xme=t2*qt(i)+t1*qt(i+1)
      xjw=t4*pt(i)+t3*pt(i-1)
      xmw=t4*qt(i)+t3*qt(i-1)
      con(i,j)=con(i,j)+(xje-f(i,j,nf)*xme)*ak2(i+1,j)
     1        -(xjw-f(i,j,nf)*xmw)*ak2(i,j)
  300 continue

      do 400 i=2,l2
      do 410 j=1,m1
      pt(j)=0.5*(fjksi(i,j)+fjksi(i+1,j))
      qt(j)=0.5*(ruksi(i,j)+ruksi(i+1,j))
  410 continue
      do 400 j=2,m2
      tem1=heta(i,j+1)+heta(i,j)
      tem2=heta(i,j-1)+heta(i,j)
      t1=heta(i,j)/tem1
      t2=heta(i,j+1)/tem1
      t3=heta(i,j)/tem2
      t4=heta(i,j-1)/tem2
      xjn=t2*pt(j)+t1*pt(j+1)
      xmn=t2*qt(j)+t1*qt(j+1)
      xjs=t4*pt(j)+t3*pt(j-1)
      xms=t4*qt(j)+t3*qt(j-1)
      con(i,j)=con(i,j)+(xjn-f(i,j,nf)*xmn)*ae2(i,j+1)
     1        -(xjs-f(i,j,nf)*xms)*ae2(i,j)
  400 continue

      return

      end

      subroutine solve_c
	include 'sc.h'
	complex*16 denom,blcm,bl

      n=nf

	nsolve_c=50
	if(itest.eq.1)nsolve_c=20
      do 999 nt=1,nsolve_c
      write(*,*)'in loop',nt
c--------------------------------------------------------

c      go to 10
      ptc(1)=0.
      qtc(1)=0.
      do 11 i=2,l2
      bl=0.
      blp=0.
      blm=0.
      blcm=0.
      do 12 j=2,m2
      bl=bl+apc(i,j)
      if(j.ne.m2) bl=bl-ajp(i,j)
      if(j.ne.2) bl=bl-ajm(i,j)
      blp=blp+aip(i,j)
      blm=blm+aim(i,j)
      blcm=blcm+con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
     &    +ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)-apc(i,j)*fc(i,j)
   12 continue
      denom=bl-ptc(i-1)*blm
      if(abs(denom/bl).lt.1.0e-10) denom=1.e30
      ptc(i)=blp/denom
      qtc(i)=(blcm+blm*qtc(i-1))/denom
   11 continue
      bl=0.
      do 13 ii=2,l2
      i=l2+2-ii
      bl=bl*ptc(i)+qtc(i)
      do 13 j=2,m2
   13 fc(i,j)=fc(i,j)+bl
c-----------------------------------------------------------------
      ptc(1)=0.
      qtc(1)=0.
      do 21 j=2,m2
      bl=0.
      blp=0.
      blm=0.
      blcm=0.
      do 22 i=2,l2
      bl=bl+apc(i,j)
      if(i.ne.l2) bl=bl-aip(i,j)
      if(i.ne.2) bl=bl-aim(i,j)
      blp=blp+ajp(i,j)
      blm=blm+ajm(i,j)
      blcm=blcm+con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
     &    +ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)-apc(i,j)*fc(i,j)
   22 continue
      denom=bl-ptc(j-1)*blm
      if(abs(denom/bl).lt.1.0e-10) denom=1.e30
      ptc(j)=blp/denom
      qtc(j)=(blcm+blm*qtc(j-1))/denom
   21 continue
      bl=0.
      do 23 jj=2,m2
      j=m2+2-jj
      bl=bl*ptc(j)+qtc(j)
      do 23 i=2,l2
   23 fc(i,j)=fc(i,j)+bl
   10 continue

c	write(*,*)'now 10'
	fmax=1d-6
	bmax=0
	imax=1
	jmax=1
	do i=2,l2
	do j=2,m2
	if(id_bound(i,j).eq.1)then
      fij=abs(fc(i,j))
	fmax=dmax1(fmax,fij)
      blcm=(con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
     &  +ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)-apc(i,j)*fc(i,j))
     1   /apc(i,j)
	if(abs(blcm).gt.bmax)then
	bmax=abs(blcm)
	imax=i
	jmax=j
	endif
	endif
	enddo
	enddo
c	write(*,*)'l2=',l2,'m2=',m2
c	write(*,*)'bmax=',bmax,'imax=',imax,'jmax=',jmax,
c     1   'x=',x(imax,jmax),'y=',y(imax,jmax)
c        write(*,*)'solve nf,bmax,fmax,nt',nf,bmax,fmax,nt


	call solve_tdma_c
  999 continue
1000	continue

	return
	end


	subroutine solve_tdma_c
	include 'sc.h'
	complex*16 denom,temp,blcm
c -- TDMA solver ----------------------------
	blc=0
	do i=2,l2
	do j=2,m2
      blcm=con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
     &    +ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)-apc(i,j)*fc(i,j)
	if(blc.le.cdabs(blcm))blc=cdabs(blcm)
	enddo
	enddo
	write(*,*)'blc=',blc
      do 90 j=2,m2
      ptc(1)=0.0
      qtc(1)=fc(1,j)
      do 70 i=2,l2
      denom=apc(i,j)-ptc(i-1)*aim(i,j)
      ptc(i)=aip(i,j)/denom
      temp=con(i,j)+ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)
      qtc(i)=(temp+aim(i,j)*qtc(i-1))/denom
c	write(*,*)denom,temp
c	pause
   70 continue
      do 80 ii=2,l2
      i=l2+2-ii
   80 fc(i,j)=fc(i+1,j)*ptc(i)+qtc(i)
   90 continue
c---------------------------------------------------
      do 190 jj=2,m3
      j=m3+2-jj
      ptc(1)=0.0
      qtc(1)=fc(1,j)
      do 170 i=2,l2
      denom=apc(i,j)-ptc(i-1)*aim(i,j)
      ptc(i)=aip(i,j)/denom
      temp=con(i,j)+ajp(i,j)*fc(i,j+1)+ajm(i,j)*fc(i,j-1)
      qtc(i)=(temp+aim(i,j)*qtc(i-1))/denom
  170 continue
      do 180 ii=2,l2
      i=l2+2-ii
  180 fc(i,j)=fc(i+1,j)*ptc(i)+qtc(i)
  190 continue
c---------------------------------------------------
      do 290 i=2,l2
      ptc(1)=0.0
      qtc(1)=fc(i,1)
      do 270 j=2,m2
      denom=apc(i,j)-ptc(j-1)*ajm(i,j)
      ptc(j)=ajp(i,j)/denom
      temp=con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
      qtc(j)=(temp+ajm(i,j)*qtc(j-1))/denom
  270 continue
      do 280 jj=2,m2
      j=m2+2-jj
  280 fc(i,j)=fc(i,j+1)*ptc(j)+qtc(j)
  290 continue
c --------------------------------------------------
      do 390 ii=2,l3
      i=l3+2-ii
      ptc(1)=0.0
      qtc(1)=fc(i,1)
      do 370 j=2,m2
      denom=apc(i,j)-ptc(j-1)*ajm(i,j)
      ptc(j)=ajp(i,j)/denom
      temp=con(i,j)+aip(i,j)*fc(i+1,j)+aim(i,j)*fc(i-1,j)
      qtc(j)=(temp+ajm(i,j)*qtc(j-1))/denom
  370 continue
      do 380 jj=2,m2
      j=m2+2-jj
  380 fc(i,j)=fc(i,j+1)*ptc(j)+qtc(j)
  390 continue
	return
	end


	subroutine solve_tdma
	include 'sc.h'

c -- TDMA solver ----------------------------

      do 90 j=2,m2
	ifst1=0
	ilst=0
	do i=2,l2
	if(id_bound(i,j).ne.0)then
c	if(abs(ap(i,j)).ge.1d-8)then
	ifst1=i
	pt(ifst1-1)=0.0
	qt(ifst1-1)=f(ifst1-1,j,nf)
	goto 50
	endif
	enddo
	goto 90
50	do i=l2,2,-1
	if(id_bound(i,j).ne.0)then
	ilst=i
	goto 55
	endif
	enddo
	goto 90
55	continue
c	write(*,*)'ifst1,ilst=',ifst,ilst
	if(ifst1.eq.0.or.ilst.eq.0)goto 90
	do 70 i=ifst1,ilst
      denom=ap(i,j)-pt(i-1)*aim(i,j)
      pt(i)=aip(i,j)/denom
      temp=con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)
      qt(i)=(temp+aim(i,j)*qt(i-1))/denom
   70 continue
      do 80 i=ilst,ifst1,-1
   80 f(i,j,nf)=f(i+1,j,nf)*pt(i)+qt(i)

   90 continue
c---------------------------------------------------
      do 190 j=m3,2,-1

	ifst1=0
	ilst=0
	do i=2,l2
	if(id_bound(i,j).ne.0)then
	ifst1=i
	pt(ifst1-1)=0.0
	qt(ifst1-1)=f(ifst1-1,j,nf)
	goto 110
	endif
	enddo
	goto 190
110	do i=l2,2,-1
	if(id_bound(i,j).ne.0)then
	ilst=i
	goto 115
	endif
	enddo
	goto 190
115	if(ifst1.eq.0.or.ilst.eq.0)goto 190
c	write(*,*)'ifst1,ilst=',ifst,ilst
      pt(ifst1-1)=0.0
      qt(ifst1-1)=f(ifst1-1,j,nf)
      do 170 i=ifst1,ilst
      denom=ap(i,j)-pt(i-1)*aim(i,j)
      pt(i)=aip(i,j)/denom
      temp=con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)
      qt(i)=(temp+aim(i,j)*qt(i-1))/denom
  170 continue
      do 180 i=ilst,ifst1,-1
  180 f(i,j,nf)=f(i+1,j,nf)*pt(i)+qt(i)
  190 continue

	do 290 i=2,l2
	jfst=0
	jlst=0
	do j=2,m2
	if(id_bound(i,j).ne.0)then
	jfst=j
	pt(jfst-1)=0.0
	qt(jfst-1)=f(i,jfst-1,nf)
	goto 250
	endif
	enddo
	goto 290
250	do j=m2,2,-1
	if(id_bound(i,j).ne.0)then
	jlst=j
	goto 255
	endif
	enddo
	goto 290
255	if(jfst.eq.0.or.jlst.eq.0)goto 290
c	write(*,*)'jfst,jlst=',jfst,jlst
	do 270 j=jfst,jlst
      denom=ap(i,j)-pt(j-1)*ajm(i,j)
      pt(j)=ajp(i,j)/denom
      temp=con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf)
      qt(j)=(temp+ajm(i,j)*qt(j-1))/denom
  270 continue
      do 280 j=jlst,jfst,-1
  280 f(i,j,nf)=f(i,j+1,nf)*pt(j)+qt(j)
  290 continue

      do 390 i=l3,2,-1

	jfst=0
	jlst=0
	do j=2,m2
	if(id_bound(i,j).ne.0)then
	jfst=j
	pt(jfst-1)=0.0
	qt(jfst-1)=f(i,jfst-1,nf)
	goto 350
	endif
	enddo
	goto 390
350	do j=m2,2,-1
	if(id_bound(i,j).ne.0)then
	jlst=j
	goto 355
	endif
	enddo
	goto 390
355	if(jfst.eq.0.or.jlst.eq.0)goto 390
c	write(*,*)'jfst,jlst=',jfst,jlst
      pt(jfst-1)=0.0
      qt(jfst-1)=f(i,jfst-1,nf)
      do 370 j=jfst,jlst
      denom=ap(i,j)-pt(j-1)*ajm(i,j)
      pt(j)=ajp(i,j)/denom
      temp=con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf)
      qt(j)=(temp+ajm(i,j)*qt(j-1))/denom
  370 continue
      do 380 j=jlst,jfst,-1
  380 f(i,j,nf)=f(i,j+1,nf)*pt(j)+qt(j)
  390 continue
	return
	end

      subroutine solve
	include 'sc.h'

      n=nf

c      if(nf.ne.5) then
c       do j=1,m1
c          do i=1,l1
c	id=id_bound(i,j)
c	if(id.eq.0)then
c          ap(i,j)=1.d20
c          aip(i,j)=0.0
c          ajp(i,j)=0.0
c          aim(i,j)=0.0
c          ajm(i,j)=0.0
c          con(i,j)=0.0
c	f(i,j,nf)=0
c       endif
c          enddo
c       enddo
c      endif    

c	write(*,*)'in solve'
     
c       do j=1,m1
c	if(lst3.le.l1)then
c          ap(lst3,j)=1.d32
c          aip(lst3,j)=0.0
c          ajp(lst3,j)=0.0
c          aim(lst3,j)=0.0
c          ajm(lst3,j)=0.0
c          con(lst3,j)=0.0
c	endif
c       enddo


      do 999 nt=1,nsolve(nf)
	write(*,*)'in loop',nt,nf
c--------------------------------------------------------
c --- lblk = .true. for block-correction no good now turn off --------

	lblk(nf)=.true.
      if( .not.lblk(nf)) go to 10
      pt(1)=0.
      qt(1)=0.
      do 11 i=2,l2
	pt(i)=0
	qt(i)=0
      bl=0.
      blp=0.
      blm=0.
      blc=0.
      do 12 j=2,m2
	if(id_bound(i,j).eq.0)goto 12
      bl=bl+ap(i,j)
      if(j.ne.m2) bl=bl-ajp(i,j)
      if(j.ne.2) bl=bl-ajm(i,j)
      blp=blp+aip(i,j)
      blm=blm+aim(i,j)
      blc=blc+con(i,j)+aip(i,j)*f(i+1,j,n)+aim(i,j)*f(i-1,j,n)
     &    +ajp(i,j)*f(i,j+1,n)+ajm(i,j)*f(i,j-1,n)-ap(i,j)*f(i,j,n)
   12 continue
	if(bl.eq.0)goto 11
      denom=bl-pt(i-1)*blm
      if(abs(denom/bl).lt.1.0e-10) denom=1.e30
      pt(i)=blp/denom
      qt(i)=(blc+blm*qt(i-1))/denom
   11 continue
      bl=0.
      do 13 ii=2,l2
      i=l2+2-ii
      bl=bl*pt(i)+qt(i)
      do 13 j=2,m2
	if(id_bound(i,j).eq.0)goto 13
	f(i,j,n)=f(i,j,n)+bl
13	continue
c-----------------------------------------------------------------
      pt(1)=0.
      qt(1)=0.
      do 21 j=2,m2
	pt(j)=0
	qt(j)=0
      bl=0.
      blp=0.
      blm=0.
      blc=0.
      do 22 i=2,l2
	if(id_bound(i,j).eq.0)goto 22
      bl=bl+ap(i,j)
      if(i.ne.l2) bl=bl-aip(i,j)
      if(i.ne.2) bl=bl-aim(i,j)
      blp=blp+ajp(i,j)
      blm=blm+ajm(i,j)
      blc=blc+con(i,j)+aip(i,j)*f(i+1,j,n)+aim(i,j)*f(i-1,j,n)
     &    +ajp(i,j)*f(i,j+1,n)+ajm(i,j)*f(i,j-1,n)-ap(i,j)*f(i,j,n)
   22 continue
	if(bl.eq.0)goto 21
      denom=bl-pt(j-1)*blm
      if(abs(denom/bl).lt.1.0e-10) denom=1.e30
      pt(j)=blp/denom
      qt(j)=(blc+blm*qt(j-1))/denom
   21 continue
      bl=0.
      do 23 jj=2,m2
      j=m2+2-jj
      bl=bl*pt(j)+qt(j)
      do 23 i=2,l2
	if(id_bound(i,j).eq.0)goto 23
	f(i,j,n)=f(i,j,n)+bl
23	continue
   10 continue

c	write(*,*)'now 10'
	fmax=1d-6
	bmax=0
	imax=1
	jmax=1
	do i=2,l2
	do j=2,m2
	if(id_bound(i,j).ne.0)then
      fij=abs(f(i,j,nf))
	fmax=dmax1(fmax,fij)
      blc=(con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf)
     &  +ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)-ap(i,j)*f(i,j,nf))
     1   /ap(i,j)
	if(abs(blc).gt.bmax)then
	bmax=abs(blc)
	imax=i
	jmax=j
	endif
	endif
	enddo
	enddo
c	write(*,*)'l2=',l2,'m2=',m2
c	write(*,*)'bmax=',bmax,'imax=',imax,'jmax=',jmax,
c     1   'x=',x(imax,jmax),'y=',y(imax,jmax)
c        write(*,*)'solve nf,bmax,fmax,nt',nf,bmax,fmax,nt


	call solve_tdma
  999 continue
1000	continue
	return
	end

	subroutine user2
	include "sc.h"
      entry reset
      do 500 j=2,m2
      do 500 i=2,l2
      con(i,j)=0.
      ap(i,j)=0.
	aip(i,j)=0
	aim(i,j)=0
	ajp(i,j)=0
	ajm(i,j)=0
  500 continue
      return

	entry reset1
c called by us_elemag and test_us_elemag 
	do i=1,l1
	do j=1,m1
	ap(i,j)=0
	con(i,j)=0
	aip(i,j)=0
	aim(i,j)=0
	ajp(i,j)=0
	ajm(i,j)=0
	apc(i,j)=0
	enddo
	enddo
	return
	end
c-------------------------------------------   awrite   ------------

	subroutine awrite(iuv)
	include "sc.h"

c --- coefficient and boundary boundary temporary store and transfer
c --- for later use by aread and awrite.
      if(iuv.eq.2) go to 600

	do 501 j=1,m1
      biu(j)=gam(1,j)
  501 blu(j)=gam(l1,j)
      do 502 i=1,l1
      bju(i)=gam(i,1) 
  502 bmu(i)=gam(i,m1)
      do 503 i=1,l1
      do 503 j=1,m1
      conu(i,j)=con(i,j)
      cofu(i,j,1)=ap(i,j)
      cofu(i,j,2)=aip(i,j)
      cofu(i,j,3)=aim(i,j)
      cofu(i,j,4)=ajp(i,j)
      cofu(i,j,5)=ajm(i,j)
  503 continue
      return
  600 continue
      do 601 j=1,m1
      biv(j)=gam(1,j)
  601 blv(j)=gam(l1,j)
      do 602 i=1,l1
      bjv(i)=gam(i,1)
  602 bmv(i)=gam(i,m1)
      do 603 i=1,l1
      do 603 j=1,m1
      conv(i,j)=con(i,j)
      cofv(i,j,1)=ap(i,j)
      cofv(i,j,2)=aip(i,j)
      cofv(i,j,3)=aim(i,j)
      cofv(i,j,4)=ajp(i,j)
      cofv(i,j,5)=ajm(i,j)
  603 continue
      return
	end


      subroutine aread(iuv)
	include "sc.h"
      if(iuv.eq.2) go to 800

      do 701 j=1,m1
      gam(1,j)=biu(j)
  701 gam(l1,j)=blu(j)
      do 702 i=1,l1
      gam(i,1)=bju(i)
  702 gam(i,m1)=bmu(i)
      do 703 i=1,l1
      do 703 j=1,m1
      con(i,j)=conu(i,j)
      ap(i,j)=cofu(i,j,1)
      aip(i,j)=cofu(i,j,2)
      aim(i,j)=cofu(i,j,3)
      ajp(i,j)=cofu(i,j,4)
      ajm(i,j)=cofu(i,j,5)
  703 continue
      return
  800 continue
      do 801 j=1,m1
      gam(1,j)=biv(j)
  801 gam(l1,j)=blv(j)
      do 802 i=1,l1
      gam(i,1)=bjv(i)
  802 gam(i,m1)=bmv(i)
      do 803 i=1,l1
      do 803 j=1,m1
      con(i,j)=conv(i,j)
      ap(i,j)=cofv(i,j,1)
      aip(i,j)=cofv(i,j,2)
      aim(i,j)=cofv(i,j,3)
      ajp(i,j)=cofv(i,j,4)
      ajm(i,j)=cofv(i,j,5)
  803 continue
      return
      end


c -------------------------------------------------------------------
      subroutine stprint
c -------------------------------------------------------------------
	include "sc.h"
   10 format(54(1h*),3x,a10,3x,53(1h*),/)
   11 format(1h1,2x,54(1h*),3x,a10,3x,53(1h*),/)
   20 format(/,5h  i =,i6,12i10)
   30 format(2h j)
   35 format(1h2)
   40 format(i3,3x,13e10.3)
c -----------------------------------------------------------

      entry stream

c -- calculate the stream fnc. ---------
      
      f(2,2,10)=0.0
      do 82 i=2,lst3-1
      if(i.ne.2) t1=0.5*(ruksi(i-1,1)+ruksi(i,1))
      if(i.ne.2) f(i,2,10)=f(i-1,2,10)-rueta(i-1,2)*ae1(i-1,2)
     1                    -t1*ae2(i-1,2)
      do 82 j=3,m1-1
      t1=0.5*(rueta(i,j-1)+rueta(i,j))
      t2=0.5*(rueta(i-1,j-1)+rueta(i-1,j))
      t3=(hksi(i,j-1)*t2+hksi(i-1,j-1)*t1)/(hksi(i,j-1)+hksi(i-1,j-1))
      f(i,j,10)=f(i,j-1,10)+ruksi(i,j-1)*ak1(i,j-1)
     1         -t3*ak2(i,j-1)
   82 continue
      return

c ----------------------------------------------------------

      entry print

      if(lpgrid) go to 71

c -- extrapolate to get pressures on boundary --
      p(1,1)=p(2,1)+p(1,2)-p(2,2)
      p(l1,1)=p(l2,1)+p(l1,2)-p(l2,2)
      p(1,m1)=p(2,m1)+p(1,m2)-p(2,m2)
      p(l1,m1)=p(l2,m1)+p(l1,m2)-p(l2,m2)
      pref=p(ipref,jpref)
      do 70 j=1,m1
      do 70 i=1,l1
   70 p(i,j)=p(i,j)-pref

      nend=14
      nbeg=1

   71 continue
      if(.not.lpgrid) go to 72
      nbeg=15
      nend=16
   72 continue

      do 999 n=nbeg,nend
      nf=n
      if(.not.lpgrid) go to 55
      if(n.eq.15) print 10,titlex
      if(n.eq.16) print 10,titley
      go to 56
   55 continue
      if(.not.lprint(nf)) go to 999
      print 11,title(nf)
   56 continue
      ifst=1
      jfst=1
      ibeg=ifst-13
  110 continue
      ibeg=ibeg+13
      iend=ibeg+12
      iend=min0(iend,l1)
      print 20,(i,i=ibeg,iend)
      print 30
      jfl=jfst+m1
      do 115 jj=jfst,m1
      j=jfl-jj
      print 40,j,(f(i,j,nf),i=ibeg,iend)
  115 continue
      print 35
      if(iend.lt.l1) go to 110
  999 continue

      return
      end

c ------------------------------------------------------------------

      subroutine diflow(vv,i1,j1)

	include 'sc.h'
c      common/coef/flow,diff,acof
      data fl0/0.0/
c ischem :1-plds,2-central difference,3-upwind,4-hybrid,5-exp.
      ischem=1
c power law
        if((nf.eq.1.or.nf.eq.2))then
        flow=flow/vv/vv
        endif

c for temperature
	if(nf.eq.5)then
	flow=flow*cp(i1,j1)/vv
	endif
c for concentration
c molecular weight of SiC2: 0.0521
	if(nf.eq.6)then
        flow=flow/rho(i1,j1)*0.0521
	endif

      acof=diff
      if(diff.eq.fl0) return
      if(ischem.ne.1) go to 2
      pp=flow/diff
      p1=abs(pp)
      acof=dmax1(1.d-10,(1.0d0-0.1*p1)**5)
      acof=acof*diff
      if(flow.lt.fl0) acof=acof-flow
      return

    2 if(ischem.ne.2) go to 3
      temp=diff-dabs(flow)*0.5
      temp=temp/diff
      acof=diff*temp
      if(flow.lt.fl0) acof=acof-flow
      return
c
    3 if(ischem.ne.3) go to 4
      temp=diff
      temp=temp/diff
      acof=diff*temp
      if(flow.lt.fl0) acof=acof-flow
      return
c
    4 if(ischem.ne.4) go to 5
      temp=diff-dabs(flow)*0.5
      if(temp.le.fl0)return
      temp=temp/diff
      acof=diff*temp
      if(flow.lt.fl0) acof=acof-flow
      return
c
    5 if(ischem.ne.5) return
      temp=dabs(flow)/(exp(dabs(flow/diff))-1.)
      temp=temp/diff
      acof=diff*temp
      if(flow.lt.fl0) acof=acof-flow
      return
      end



      subroutine adaptive 
	include "sc.h"
      common/xyb/xbb(ni),ybb(ni)
      dimension pb1(ni,nj),pd(ni,nj),pc1(ni,nj),pc2(ni,nj)
      dimension pc3(ni,nj),pc3d(ni,nj)
      dimension dxc(ni,nj),dyc(ni,nj)

      fmax=0.0
      epsil=0.0001
      smooth=1.0
      orth=5.0
      cvn=20.0
      relgrd=0.5
      ngd=1

c non-dimensionless variables in MAGG

      smooth=1.0e-4
      orth=orth*1.0e4
      cvn=cvn*1.0e4
      do 30 i=1,l1
      do 30 j=1,m1
      dxc(i,j)=0.0
   30 dyc(i,j)=0.0

      do 40 i=2,l1
      do 40 j=2,m1
      if(i.le.lst) pc1(i,j)=0.9*((i-(2+lst)/2)/10.0)**2+0.1
      if(i.ge.lst) pc1(i,j)=0.8*((i-l1)/12.0)**2+0.2
      pb1(i,j)=pc1(i,j)
   40 pd(i,j)=1.0*pb1(i,j)
 
c----- find points to generate symmetry boundaries----
 
      ilast=ngd*150
      do 1000 nin=1,ilast
      if(mod(nin,50).eq.0)print 977,nin,fmax,xc(20,20)
     1 ,yc(20,20),xc(3,m1),yc(3,m1)
  977 format(1x,'MAGG procedure',3x,i5,2x,e9.3,2x,4e9.3)
      ssmot=0.
      ssweg=0.

      do 200 j=3,m2
      do 200 i=3,l2
c --- i=2,l1  if it is a line of symmetry otherwise i=3,l1 -----
      xt=0.5*(xc(i+1,j)-xc(i-1,j))
      xn=0.5*(xc(i,j+1)-xc(i,j-1))
      yt=0.5*(yc(i+1,j)-yc(i-1,j))
      yn=0.5*(yc(i,j+1)-yc(i,j-1))
      xts=xt*xt 
      xns=xn*xn 
      yts=yt*yt 
      yns=yn*yn 
      xtyn=xt*yn
      xtyt=xt*yt
      xnyn=xn*yn
      xnyt=xn*yt
      xtxn=xt*xn
      ytyn=yt*yn
      cacbij=xtyn-xnyt
      rjcbij=1.0/cacbij 
      rjcb3=2.0*rjcbij**3 
      xtt=xc(i+1,j)-2.0*xc(i,j)+xc(i-1,j) 
      xnn=xc(i,j+1)-2.0*xc(i,j)+xc(i,j-1) 
      ytt=yc(i+1,j)-2.0*yc(i,j)+yc(i-1,j) 
      ynn=yc(i,j+1)-2.0*yc(i,j)+yc(i,j-1) 
      xtn=0.25*(xc(i+1,j+1)+xc(i-1,j-1)-xc(i-1,j+1)-xc(i+1,j-1))
      ytn=0.25*(yc(i+1,j+1)+yc(i-1,j-1)-yc(i-1,j+1)-yc(i+1,j-1))
      aa=xtyt+xnyn
      bb=yts+yns
      cc=xts+xns
      alf=rjcb3*(xns+yns) 
      bt=rjcb3*(xtxn+ytyn)
      gama=rjcb3*(xts+yts)
 
      a1=-aa*alf*smooth
      a2=2.0*aa*bt*smooth
      a3=-aa*gama*smooth
      b1=bb*alf*smooth
      b2=-2.0*bb*bt*smooth
      b3=bb*gama*smooth
      c1=cc*alf*smooth
      c2=-2.0*cc*bt*smooth
      c3=cc*gama*smooth
 
      othij=orth*pd(i,j)
      rrjc3=1/rjcb3 
      bt=bt*rrjc3 
      a1=a1+othij*xn*yn 
      a2=a2+othij*(xt*yn+xn*yt) 
      a3=a3+othij*xt*yt 
      b1=b1+othij*xn*xn 
      b2=b2+othij*(bt+xt*xn)*2.0
      b3=b3+othij*xt*xt 
      c1=c1+othij*yn*yn 
      c2=c2+othij*(bt+yn*yt)*2.0
      c3=c3+othij*yt*yt 
      gt=0.5*(pd(i+1,j)-pd(i-1,j))
      gn=0.5*(pd(i,j+1)-pd(i,j-1))
      gx=orth*0.5*(gt*yn-gn*yt)*bt**2*rjcbij
      gy=orth*0.5*(gn*xt-gt*xn)*bt**2*rjcbij
 
      pij=2.0*cvn*pb1(i,j) 
      a1=a1 -xnyn*pij 
      a2=a2+pij*(xtyn+xnyt) 
      a3=a3 -xtyt*pij 
      b1=b1+yns*pij 
      b2=b2 -2.0*ytyn*pij 
      b3=b3+yts*pij 
      c1=c1+xns*pij 
      c2=c2 -2.0*xtxn*pij 
      c3=c3+xts*pij 
      wt=0.5*(pb1(i+1,j)-pb1(i-1,j))
      wn=0.5*(pb1(i,j+1)-pb1(i,j-1))
      cff=cacbij*cvn
      wx=cff*(yn*wt-yt*wn)
      wy=cff*(xt*wn-xn*wt)
      res1=b1*xtt+b2*xtn+b3*xnn+a1*ytt+a2*ytn+a3*ynn+gx+wx
      res2=a1*xtt+a2*xtn+a3*xnn+c1*ytt+c2*ytn+c3*ynn+gy+wy
 
      a11=-2.0*(b1+b3)
      a12=-2.0*(a1+a3)
      a22=-2.0*(c1+c3)
      a21=a12 
 
      det=a11*a22-a12*a21 
      rdet=1./det 
 
      dx=(res2*a21-res1*a22)*rdet 
      dy=(res1*a12-res2*a11)*rdet 
 
c---------- limit the mesh excursions using projections---------
 
      ratio=.25d0 
      xij=xc(i,j) 
      yij=yc(i,j) 
 
      xdisp=xc(i+1,j)-xij 
      ydisp=yc(i+1,j)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i+1,j+1)-xij 
      ydisp=yc(i+1,j+1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i,j+1)-xij 
      ydisp=yc(i,j+1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i-1,j+1)-xij 
      ydisp=yc(i-1,j+1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i-1,j)-xij 
      ydisp=yc(i-1,j)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i-1,j-1)-xij 
      ydisp=yc(i-1,j-1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i,j-1)-xij 
      ydisp=yc(i,j-1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      xdisp=xc(i+1,j-1)-xij 
      ydisp=yc(i+1,j-1)-yij 
      vmag2=xdisp*xdisp+ydisp*ydisp 
      dot=dx*xdisp+dy*ydisp 
      ratio=dmax1(dot/vmag2,ratio)
 
      scale=dmin1(relgrd,.25d0/ratio) 
c note : rho & gam were used as excursion arrays. 
      dxc(i,j)=scale*dx 
      dyc(i,j)=scale*dy 
 
      smot=(xts+xns+yts+yns)/cacbij 
      ssmot=ssmot+smot
      cacb2=cacbij*cacbij 
      weigt=pb1(i,j)*cacb2*cvn 
      ssweg=ssweg+weigt 
  200 continue
 
c------find the change in boundary points---- 
 
      do 220 j=3,m2 
      do 220 i=3,l2 
      aaa=1.0
      if(i.eq.lst3)aaa=0.0
      xc(i,j)=xc(i,j)+dxc(i,j)*aaa
      yc(i,j)=yc(i,j)+dyc(i,j)
  220 continue

      do 222 j=2,m2
      yc(2,j)=yc(3,j)
  222 yc(l1,j)=yc(l2,j)
      do 223 i=3,l1 
      xc(i,2)=xc(i,3)
  223 xc(i,m1)=xc(i,m2)
 
      fmax=1.e-8
      do 400 j=3,m2 
      do 400 i=3,l2 
      fij1=abs(dxc(i,j))
      fij2=abs(dyc(i,j))
  400 fmax=dmax1(fmax,fij1,fij2)
      if(fmax.le.epsil) goto 1001
 1000 continue
 1001 continue
      return
      end 


	real*8 function rinterp0(xxx,yyy,m1,n1,x,y,u)
	implicit real*8 (a-h,o-z)
	dimension x(m1,n1),y(m1,n1),u(m1,n1)
	do 100 i1=1,m1
	do 100 j1=1,n1
	if(xxx.ge.x(i1,j1).and.xxx.le.x(i1+1,j1+1)
     1   .and.yyy.ge.y(i1,j1).and.yyy.le.y(i1+1,j1+1))then
	goto 200
	endif
100	continue
	stop 'err in rinterp'
200	delx=x(i1+1,j1+1)-x(i1,j1)
	dely=y(i1+1,j1+1)-y(i1,j1)
	x1=(xxx-x(i1,j1))/delx
	y1=(yyy-y(i1,j1))/dely
	uuu=u(i1,j1)*(1-x1)*(1-y1)
     1     +u(i1+1,j1)*x1*(1-y1)
     1     +u(i1,j1+1)*(1-x1)*y1
     1     +u(i1+1,j1+1)*x1*y1
	rinterp0=uuu
	return
	end

      subroutine lineinp(xnew,ynew)
	include "sc.h"
      dimension xl(ni),yl(ni),tt(ni),tm(ni)
	do i=1,ni
	xl(i)=0
	yl(i)=0
	tt(i)=0
	tm(i)=0
	enddo

      l=inum+2
      xl(1)=xin(1)+(xin(1)-xin(3))
      xl(2)=xin(2)+(xin(1)-xin(3))
      do 10 i=3,l
      xl(i)=xin(i-2)
   10 yl(i)=yin(i-2)
      xl(l+1)=xl(l-1)+(xl(l)-xl(l-2))
      xl(l+2)=xl(l)+(xl(l)-xl(l-2))
      shift=(yl(3)-yl(4))/(xl(3)-xl(4))-(yl(4)-yl(5))/(xl(4)-xl(5))
      yl(2)=yl(3)+(xl(2)-xl(3))*(shift+(yl(3)-yl(4))/(xl(3)-xl(4)))
      yl(1)=yl(2)+(xl(1)-xl(2))*(shift+(yl(2)-yl(3))/(xl(2)-xl(3)))
      shift=(yl(l)-yl(l-1))/(xl(l)-xl(l-1))-(yl(l-1)-yl(l-2))
     1 /(xl(l-1)-xl(l-2))
      yl(l+1)=yl(l)+(xl(l+1)-xl(l))*(shift+(yl(l)-yl(l-1))
     1 /(xl(l)-xl(l-1)))
      yl(l+2)=yl(l+1)+(xl(l+2)-xl(l+1))*(shift+(yl(l+1)-yl(l))
     1 /(xl(l+1)-xl(l)))
      do 20 i=1,l+1
   20 tm(i)=(yl(i+1)-yl(i))/(xl(i+1)-xl(i))
      do 30 i=3,l
      dome=abs(tm(i+1)-tm(i))+abs(tm(i-1)+tm(i-2))
      if(dome.lt.1.e-10) then
      tt(i)=0.0
      else
      tt(i)=(abs(tm(i+1)-tm(i))*tm(i-1)+abs(tm(i-1)-tm(i-2))*tm(i))/dome
      endif
   30 continue

      xmin=1.e9
      do 40 i=3,l
      if(abs(xl(i)-xnew).lt.xmin) then
      ist=i
      xmin=abs(xl(i)-xnew)
      endif
   40 continue

      if((xl(ist)-xnew).gt.0) then
      x1=xl(ist-1)
      y1=yl(ist-1)
      x2=xl(ist)
      y2=yl(ist)
      t1=tt(ist-1)
      t2=tt(ist)
      else
      x1=xl(ist)
      y1=yl(ist)
      x2=xl(ist+1)
      y2=yl(ist+1)
      t1=tt(ist)
      t2=tt(ist+1)
      endif
      p0=y1
      p1=t1
      p2=(3.*(y2-y1)/(x2-x1)-2.*t1-t2)/(x2-x1)
      p3=(t1+t2-2.*(y2-y1)/(x2-x1))/(x2-x1)**2
      ynew=p0+p1*(xnew-x1)+p2*(xnew-x1)**2+p3*(xnew-x1)**3
      dydx=p1+2.0*p2*(xnew-x1)+3.0*p3*(xnew-x1)**2
      dy2dx2=2.0*p2+6.0*p3*(xnew-x1)
      rrnew=dy2dx2/(1.0+dydx**2)**1.5
      return
      end
