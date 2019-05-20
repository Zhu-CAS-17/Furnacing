c these subroutines are copied from sc_new.f
c in the future we may delete the subroutines in the sc_new.f
c and use these instead
c the program sc_advan.f is also using these subroutines
c Nov. 11, 1999
c cal_ex1 is used by sc_new.f
c cal_gr_rate is used by sc_advan.f
c cal_d0 estimates the coefficient used to calculate the diffusivity
c taken as 1d-5 for high temperature and 0.5d-5 for low temperature.

c calculate growth rate by sticking coefficient model
c      call sticking

c      stop
c      end



        block data
        include 'sc.h'

      data title(1),title(2),title(3),titlex,titley
     1    /' vel u ',' vel v ','presure','   xc  ','   yc  '/
      data title(5),title(6),title(10)
     1    /' temp  ','rot. v.','stream '/
      data title(11),title(12)/' rho  ','  gam  '/
      data nfmax,np,npc/6,3,4/
      data lstop,lsolve,lprint/25*.false./
      data lblk/10*.false./
      data iter,ndt/0,0/
      data relax,ntimes,nsolve/14*1.0,10*1,10*1/
      data ipref,jpref,rhocon/11,30,1.0/
      data cfo,dtm,epsil,tinit/0.0001,5.,.0001,0.0/
	data id_graphite,id_soft,id_sic,id_argon,id_hard,id_foam,id_max
     1   /1,2,3,4,5,6,6/
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: griphe foam

        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi

c at temperatures above 2900 K, the vapor pressure of silicon is less
c then that of sic2. The following uses vapor pressure of sic2 at
c temperature below 2900 K and that of silicon above 2900 K.
c from Lilov
      data p_a/5.131d-8,1.372d-6,2.473d-5,3.207d-4,
     1    3.153d-3,2.452d-2,1.559d-1,8.352d-1,
     1   3.839,15.47,55.56,96.492,170.02,456.3,1137.9,2366,
     1   4419,7922,10477/

      data t_a/1500,1600,1700,1800,1900,2000,2100,2200,
     1  2300,2400,2500,2546,2600,2700,2800,2900,
     1   3000,3100,3150/

c from Lilov
      data psic/6.2d-13,3.2d-11,1.0d-9,2.3d-8,3.7d-7,4.4d-6,
     1  4.2d-5,3.2d-4,2.0d-3,1.11d-2,5.23d-2,1.02d-1,0.2188,
     1  0.8195,2.784,8.673,24.98,67.19,107.61/

c vapor pressure of sic2 for reaction 2sic=sic2+si
c from Lilov
       data psic2/5.131d-8,1.372d-6,2.473d-5,3.207d-4,
     1    3.153d-3,2.452d-2,1.559d-1,8.352d-1,
     1   3.839,15.47,55.56,96.492,170.02,456.3,1137.9,2659,
     1   5855,12250,17418/

c for t ge 2900 si vapor pressures are taken as satuated silicon pressure
c times activity of silicon.
c from Lilov
      data psi/3.00d-6,4.23d-5,4.35d-4,3.451d-3,2.191d-2,1.152d-1,
     1  0.5163,2.012,6.96,21.67,61.44,96.522,170.02,456.3,
     1  1137.9,2366,4419,7922,10477/

c from Lilov
      data psi2c/2.52d-8,6.85d-7,1.24d-5,1.63d-4,1.61d-3,
     1  1.256d-2,8.006d-2,4.281d-1,1.966,7.913,28.32,49.12,
     1   91.48,269.8,733.6,1853,4384,9811,14398/

      data t_a1/1500,1600,1700,1800,1900,2000,2100,2200,
     1  2300,2400,2500,2600,2700,2800,2900,
     1   3000,3100,3200/

c by self calculation, may have errors
c      data p_a/1.346d-8,4.001d-7,7.926d-6,1.119d-4,
c     1   1.189d-3,0.00992,0.0597,0.3815,
c     1   1.855,7.875, 29.65, 62.241, 100.5, 310.2, 908.0,2786,
c     1   4861, 8177, 12780/
c from self calculation
c        data psic/.6438E-12,.3337E-10,.1081E-08,.2368E-07,.3731E-06,
c     1 .4443E-05,.4163E-04,.3171E-03,.2018E-02,.1097E-01,.5195E-01,
c     1 .2177E+00,.8180E+00,.2791E+01,.8728E+01,.2524E+02,.6806E+02,
c     1 .1721E+03/
c from self calculation
c        data psic2/.1346E-07,.4001E-06,.7926E-05,.1119E-03,.1189E-02,
c     1 .9920E-02,.5970E-01,.3815E+00,.1855E+01,.7875E+01,.2965E+02,
c     1 .1005E+03,.3102E+03,.9080E+03,.2798E+04,.7980E+04,.2122E+05,
c     1 .5516E+05/
c from self calculation
c        data psi/.7681E-05,.1026E-03, .1007E-02, .7657E-02, .4688E-01,
c     1  .2390E+00,.1174E+01,.3967E+01, .1340E+02, .4089E+02, .1139E+03,
c     1  .2932E+03,.7025E+03,.1533E+04, .2786E+04, .4861E+04, .8177E+04,
c     1  .1278E+05/
c from self calculation
c        data psi2c/.1164E-07,.3256E-06,.6077E-05,.8138E-04,
c     1 .8230E-03,.6565E-02,.4817E-01,.2335E+00,.1093E+01,
c     1 .4485E+01,.1637E+02,.5390E+02,.1619E+03,.4482E+03,
c     1 .1154E+04,.2780E+04,.6338E+04,.1374E+05/


        end


	function cal_parti(t0)
	implicit real*8 (a-h,o-z)

        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi


	n=19
      pp=rinterp(t_a,p_a,n,t0)
	goto 100
c calculate vapor pressure of sic2 for reaction 2sic=sic2+si
        if(t0.lt.1500)then
        pp=5.131d-8
        goto 100
        else if(t0.gt.3150)then
        pp=17418
        goto 100
        endif
        if(t0.ge.1500.and.t0.le.2000)then
        a1=-34075.8
        b1=15.4274
        else if(t0.ge.2000.and.t0.le.2546)then
        a1=-33526.61
        b1=15.1528
        else if(t0.ge.2546.and.t0.le.2900)then
        a1=-33012.65
        b1=14.9510
        else if(t0.ge.2900.and.t0.le.3150)then
        a1=-32734.33
        b1=14.8550
	endif
	pp=10**(a1/t0+b1)
100	cal_parti=pp
	return
	end

      subroutine cal_gr_rate
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: graphite foam 
c i_rt, i_rb, j_rr, pres_sys

	include 'sc.h'
c partial pressure of SiC2
	common /cal_d0_v/cal_d0_v
	cal_d0_v=1d-5

	do i=1,l1
	do j=1,m1

	id=id_bound(i,j)

c	if(id.eq.0)then
c	do nf=1,4
c	f(i,j,nf)=0
c	enddo
c	f(i,j,6)=0
c	endif

	id=id_ele(i,j)
	parti_p(i,j)=0
	ex1(i,j)=0
	scargon(i,j)=0
	rho_s(i,j)=rho_cal(i,j)

	rho(i,j)=0.44
c	if(id.eq.3.or.id.eq.4)rho(i,j)=ccen(i,j)*0.040097
c     1    +(200/8.31/t(i,j)-ccen(i,j))*0.039948

c id 4 and 14 all for argon gas
        if(id.ne.3.and.id.ne.4.and.id.ne.14)then
	goto 200
	endif

c	goto 200

	t0=t(i,j)
c      pp=cal_parti(t0)
c using spline interpolation
      pp=cal_parti_s(t0)


c	fo(i,j,1)=0
c	fo(i,j,2)=0

	ccen(i,j)=pp/8.31/t(i,j)
	fo(i,j,6)=ccen(i,j)


100	parti_p(i,j)=pp

	if(id.eq.3)then
c powder diameter 1600 micro meter=0.0016m
c molecular weight 0.040097 kg/mol
c sour: mole flux of dissociation of SiC
c ex1: mol/s
c see formulation in the powder charge
	pp=8.31*t(i,j)*ccen(i,j)
	delp=pp-parti_p(i,j)
      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))

c     sour=6/0.0016*xapa*
c     1    (parti_p(i,j)-pp)*vol(i,j)
c  produce of SiC2, kg/m3/t
      sour=6/0.0016*xapa*
     1    parti_p(i,j)*0.0521*vol(i,j)
c	if(j.eq.j_rr)sour=0
c scargon of SiC2: mol/s
	scargon(i,j)=sour
c ex1 of total mass: kg/s
      ex1(i,j)=6/0.0016*xapa*(parti_p(i,j)-pp)
     1      *vol(i,j)*0.0521
	endif

200	continue
	enddo
	enddo

	if(iter.eq.0.and.mod(ndt,10).eq.0)then
	open(1,file='gr_rate_r')
	endif

	
	mm=j_rr
	do j=2,mm-1

c top

	i=i_rt
	pp=8.31*t(i,j)*ccen(i,j)
c growthsp: m/s
c molecular weight 0.040097 kg/mol
c density: 3160 kg/m3
c ccen : mol/m3
	growthsp(j)=0
	if(ccen(i,j).eq.0)then
	growthsp(j)=0
	else

c	fo(i,j,6)=200/8.31/t(i,j)
c	fo(i,j,1)=0
c	fo(i,j,2)=0
c      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))
c     growthsp(j)=xapa*
c     1    (pp)*0.016/(parti_p(i,j)/8.31/t(i,j))
c use striking model
c     growthsp(j)=xapa*
c     1    0.016*parti_p(i,j)/(parti_p(i,j)/8.31/t(i,j))
c divided by local species concentration (mol/m3) to get species speed m/s
c     growthsp(j)=xapa*
c     1    (pref-parti_p(i,j))/(parti_p(i,j)/8.31/t(i,j))

	pa0=parti_p(i_rb,j_rr)
	ca0=pa0/8.31/t(i_rb,j_rr)
      ta0=t(i_rb,j_rr)

c z_pres is set in radi_ini
      i_rb1=interid(z_pres,x,l1)
      i_rb2=interid(z_pres2,x,l1)
      pa0=parti_p(i_rb1,j_rr)
      ta0=t(i_rb1,j_rr)
      write(*,*)'pa0 ',pa0
      do i_rb12=i_rb1,i_rb2
      pa0_2=parti_p(i_rb12,j_rr)
      if(pa0_2.gt.pa0)then
      pa0=pa0_2
      ta0=t(i_rb12,j_rr)
      endif
      enddo
      write(*,*)'pa0 changed to ',pa0,'after search for maximum ta0'


      pres_sys_m=pres_sys
      total_vapor=cal_parti_sic(ta0)+cal_parti_si(ta0)
     1  +cal_parti_sic2(ta0)+cal_parti_si2c(ta0)

      if(total_vapor.gt.pres_sys)then
      pres_sys_m=total_vapor
      endif


      c_all=pres_sys_m/8.31/t(i_rb,j_rr)
c pres_sys: system pressure

	if(iter.eq.0)then
	write(*,*)'pa0=',pa0,'c_all',c_all
	endif

	uref=0

10	uref1=uref
c ca1 concentration of A at top
c ca0 concentration of A at bottom
c assume c_all and ca0 are fixed
c x=0.606/sqrt(t(i,j))

      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))

      pa1=parti_p(i,j)+uref*pres_sys_m/2/8.31/t(i,j)/xapa
      cred=(pres_sys_m-2*pa1)/(pres_sys_m-2*pa0)
c zt0, zb0 from radi_ini
      d_ab=cal_d0()*(t(i,j)/273)**1.8*(1d5/pres_sys_m)

      uref=log(cred)*d_ab/(zt0-zb0)
	if(iter.eq.0)then
	write(*,*)'uref',uref,'ca0',ca0,'ca1',ca1,'pa0',pa0
	endif

	if(abs(uref-uref1).ge.1d-5)goto 10

c growthsp: m/s
c ccen: mol/m3
c molecular weight 0.040097 kg/mol
c density: 3160 kg/m3
	growthsp(j)=uref
c	write(*,*)pp,'pascal',growthsp(j)
	if(iter.eq.0)then
	write(*,*)pp,'pascal',growthsp(j),'m/s',ccen(i,j),'mol/m3',t(i,j)
	endif

	endif
c	scargon(i,j)=-growthsp(j)*ccen(i,j)*ak1(i,j)
c 1 mol of SiC2 corresponds to 2 mol SiC
      ghk=uref*pres_sys_m/8.31/t(i,j)*0.040097/3160*3600*1000

	if(iter.eq.0)then
	write(*,*)'growth rate',ghk,'mm/hour',pp,'Pa'
	endif

	if(iter.eq.0.and.mod(ndt,10).eq.0)then
      write(1,'(5f16.6)')r(i,j),ghk,t(i,j),total_vapor,pres_sys
	endif

	enddo
	if(iter.eq.0.and.mod(ndt,10).eq.0)then
	close(1)
	endif

	return
	end

      function cal_d0()
	implicit real*8 (a-h,o-z)
	common /cal_d0_v/cal_d0_v
c for low pressure
c	cal_d0=5d-6

c for high pressure as used by Sterling
c      cal_d0=0.1d-4
	cal_d0=cal_d0_v
      return
      end

c this subroutine is used in sc_new.f
	subroutine cal_ex1
c 1:graph 2:soft felt 3: sic powder 4: Argon 5: hard felt 6: graphite foam 
	include 'sc.h'
c partial pressure of SiC2
	common /cal_d0_v/cal_d0_v
	cal_d0_v=1d-5

	do i=1,l1
	do j=1,m1

	id=id_bound(i,j)

c	if(id.eq.0)then
c	do nf=1,4
c	f(i,j,nf)=0
c	enddo
c	f(i,j,6)=0
c	endif

	id=id_ele(i,j)
	parti_p(i,j)=0
	ex1(i,j)=0
	scargon(i,j)=0
	rho_s(i,j)=rho_cal(i,j)

	rho(i,j)=0.44
c	if(id.eq.3.or.id.eq.4)rho(i,j)=ccen(i,j)*0.040097
c     1    +(200/8.31/t(i,j)-ccen(i,j))*0.039948

c id 4 and 14 all for argon gas
        if(id.ne.3.and.id.ne.4.and.id.ne.14)then
	goto 200
	endif

c	goto 200

	t0=t(i,j)
c      pp=cal_parti(t0)
c using spline interpolation
      pp=cal_parti_s(t0)


c	fo(i,j,1)=0
c	fo(i,j,2)=0

	ccen(i,j)=pp/8.31/t(i,j)
	fo(i,j,6)=ccen(i,j)
	goto 100


100	parti_p(i,j)=pp

	if(id.eq.3)then
c powder diameter 1600 micro meter=0.0016m
c molecular weight 0.040097 kg/mol
c sour: mole flux of dissociation of SiC
c ex1: mol/s
c see formulation in the powder charge
	pp=8.31*t(i,j)*ccen(i,j)
	delp=pp-parti_p(i,j)

      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))

c     sour=6/0.0016*xapa*
c     1    (parti_p(i,j)-pp)*vol(i,j)
c  produce of SiC2, kg/m3/t
      sour=6/0.0016*xapa*
     1    parti_p(i,j)*0.0521*vol(i,j)
c	if(j.eq.j_rr)sour=0
c scargon of SiC2: mol/s
	scargon(i,j)=sour
c ex1 of total mass: kg/s

      ex1(i,j)=6/0.0016*xapa*(parti_p(i,j)-pp)
     1      *vol(i,j)*0.0521
	endif

200	continue
	enddo
	enddo


	
	if(iter.eq.0.and.mod(ndt,10).eq.0)then
	open(1,file='gr_rate_r')
	endif

	mm=j_rr

        do j=1,mm-1
c top
	i=i_rt
	pp=8.31*t(i,j)*ccen(i,j)
c growthsp: m/s
c molecular weight 0.040097 kg/mol
c density: 3160 kg/m3
c ccen : mol/m3
	growthsp(j)=0
	if(ccen(i,j).eq.0)then
	growthsp(j)=0
	else

c	fo(i,j,6)=200/8.31/t(i,j)
c	fo(i,j,1)=0
c	fo(i,j,2)=0
c      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))
c     growthsp(j)=xapa*
c     1    (pp)*0.016/(parti_p(i,j)/8.31/t(i,j))
c use striking model
c     growthsp(j)=xapa*
c     1    0.016*parti_p(i,j)/(parti_p(i,j)/8.31/t(i,j))
c divided by local species concentration (mol/m3) to get species speed m/s
c     growthsp(j)=xapa*
c     1    (pref-parti_p(i,j))/(parti_p(i,j)/8.31/t(i,j))

c pa0 needs to be replaced by the vapor pressure in the middle of powder


c z_pres is set in radi_ini
      i_rb1=interid(z_pres,x,l1)
      i_rb2=interid(z_pres2,x,l1)


      vol_powder=0
      i_totalp=0
      do j1=1,j_rr
      do i1=i_rb1,i_rb2
      if(i_powder(i1,j1).eq.1)then
      vol_powder=vol_powder+vol(i1,j1)
      i_totalp=i_totalp+1
      endif
      enddo
      enddo

      write(*,*)'total powder element',i_totalp
50    continue

      if(vol_powder.gt.crys_len*rc*rc*pi*2.5)goto 80
      tmax1=0
      i1=0
      j1=0
      do j1=1,j_rr-1
      do i1=i_rb1+1,i_rb2-1
      if(t(i1,j1).gt.tmax1.and.i_powder(i1,j1).eq.0)then
      tmax1=t(i1,j1)
      imaxp=i1
      jmaxp=j1
      endif
      enddo
      enddo
      write(*,*)'tmax1',tmax1,imaxp,jmaxp,i_rb1,i_rb2,j_rr
      if(i1.eq.0.or.j1.eq.0)goto 80
c      i_powder(imaxp,jmaxp)=1
c      vol_powder=vol_powder+vol(imaxp,jmaxp)
c      i_totalp=i_totalp+1
c      goto 50
80    continue

      write(*,*)'total powder element',i_totalp

      write(*,*)'powder size',vol_powder,'imaxp',imaxp,jmaxp
      write(*,*)'crystal size',crys_len*rc*rc*pi
      pa0=parti_p(imaxp,jmaxp)
      ta0=t(imaxp,jmaxp)

c      pa0=parti_p(i_rb1,j_rr)
c      ta0=t(i_rb1,j_rr)
c      write(*,*)'pa0 ',pa0
c      do i_rb12=i_rb1,i_rb2
c      pa0_2=parti_p(i_rb12,j_rr)
c      if(pa0_2.gt.pa0)then
c      pa0=pa0_2
c      ta0=t(i_rb12,j_rr)
c      endif
c      enddo

      write(*,*)'pa0 changed to ',pa0,'after search for maximum ta0'
      write(*,*)'ta0', ta0

c ta0 only used for calculating total pressure
      pres_sys_m=pres_sys
      total_vapor=cal_parti_sic(ta0)+cal_parti_si(ta0)
     1  +cal_parti_sic2(ta0)+cal_parti_si2c(ta0)

      if(total_vapor.gt.pres_sys)then
      pres_sys_m=total_vapor
      endif


	ca0=pa0/8.31/t(i_rb,j_rr)
      c_all=pres_sys_m/8.31/t(i_rb,j_rr)
c pres_sys: system pressure
	if(iter.eq.0)then
	write(*,*)'pa0=',pa0,'c_all',c_all
	endif

	uref=0

10	uref1=uref

c ca1 concentration of A at top
c ca0 concentration of A at bottom
c assume c_all and ca0 are fixed
c x=0.606/sqrt(t(i,j))

      xapa=1/sqrt(2*3.1415*0.040097*8.31*t(i,j))

      pa1=parti_p(i,j)+uref*pres_sys_m/2/8.31/t(i,j)/xapa
      cred=(pres_sys_m-2*pa1)/(pres_sys_m-2*pa0)

c zt0, zb0 from radi_ini
      d_ab=cal_d0()*(t(i,j)/273)**1.8*(1d5/pres_sys_m)

      uref=log(cred)*d_ab/(zt0-zb0)
	if(iter.eq.0)then
	write(*,*)'uref',uref,'ca0',ca0,'ca1',ca1,'pa0',pa0
	endif

	if(abs(uref-uref1).ge.1d-5)goto 10

c growthsp: m/s
c ccen: mol/m3
c molecular weight 0.040097 kg/mol
c density: 3160 kg/m3
      growthsp(j)=uref


c	write(*,*)pp,'pascal',growthsp(j)
c	growthsp(j)=-0.009*(ccen(i,j)-ccen(i-1,j))/(x(i,j)-x(i-1,j))
c     1    /(200/8.31/t(i,j)-ccen(i,j))
	if(iter.eq.0)then
	write(*,*)pp,'pascal',growthsp(j),'m/s',ccen(i,j),'mol/m3',t(i,j)
	endif

	endif
c	scargon(i,j)=-growthsp(j)*ccen(i,j)*ak1(i,j)
c 1 mol of SiC2 corresponds to 2 mol SiC
c      ghk=2*growthsp(j)*ca1*0.040097/3160*3600*1000
      ghk=uref*pres_sys_m/8.31/t(i,j)*0.040097/3160*3600*1000

	if(iter.eq.0)then
	write(*,*)'growth rate',ghk,'mm/hour',pp,'Pa'
	endif

	if(iter.eq.0.and.mod(ndt,10).eq.0)then
      write(1,'(5f16.6)')r(i,j),ghk,t(i,j),total_vapor,pres_sys
	endif

	enddo
	if(iter.eq.0.and.mod(ndt,10).eq.0)then
	close(1)
	endif
	return
	end


        function cal_parti_sic(t0)
	implicit real*8 (a-h,o-z)
        real*8 y2(19)

        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi

        n=19
        call spline(t_a,psic,19,1.d30,1.d30,y2)
        call splint(t_a,psic,y2,19,t0,pp)
100     cal_parti_sic=pp
	return
	end

        function cal_parti_sic2(t0)
	implicit real*8 (a-h,o-z)

        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi

        real*8 y2(19)


        n=19
        call spline(t_a,psic2,n,1.d30,1.d30,y2)
        call splint(t_a,psic2,y2,n,t0,pp)
100     cal_parti_sic2=pp
	return
	end

        function cal_parti_si(t0)
	implicit real*8 (a-h,o-z)
        real*8 y2(19)

        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi

        n=19
        call spline(t_a,psi,n,1.d30,1.d30,y2)
        call splint(t_a,psi,y2,n,t0,pp)
100     cal_parti_si=pp
	return
	end

        function cal_parti_si2c(t0)
	implicit real*8 (a-h,o-z)
        real*8 y2(19)
        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi

        n=19
        call spline(t_a,psi2c,n,1.d30,1.d30,y2)
        call splint(t_a,psi2c,y2,n,t0,pp)
100     cal_parti_si2c=pp
	return
	end


        function cal_parti_s(t0)
	implicit real*8 (a-h,o-z)

        real*8 y2(19)
        real*8 p_a(19),t_a(19),t_a1(18),psic(19),psic2(19),
     1     psi2c(19),psi(19)
        common /pressure_bock/t_a,p_a,t_a1,psic,psic2,psi2c,psi


        n=19
        call spline(t_a,p_a,n,1.d30,1.d30,y2)
c        pp=rinterp(t_a,p_a,n,t0)
        call splint(t_a,p_a,y2,n,t0,pp)

        goto 100
        open(5,file='pp_spline')
        do i=1,170
        tt=1500+i*10
        call splint(t_a,p_a,y2,n,tt,pp)
        write(5,*)tt,pp
        enddo
        stop

c calculate vapor pressure of sic2 for reaction 2sic=sic2+si
        if(t0.lt.1500)then
        pp=5.131d-8
        goto 100
        else if(t0.gt.3150)then
        pp=17418
        goto 100
        endif
        if(t0.ge.1500.and.t0.le.2000)then
        a1=-34075.8
        b1=15.4274
        else if(t0.ge.2000.and.t0.le.2546)then
        a1=-33526.61
        b1=15.1528
        else if(t0.ge.2546.and.t0.le.2900)then
        a1=-33012.65
        b1=14.9510
        else if(t0.ge.2900.and.t0.le.3150)then
        a1=-32734.33
        b1=14.8550
	endif
	pp=10**(a1/t0+b1)
100     cal_parti_s=pp
	return
	end


      function speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)
	implicit real*8 (a-h,o-z)
c seed temperature t and source temerature t+del_t

      ta0=t+del_t
	pa0=cal_parti_s(ta0)

      pres_sys_m=pres_sys


      total_vapor=cal_parti_sic(ta0)+cal_parti_si(ta0)
     1  +cal_parti_sic2(ta0)+cal_parti_si2c(ta0)


      if(total_vapor.gt.pres_sys)then
      pres_sys_m=total_vapor

c      pa_seed=pa0
c      xapa=1/sqrt(2*3.1415*0.040097*8.31*t)
c      flux=xapa*
c     1    (pres_sys_m-2*cal_parti_s(t))
c      growthsp=flux*8.31*t/pres_sys_m
c      ghk=flux*0.040097/3160*3600*1000
c      speed_fun=ghk
c      return

      endif


	uref=0

1     uref1=uref

c input t
c output growthsp

c ca1 concentration of A at top: including both of SiC2 and Si
c ca0 concentration of A at bottom
c assume c_all and ca0 are fixed

      xapa=1/sqrt(2*3.1415*0.040097*8.31*t)
      pa1=cal_parti_s(t)+uref*pres_sys_m/2/xapa/8.31/t
      cred=(pres_sys_m-2*pa1)/(pres_sys_m-2*pa0)
      d_ab=cal_d0()*(t/273)**1.8*(1d5/pres_sys_m)
      uref=log(cred)*d_ab/(dist)
	uref=uref1+(uref-uref1)*0.1



      if(abs(uref-uref1).ge.1d-6)goto 1

c growthsp: m/s
c molecular weight of sic 0.040097 kg/mol
c density: 3160 kg/m3
c ccen : mol/m3

      growthsp=uref
      pa_seed=pa1
c        write(*,*)pa0,'pascal',growthsp,'m/s',ccen,'mol/m3',t

      ghk=uref*pres_sys_m/8.31/t*0.040097/3160*3600*1000
	speed_fun=ghk
	return
	end

c this subroutine is used to obtain growth rate curves with 
c  pressures for the paper to JCG

	subroutine cal_gr_ster
	implicit real*8 (a-h,o-z)
      character*20 fn
	common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

        write(*,*)'in cal_gr_ster'
      p=200
	t=3000
c temperature at seed
c for Sterling del_t=10
c     ATMI  del_t=100
c ipres=1 : pressure fixed
c       0 : temperature   fixed
c       other: interactive

      ipres=0

      if(ipres.eq.2)then
      write(*,*)'input P(Torr)'
      read(*,*)p
      write(*,*)'input t(K)'
      read(*,*)t
      del_t=100
	pres_sys=133.3*p
      dist=2*0.0254
      ghk=speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)
c        write(*,*)'growth speed',ghk,'mm/hour'

      stop
      endif


      do 2000 id=1,4
      if(id.eq.1)del_t=10
      if(id.eq.2)del_t=20
      if(id.eq.3)del_t=50
      if(id.eq.4)del_t=100
      do 2000 loop=1,9


      if(del_t.eq.10)then
      if(loop.eq.1)then
      fn='gr_2300k10k'
      t=2300
      else if(loop.eq.2)then
      fn='gr_2400k10k'
      t=2400
      else if(loop.eq.3)then
      fn='gr_2500k10k'
      t=2500
      else if(loop.eq.4)then
      fn='gr_2600k10k'
      t=2600
      else if(loop.eq.5)then
      fn='gr_2700k10k'
      t=2700
      else if(loop.eq.6)then
      fn='gr_2800k10k'
      t=2800
      else if(loop.eq.7)then
      fn='gr_2900k10k'
      t=2900
      else if(loop.eq.8)then
      fn='gr_3000k10k'
      t=3000
      else if(loop.eq.9)then
      fn='gr_3100k10k'
      t=3100
      endif

      else if(del_t.eq.50)then
      if(loop.eq.1)then
      fn='gr_2000k50k'
      t=2000
      else if(loop.eq.2)then
      fn='gr_2100k50k'
      t=2100
      else if(loop.eq.3)then
      fn='gr_2200k50k'
      t=2200
      else if(loop.eq.4)then
      fn='gr_2300k50k'
      t=2300
      else if(loop.eq.5)then
      fn='gr_2400k50k'
      t=2400
      else if(loop.eq.6)then
      fn='gr_2500k50k'
      t=2500
      else if(loop.eq.7)then
      fn='gr_2600k50k'
      t=2600
      else if(loop.eq.8)then
      fn='gr_2700k50k'
      t=2700
      else if(loop.eq.9)then
      fn='gr_2800k50k'
      t=2800
      endif

      else if(del_t.eq.100)then

      if(loop.eq.1)then
      fn='gr_2000k100k'
      t=2000
      else if(loop.eq.2)then
      fn='gr_2100k100k'
      t=2100
      else if(loop.eq.3)then
      fn='gr_2200k100k'
      t=2200
      else if(loop.eq.4)then
      fn='gr_2300k100k'
      t=2300
      else if(loop.eq.5)then
      fn='gr_2400k100k'
      t=2400
      else if(loop.eq.6)then
      fn='gr_2500k100k'
      t=2500
      else if(loop.eq.7)then
      fn='gr_2600k100k'
      t=2600
      else if(loop.eq.8)then
      fn='gr_2700k100k'
      t=2700
      else if(loop.eq.9)then
      fn='gr_2800k100k'
      t=2800
      else if(loop.eq.10)then
      fn='gr_2900k100k'
      t=2900
      endif

      endif


        open(1,file=fn)
        write(*,*)'file name',fn
	ibegin=0
      do 1000 loop1=1,50

      if(loop1.le.10)then
      p=loop1*1
      else if(loop1.le.19)then
      p=10+(loop1-10)*10
      else if(loop1.le.1+3*9)then
      p=100+(loop1-19)*100
      else
      p=1000+(loop1-28)*1000
      endif


      pres_sys=p

      dist=2*0.0254

      ghk=speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)

c        write(*,*)'growth speed',ghk,'mm/hour'

      if(ghk.ge.0.and.ghk.le.1000)goto 900
	if(ibegin.eq.0)then
	goto 1000
	else 
	goto 1001
	endif

900	ibegin=1

        write(1,'(4f16.6)')pres_sys,ghk,growthsp,pa_seed

1000    continue
1001    close(1)
2000    continue
	return
	end

c for pressure
      subroutine cal_gr_ster_pres
	implicit real*8 (a-h,o-z)
      character*20 fn
		common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

      write(*,*)'in cal_gr_ster_pres'
      p=200
	t=3000
c temperature at seed
c for Sterling del_t=10
c     ATMI  del_t=100
c ipres=1 : pressure fixed
c       0 : temperature   fixed
c       other: interactive



      do 2000 id=1,5
      if(id.eq.1)del_t=10
      if(id.eq.2)del_t=20
      if(id.eq.3)del_t=50
      if(id.eq.4)del_t=100
      if(id.eq.5)del_t=250
      do 2000 loop=1,7


      if(del_t.eq.10)then

      if(loop.eq.1)then
      fn='gr_5torr10k'
      p=5
      else if(loop.eq.2)then
      fn='gr_10torr10k'
      p=10
      else if(loop.eq.3)then
      fn='gr_20torr10k'
      p=20
      else if(loop.eq.4)then
      fn='gr_50torr10k'
      p=50
      else if(loop.eq.5)then
      fn='gr_100torr10k'
      p=100
      else if(loop.eq.6)then
      fn='gr_200torr10k'
      p=200
      else if(loop.eq.7)then
      fn='gr_300torr10k'
      p=300
      endif

      else if(del_t.eq.20)then
      if(loop.eq.1)then
      fn='gr_5torr20k'
      p=5
      else if(loop.eq.2)then
      fn='gr_10torr20k'
      p=10
      else if(loop.eq.3)then
      fn='gr_20torr20k'
      p=20
      else if(loop.eq.4)then
      fn='gr_50torr20k'
      p=50
      else if(loop.eq.5)then
      fn='gr_100torr20k'
      p=100
      else if(loop.eq.6)then
      fn='gr_200torr20k'
      p=200
      else if(loop.eq.7)then
      fn='gr_300torr20k'
      p=300
      endif

      else if(del_t.eq.50)then
      if(loop.eq.1)then
      fn='gr_5torr50k'
      p=5
      else if(loop.eq.2)then
      fn='gr_10torr50k'
      p=10
      else if(loop.eq.3)then
      fn='gr_20torr50k'
      p=20
      else if(loop.eq.4)then
      fn='gr_50torr50k'
      p=50
      else if(loop.eq.5)then
      fn='gr_100torr50k'
      p=100
      else if(loop.eq.6)then
      fn='gr_200torr50k'
      p=200
      else if(loop.eq.7)then
      fn='gr_300torr50k'
      p=300
      endif


      else if(del_t.eq.100)then

      if(loop.eq.1)then
      fn='gr_5torr100k'
      p=5
      else if(loop.eq.2)then
      fn='gr_10torr100k'
      p=10
      else if(loop.eq.3)then
      fn='gr_20torr100k'
      p=20
      else if(loop.eq.4)then
      fn='gr_50torr100k'
      p=50
      else if(loop.eq.5)then
      fn='gr_100torr100k'
      p=100
      else if(loop.eq.6)then
      fn='gr_200torr100k'
      p=200
      else if(loop.eq.7)then
      fn='gr_300torr100k'
      p=300
      endif
      else if(del_t.eq.250)then

      if(loop.eq.1)then
      fn='gr_5torr250k'
      p=5
      else if(loop.eq.2)then
      fn='gr_10torr250k'
      p=10
      else if(loop.eq.3)then
      fn='gr_20torr250k'
      p=20
      else if(loop.eq.4)then
      fn='gr_50torr250k'
      p=50
      else if(loop.eq.5)then
      fn='gr_100torr250k'
      p=100
      else if(loop.eq.6)then
      fn='gr_200torr250k'
      p=200
      else if(loop.eq.7)then
      fn='gr_300torr250k'
      p=300
      endif


      endif



        open(1,file=fn)
	ibegin=0
      do 1000 loop1=0,200
      t=2000+loop1*10
      if(t.gt.3200)goto 1001

	pres_sys=133.3*p

      dist=2*0.0254

      ghk=speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)

c        write(*,*)'growth speed',ghk,'mm/hour'
      if(ghk.ge.0.and.ghk.le.1000)goto 900
	if(ibegin.eq.0)then
	goto 1000
	else 
	goto 1001
	endif

900	ibegin=1

        write(1,'(4f16.6)')10000/t,ghk,growthsp,pa_seed

1000    continue
1001    close(1)
2000    continue
	return
	end


c for temperature gradient
      subroutine cal_gr_ster_gradient
	implicit real*8 (a-h,o-z)
      character*20 fn
		common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

      write(*,*)'in cal_gr_ster_gradient'
      p=200
	t=3000
c temperature at seed
c for Sterling del_t=10
c     ATMI  del_t=100
c ipres=1 : pressure fixed
c       0 : temperature   fixed
c       other: interactive


      do 2000 id=1,6
c assuming distant of 2 inch
      if(id.eq.1)del_t=5
      if(id.eq.2)del_t=10
      if(id.eq.3)del_t=25
      if(id.eq.4)del_t=50
      if(id.eq.5)del_t=100
      if(id.eq.6)del_t=250

      do 2000 loop=1,2

      if(loop.eq.1)then
      t=2400
      if(id.eq.1)then
      fn='gr_2400k5k'
      else if(id.eq.2)then
      fn='gr_2400k10k'
      else if(id.eq.3)then
      fn='gr_2400k25k'
      else if(id.eq.4)then
      fn='gr_2400k50k'
      else if(id.eq.5)then
      fn='gr_2400k100k'
      else if(id.eq.6)then
      fn='gr_2400k250k'
      endif
      else
      t=2500
      if(id.eq.1)then
      fn='gr_2500k5k'
      else if(id.eq.2)then
      fn='gr_2500k10k'
      else if(id.eq.3)then
      fn='gr_2500k25k'
      else if(id.eq.4)then
      fn='gr_2500k50k'
      else if(id.eq.5)then
      fn='gr_2500k100k'
      else if(id.eq.6)then
      fn='gr_2500k250k'
      endif
      endif

        open(1,file=fn)

	ibegin=0

      do 1000 loop1=1,50
      if(loop1.le.10)then
      p=loop1*1
      else if(loop1.le.19)then
      p=10+(loop1-10)*10
      else if(loop1.le.1+3*9)then
      p=100+(loop1-19)*100
      else
      p=1000+(loop1-28)*1000
      endif


      pres_sys=p

      dist=2*0.0254

      ghk=speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)

c        write(*,*)'growth speed',ghk,'mm/hour'
      if(ghk.ge.0.and.ghk.le.1000)goto 900
	if(ibegin.eq.0)then
	goto 1000
	else 
	goto 1001
	endif

900	ibegin=1

        write(1,'(4f16.6)')pres_sys,ghk,growthsp,pa_seed

1000    continue
1001    close(1)
2000    continue
	return
	end

c for temperature gradient effect on Gr-1/T
      subroutine cal_gr_ster_gradient_t
	implicit real*8 (a-h,o-z)
      character*20 fn
		common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

      write(*,*)'in cal_gr_ster_gradient_t'
      p=100
	t=3000
c temperature at seed
c for Sterling del_t=10
c     ATMI  del_t=100
c ipres=1 : pressure fixed
c       0 : temperature   fixed
c       other: interactive


      do 2000 id=1,6
c assuming distant of 2 inch
      if(id.eq.1)del_t=5
      if(id.eq.2)del_t=10
      if(id.eq.3)del_t=25
      if(id.eq.4)del_t=50
      if(id.eq.5)del_t=100
      if(id.eq.6)del_t=250

      do 2000 loop=1,6

c system pressure in Torr

      if(loop.eq.1)then
      p=5
      if(id.eq.1)then
      fn='gr_5torr5k'
      else if(id.eq.2)then
      fn='gr_5torr10k'
      else if(id.eq.3)then
      fn='gr_5torr25k'
      else if(id.eq.4)then
      fn='gr_5torr50k'
      else if(id.eq.5)then
      fn='gr_5torr100k'
      else if(id.eq.6)then
      fn='gr_5torr250k'
      endif
      else if(loop.eq.2)then
      p=10
      if(id.eq.1)then
      fn='gr_10torr5k'
      else if(id.eq.2)then
      fn='gr_10torr10k'
      else if(id.eq.3)then
      fn='gr_10torr25k'
      else if(id.eq.4)then
      fn='gr_10torr50k'
      else if(id.eq.5)then
      fn='gr_10torr100k'
      else if(id.eq.6)then
      fn='gr_10torr250k'
      endif

      else if(loop.eq.3)then
      p=20
      if(id.eq.1)then
      fn='gr_20torr5k'
      else if(id.eq.2)then
      fn='gr_20torr10k'
      else if(id.eq.3)then
      fn='gr_20torr25k'
      else if(id.eq.4)then
      fn='gr_20torr50k'
      else if(id.eq.5)then
      fn='gr_20torr100k'
      else if(id.eq.6)then
      fn='gr_20torr250k'
      endif

      else if(loop.eq.4)then
      p=50
      if(id.eq.1)then
      fn='gr_50torr5k'
      else if(id.eq.2)then
      fn='gr_50torr10k'
      else if(id.eq.3)then
      fn='gr_50torr25k'
      else if(id.eq.4)then
      fn='gr_50torr50k'
      else if(id.eq.5)then
      fn='gr_50torr100k'
      else if(id.eq.6)then
      fn='gr_50torr250k'
      endif

      else  if(loop.eq.5)then
      p=100
      if(id.eq.1)then
      fn='gr_100torr5k'
      else if(id.eq.2)then
      fn='gr_100torr10k'
      else if(id.eq.3)then
      fn='gr_100torr25k'
      else if(id.eq.4)then
      fn='gr_100torr50k'
      else if(id.eq.5)then
      fn='gr_100torr100k'
      else if(id.eq.6)then
      fn='gr_100torr250k'
      endif

      else if(loop.eq.6)then
      p=200
      if(id.eq.1)then
      fn='gr_200torr5k'
      else if(id.eq.2)then
      fn='gr_200torr10k'
      else if(id.eq.3)then
      fn='gr_200torr25k'
      else if(id.eq.4)then
      fn='gr_200torr50k'
      else if(id.eq.5)then
      fn='gr_200torr100k'
      else if(id.eq.6)then
      fn='gr_200torr250k'
      endif

      endif

        open(1,file=fn)

	ibegin=0

      do 1000 loop1=0,200
      t=2000+loop1*10


	pres_sys=133.3*p

      dist=2*0.0254

      ghk=speed_fun(t,del_t,pres_sys,dist,growthsp,pa_seed)

c        write(*,*)'growth speed',ghk,'mm/hour'
      if(ghk.ge.0.and.ghk.le.100000)goto 900
	if(ibegin.eq.0)then
	goto 1000
	else 
	goto 1001
	endif

900	ibegin=1

        write(1,'(4f16.6)')10000/t,ghk,growthsp,pa_seed

1000    continue
1001    close(1)
2000    continue
	return
	end



      subroutine sticking
      implicit real*8 (a-h,o-z)
		common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

      open(1,file='hk.dat')
      open(2,file='sticking.dat')

      stick=0.016
      del_t=50
      pres_sys=7000
      do t_inv=3.5,5,0.02
      t=10000/t_inv
      gr=stick*cal_parti_s(t)*0.606/sqrt(t)
     1  *0.040097/3160*3600*1000
c change to mm/hr
      ghk=speed_fun(t,del_t,pres_sys,0.005d0,growthsp,pa_seed)
      write(1,*)t_inv,gr/stick
      write(2,*)t_inv,gr
      enddo
      close(1)
      close(2)
	return
      end


	subroutine cal_gr_advan
	implicit real*8 (a-h,o-z)
      character*20 fn
	common /cal_d0_v/cal_d0_v
	cal_d0_v=5d-6

	t=2000
      write(*,*)'input p(torr) t(K)'
      read(*,*)p,t
      pres_sys=133.3*p
	ta0=t+150
	pa0=cal_parti_s(ta0)
	ca0=pa0/8.31/ta0
	c_all=pres_sys/8.31/ta0
	write(*,*)'pa0=',pa0,'c_all',c_all

      pres_sys_m=pres_sys

      total_vapor=cal_parti_sic(ta0)+cal_parti_si(ta0)
     1  +cal_parti_sic2(ta0)+cal_parti_si2c(ta0)

      if(total_vapor.gt.pres_sys)then
      pres_sys_m=total_vapor
      endif


	uref=0

10	uref1=uref
c input t
c output growthsp

c ca1 concentration of A at top: including both of SiC2 and Si
c ca0 concentration of A at bottom
c assume c_all and ca0 are fixed

      dist=2*0.0254

      xapa=1/sqrt(2*3.1415*0.040097*8.31*t)
      pa1=cal_parti_s(t)+uref*pres_sys_m/2/xapa/8.31/t
      cred=(pres_sys_m-2*pa1)/(pres_sys_m-2*pa0)
      d_ab=cal_d0()*(t/273)**1.8*(1d5/pres_sys_m)
      uref=log(cred)*d_ab/(dist)
	uref=uref1+(uref-uref1)*0.1


	write(*,*)'uref',uref,'ca0',ca0,'ca1',ca1,'pa0',pa0
	if(abs(uref-uref1).ge.1d-6)goto 10

c growthsp: m/s
c molecular weight 0.040097 kg/mol
c density: 3160 kg/m3
c ccen : mol/m3

c	growthsp=0.606/sqrt(t)*(8.31*t
c     1   -cal_parti(t)/ca1)
	growthsp=uref
	write(*,*)pa0,'pascal',growthsp,'m/s',ccen,'mol/m3',t

c	ghk=2*growthsp*ca1*0.040097/3160*3600*1000
      ghk=uref*pres_sys_m/8.31/t*0.040097/3160*3600*1000

c      write(*,*)'growth speed',ghk,'mm/hour'
c      write(*,*)'total pressure',total_vapor,'system pres',pres_sys

	return
	end
