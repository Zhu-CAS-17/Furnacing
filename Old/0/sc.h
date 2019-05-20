c   idesign: 1- 2 inch,(2/1998) 2 -- 3 inch (5/1998) 3 -- 2 inch  (1/1999)
c   4 -- 3 inch   (12/1999)
c   5 - 5''
c   6 - 3'' 02 Furnace 3'' crucible and heat pack as of 16 Feb 2000


cdif: 0 1st order for time integration
c     1 2nd order for time integration
      implicit real*8 (a-h,o-z)
c (f1,1,8),fc
      parameter(ni=300,nj=220,nk=14,nmaxij=300,ns=10,iprint=1,
     1    itest=1,cdif=0,nnn=300)
      parameter(pai=3.1415926535857932,pi=3.14159265)

       parameter (nfch=5)
      common /changegrid/xo(ni,nj),yo(ni,nj),
     1  uch(ni,nj,nfch)


      logical lsolve,lprint,lblk,lstop,lortho,lpgrid,lconv,lcon
      common /idesign/idesign
      common /basic/f(ni,nj,ns),rho(ni,nj),rho_o(ni,nj),
     1   gam(ni,nj),xc(ni,nj),yc(ni,nj),
     1   x(ni,nj),y(ni,nj),r(ni,nj),fo(ni,nj,ns),con(ni,nj),ap(ni,nj),
     2   aip(ni,nj),aim(ni,nj),ajp(ni,nj),ajm(ni,nj),ak1(ni,nj),
     3   ak2(ni,nj),ae1(ni,nj),ae2(ni,nj),vol(ni,nj),hksi(ni,nj),
     4   heta(ni,nj),ruksi(ni,nj),rueta(ni,nj),fjksi(ni,nj),
     5   fjeta(ni,nj),con0(ni,nj),ap0(ni,nj),parti_p(ni,nj),
     1   ex1(ni,nj),scargon(ni,nj),rho_s(ni,nj),cp(ni,nj),
     1   smass(ni,nj),growthsp(nj)
c 6: concentration
c 7: w
c 8,9: magnetic potential
	common /radi_co/rw(nnn),zw(nnn),tw(nnn),qw(nnn),pw(nnn),
     1    rs(nnn),zs(nnn),qsurf(nnn),tsurf(nnn),ps(nnn),
     1    rb(nnn),zb(nnn),tb(nnn),qb(nnn),pb(nnn),
     1    arb(nnn),arw(nnn),ars(nnn)
	common /latent/gr_lat,zb_lat,zt_lat

c i_rt: grid point just under the upper lid of the chamber
c i_rb: lower point of the chamber
c j_rr: right point of the chamber
	common /radi_con/i_rt,i_rb,j_rr
c z_pres-z_pres2 where the maximum pressure is taken
c z_incr used to add crucible height

      common /radi_con1/rc,zb0,zt0,z_pres,z_pres2,z_incr,
     1  tmax,tmin,
     1    rr_c_t,zb_c_t,zt_c_t,rr_c_l,zb_c_l,zt_c_l,
     1    tos_a,tos_d,tos_h1,tos_h2,tos_h
c coil_move: move coil position vertically (m)
c felt_incr: increase or decrease thickness of the felt above the susceptorm (m)
      common /coils/coil_xx0,coil_yy0,coil_dxx,width_coil,
     1   xx_r,yy_r,
     1  coil_move,felt_incr,t_ref,t_top,t_bot,p_int

        common /geom/xx_0(50),xx_1(50),
     1   yy_0(50),yy_1(50),id_el0(50),air_len0(50),id_total

	common/id_material/id_graphite,id_hard,id_soft,id_argon,id_sic,
     1   id_powder,id_foam,id_max,
c imaxp,jmaxp, where powder source has highest T, i_powder 0: powder, 1:carbon 
     1   imaxp,jmaxp,icontrol
      common /i_powder/i_powder(ni,nj),id_bound(ni,nj)

      common /iid/iid(6000),jjd(6000)

      common/indx1/nf,nfmax,np,npc,l1,l2,l3,m1,m2,m3,mode,iter,n,nt
      common/indx2/res(ns),relax(nk),relax_b,rhocon,
     1     iout_dat,icom_mag,
     1    ipref,jpref,ndt,nsolve(ns),
     1       lsolve(ns),lprint(nk),lblk(ns),ntimes(ns),lortho,lpgrid


      common/bndj/fjbi1(nj,ns),fjbj1(ni,ns),fjbl1(nj,ns),fjbm1(ni,ns)
      common/cntl/lstop,lconv,lcon(ns)
      common/sorc/smax,ssum,epsil,cfo,dtm,tinit,tpass,tend,
     1   blc
      common/grd/rdtm,xlen,ylen
      common/const/re,rre,pr,rpr,ra,stes,stel,fksl,fcsl,frsl,hf,tf
      common/cons2/tin,fr,rbo,tanca,tanca2,fma,hamag,sigma,grash,
     1   darcy,pres_sys,crys_len,air_len
c pres_sys: system pressure (Pa)
c crys_len: crystal length
      common/lst/b1jbl(ni),b2jbl(ni),lst,m0,lst2,l0,izone,lst3,
     1  ifst1,ilst
      common/in/xin(ni),yin(ni),inum
      common/ini/init
      dimension u(ni,nj),v(ni,nj),w(ni,nj),pc(ni,nj),p(ni,nj),t(ni,nj),
     1   ccen(ni,nj)
      equivalence(f(1,1,1),u(1,1)),(f(1,1,2),v(1,1)),(f(1,1,3),p(1,1)),
     1   (f(1,1,4),pc(1,1)),(f(1,1,5),t(1,1))
      equivalence (f(1,1,6),ccen(1,1)),(f(1,1,7),w(1,1))
      common /dhdy/ dhdy(nj),dhdym(nj),dhdyp(nj),d2hdy2(nj),rr(nj),
     1  bx(ni),by(ni),
     1   tau11(ni,nj),tau12(ni,nj),tau22(ni,nj),
     1    pt(nmaxij),qt(nmaxij),
     1   fn(ni,nj,ns),
     1   hutop(ni,nj),hvtop(ni,nj),tmp(ni,nj),
     1   ptc,qtc,apc

      complex*16 ptc(nmaxij),qtc(nmaxij),fc(ni,nj),fco(ni,nj),apc(ni,nj)

      common/xtrms/dudx(ni,nj),dudy(ni,nj),dvdx(ni,nj),dvdy(ni,nj),
     1   dwdx(ni,nj),dwdy(ni,nj)
      common/pmi/b3jbl(ni),fluxm1,fluxm2
      common/save/cofu(ni,nj,5),cofv(ni,nj,5),conu(ni,nj),conv(ni,nj),
     1      biu(nj),biv(nj),bju(ni),bjv(ni),blu(nj),blv(nj),
     2      bmu(ni),bmv(ni)

      character title*7,titlex*7,titley*7
      common/titl/title(nk),titlex,titley

	dimension a_mag(2,ni,nj)
	equivalence (f(1,1,8),fc(1,1)),(fo(1,1,8),fco(1,1))
	equivalence(f(1,1,8),a_mag(1,1,1))
	common/coef/flow,diff,acof

      common /elemag/xmiu,epsilon,omega,cden0,cden,xxl,yyl,iphase
	common /porous/qth(ni,nj),b_por(ni,nj),void(ni,nj)

      common /pksi/pksi(ni,nj),peta(ni,nj),hksip(ni,nj),hetap(ni,nj)


      common /grid_info/nxx,nyy,
     1   xxd(50),yyd(50),iddx(50),iddy(50)

