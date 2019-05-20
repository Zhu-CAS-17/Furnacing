      Module meshElectroMagnetic
      Implicit None
c mode=0 rectangular =1 cylinder coordinate
      Integer(kind=4),Parameter::mode=1
      Real(kind=8),Parameter::pai=3.1415926,pi=3.1415926
c input information of design mesh
      Integer(kind=4)::nxx,nyy,num_xtem,num_ytem
      Real(kind=8)::long_x,long_y
      Real(kind=8),Allocatable::i_dx(:),i_dy(:),x_dx(:),y_dy(:)
c calculation zone and element info
      Real(kind=8),Allocatable::xc(:,:),yc(:,:),x(:,:),y(:,:)
      Real(kind=8),Allocatable::hksi(:,:),heta(:,:),vol(:,:)
      Real(kind=8),Allocatable::ak1(:,:),ak2(:,:),ae1(:,:),ae2(:,:)
      Real(kind=8),Allocatable::r(:,:)
      Contains
              Subroutine readInMeshData()
              Implicit None
              Integer(kind=4)::i
              Write(*,*)'readin mesh data from grid.dat'
              Open(10,file='grid.dat')
              Read(10,*)
              Read(10,*)nxx
              Allocate(i_dx(nxx),x_dx(nxx))
              Do i=1,nxx
              Read(10,*)i_dx(i),x_dx(i)
              Enddo
              Read(10,*)
              Read(10,*)nyy
              Allocate(i_dy(nyy),y_dy(nyy))
              Do i=1,nyy
              Read(10,*)i_dy(i),y_dy(i)
              Enddo
              Close(10)
              Return
              End Subroutine

              Subroutine creatTemMesh()
              Implicit None
              Integer(kind=4)::i,j,k,id,idx,idy
              Write(*,*)'creat initial mesh'
              long_x=x_dx(nxx)-x_dx(1)
              num_xtem=2
              Do i=2,nxx
              If(x_dx(i).le.x_dx(i-1))Stop 'error in grid:x(i)>x(i+1)'
              num_xtem=num_xtem+i_dx(i)
              Enddo
              long_y=y_dy(nyy)-y_dy(1)
              num_ytem=2
              Do i=2,nyy
              If(y_dy(i).le.y_dy(i-1))Stop 'error in grid:y(i)>y(i+1)'
              num_ytem=num_ytem+i_dy(i)
              Enddo
              
              Allocate(xc(num_xtem,num_ytem),yc(num_xtem,num_ytem))

              i=2
              Do j=2,num_ytem
              xc(2,j)=x_dx(1)
              Enddo

              Do id=2,nxx
              idx=i_dx(id)
              If(idx.lt.1)idx=1
              Write(*,*)idx
              If(id.eq.nxx)idx=num_xtem-i
              Do k=1,idx
              Do j=2,num_ytem
              xc(i+1,j)=x_dx(id-1)+(x_dx(id)-x_dx(id-1))*k/idx
              Enddo
              i=i+1
              Enddo
              Enddo
              
              j=2
              Do i=2,num_xtem
              yc(i,2)=y_dy(1)
              Enddo

              Do id=2,nyy
              idy=i_dy(id)
              Write(*,*)idy
              If(idy.lt.1)idy=1
              If(id.eq.nyy)idy=num_ytem-j
              Do k=1,idy
              Do i=2,num_xtem
              yc(i,j+1)=y_dy(id-1)+(y_dy(id)-y_dy(id-1))*k/idy
              Enddo
              j=j+1
              Enddo
              Enddo
              Return
              End Subroutine

              Subroutine electroCalculationInfo()
              Implicit None
              Integer(kind=4)::i,j
              Real(kind=8)::xp,yp,xm,ym,xu,yu,xd,yd
              Real(kind=8)::dxdksi,dydksi,dxdeta,dydeta
              Real(kind=8)::t1,t2,t3,xjacb
              Write(*,*)'calculate nodes,areas,volumes'
              Allocate(x(num_xtem,num_ytem),y(num_xtem,num_ytem))
              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
              x(i,j)=0.25*(xc(i,j)+xc(i+1,j)+xc(i,j+1)+xc(i+1,j+1))
              y(i,j)=0.25*(yc(i,j)+yc(i+1,j)+yc(i,j+1)+yc(i+1,j+1))
              Enddo
              Enddo
              Do i=2,num_xtem-1
              x(i,1)=0.5*(xc(i,2)+xc(i+1,2))
              y(i,1)=0.5*(yc(i,2)+yc(i+1,2))
              x(i,num_ytem)=0.5*(xc(i,num_ytem)+xc(i+1,num_ytem))
              y(i,num_ytem)=0.5*(yc(i,num_ytem)+yc(i+1,num_ytem))
              Enddo
              Do j=2,num_ytem-1
              x(1,j)=0.5*(xc(2,j)+xc(2,j+1))
              y(1,j)=0.5*(yc(2,j)+yc(2,j+1))
              x(num_xtem,j)=0.5*(xc(num_xtem,j)+xc(num_xtem,j+1))
              y(num_xtem,j)=0.5*(yc(num_xtem,j)+yc(num_xtem,j+1))
              Enddo
              x(1,1)=xc(2,2)
              y(1,1)=yc(2,2)
              x(num_xtem,num_ytem)=xc(num_xtem,num_ytem)
              y(num_xtem,num_ytem)=xc(num_xtem,num_ytem)
              x(1,num_ytem-1)=xc(2,num_ytem-1)
              y(1,num_ytem-1)=yc(2,num_ytem-1)
              x(num_xtem,1)=xc(num_xtem,2)
              y(num_xtem,1)=xc(num_xtem,2)
c calculate local scale factors and volumes
c hksi heta vol are all basing on the domain control volume
              Allocate(hksi(num_xtem,num_ytem),heta(num_xtem,num_ytem))
              Allocate(vol(num_xtem,num_ytem))
              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
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
              Enddo
              Enddo

              Do i=1,num_xtem
              heta(i,1)=0.0
              heta(i,num_ytem)=0.0
              vol(i,1)=0.0
              vol(i,num_ytem)=0.0
              If(i.eq.1.or.i.eq.num_xtem)then
              hksi(i,1)=0.0
              hksi(i,num_ytem)=0.0
              else
      hksi(i,1)=sqrt((xc(i+1,2)-xc(i,2))**2+(yc(i+1,2)-yc(i,2))**2)
      hksi(i,num_ytem)=sqrt((xc(i+1,num_ytem)-xc(i,num_ytem))**2+
     1 (yc(i+1,num_ytem)-yc(i,num_ytem))**2)
              Endif
              Enddo

              Do j=1,num_ytem
              hksi(1,j)=0.0
              hksi(num_xtem,j)=0.0
              vol(1,j)=0.0
              vol(num_xtem,j)=0.0
              If(j.eq.1.or.j.eq.num_ytem)then
              heta(1,j)=0.0
              heta(num_xtem,j)=0.0
              else
      heta(1,j)=sqrt((xc(2,j)-xc(2,j+1))**2+(yc(2,j)-yc(2,j+1))**2)
      heta(num_xtem,j)=sqrt((xc(num_xtem,j)-xc(num_xtem,j+1))**2+
     1 (yc(num_xtem,j)-yc(num_xtem,j+1))**2)
              Endif
              Enddo
c calcaulate areas on the control-volume faces(staggered grids)
              Allocate(ak1(num_xtem,num_ytem),ak2(num_xtem,num_ytem))
              Allocate(ae1(num_xtem,num_ytem),ae2(num_xtem,num_ytem))
              Allocate(r(num_xtem,num_ytem))
              Do j=2,num_ytem-1
              Do i=2,num_xtem
              dxdeta=xc(i,j+1)-xc(i,j)
              dydeta=yc(i,j+1)-yc(i,j)
              dxdksi=x(i,j)-x(i-1,j)
              dydksi=y(i,j)-y(i-1,j)
              If(i.eq.2.or.i.eq.num_xtem)dxdksi=dxdksi*2.0
              If(i.eq.2.or.i.eq.num_xtem)dydksi=dydksi*2.0
              t1=dxdeta**2+dydeta**2
              t2=dxdksi**2+dydksi**2
              xjacb=dxdksi*dydeta-dxdeta*dydksi
              t3=dxdksi*dxdeta+dydeta*dydksi
              ak1(i,j)=sqrt(t2)*t1/xjacb
              ak2(i,j)=sqrt(t1)*t3/xjacb
              Enddo
              Enddo

              Do i=2,num_xtem-1
              Do j=2,num_ytem
              dxdeta=x(i,j)-x(i,j-1)
              dydeta=y(i,j)-y(i,j-1)
              dxdksi=xc(i+1,j)-xc(i,j)
              dydksi=yc(i+1,j)-yc(i,j)
              If(j.eq.2.or.j.eq.num_ytem)dxdeta=dxdeta*2.0
              If(j.eq.2.or.j.eq.num_ytem)dydeta=dydeta*2.0
              t1=dxdeta**2+dydeta**2
              t2=dxdksi**2+dydksi**2
              xjacb=dxdksi*dydeta-dxdeta*dydksi
              t3=dxdksi*dxdeta+dydeta*dydksi
              ae1(i,j)=sqrt(t1)*t2/xjacb
              ae2(i,j)=sqrt(t2)*t3/xjacb
              Enddo
              Enddo
              
              Do i=1,num_xtem
              Do j=1,num_ytem
              If(mode.eq.0) r(i,j)=1.0
              If(mode.eq.1) r(i,j)=y(i,j)+1.0d-8
              If(i.ge.2)then
              ak1(i,j)=ak1(i,j)*r(i,j)*2*pai
              ak2(i,j)=ak2(i,j)*r(i,j)*2*pai
              Endif
              If(j.ge.2)then
              ae1(i,j)=ae1(i,j)*yc(i,j)*2*pai
              ae2(i,j)=ae2(i,j)*yc(i,j)*2*pai
              Endif
              vol(i,j)=vol(i,j)*r(i,j)*2*pai
              Enddo
              Enddo
              Return
              End Subroutine

      End Module meshElectroMagnetic

      Module fieldPhysical
      Use meshElectroMagnetic
      Implicit None
c magnetic potential A f: far Real;fai Image
      Real(kind=8),Allocatable::t(:,:),far(:,:),fai(:,:)
      Contains
              Subroutine fieldZeroAndInitial
              Implicit None
              Integer(kind=4)::i,j
              Allocate(far(num_xtem,num_ytem),fai(num_xtem,num_ytem))
              Do i=1,num_xtem
              far(i,1)=0
              fai(i,1)=0
              far(i,num_ytem)=0
              fai(i,num_ytem)=0
              Enddo
              Do j=1,num_ytem
              far(1,j)=0
              fai(1,j)=0
              far(num_xtem,j)=0
              fai(num_xtem,j)=0
              Enddo
              Return
              End Subroutine
      End Module fieldPhysical

      Module materialElectroMagnetic
      Use meshElectroMagnetic
      Use fieldPhysical
      Implicit None
      Real(kind=8),Parameter::xmiu=4*3.1415926*1.0e-7
      Real(kind=8),Parameter::epsil=8.8542e-12
      Real(kind=8)::freq,omega
      Integer(kind=4)::id_total
      Real(kind=8)::coil_xx0,coil_yy0,coil_dxx,width_coil,xx_r,yy_r
      Real(kind=8),Allocatable::xx_0(:),xx_1(:),yy_0(:),yy_1(:),
     1 air_len0(:)
      Integer(kind=4),Allocatable::id_el0(:)
      Real(kind=8)::air_len
      Real(kind=8),Allocatable::conduc_d(:)
      Real(kind=8),Allocatable::con_sict(:),con_sic(:)
      Real(kind=8),Allocatable::con_grat(:),con_gra(:)
      Contains
              Subroutine readInCoil()
              Implicit None
              Open(10,file='coil.dat')
              Read(10,*)freq
              Read(10,*)coil_xx0,coil_yy0,coil_dxx,width_coil,xx_r,yy_r
              Close(10)
              omega=freq*2*3.1415926
              Return
              End Subroutine

              Subroutine readInFurnace
              Implicit None
              Integer(kind=4)::i
              Open(11,file='geom.dat')
              Do i=1,10
              Read(11,*)
              Enddo
              Read(11,*)id_total
              Allocate(xx_0(id_total),xx_1(id_total),yy_0(id_total))
              Allocate(yy_1(id_total),id_el0(id_total))
              ALlocate(air_len0(id_total))
              Do i=1,id_total
              Read(11,*)
              Read(11,*)xx_0(i),xx_1(i),yy_0(i),yy_1(i),id_el0(i)
              If(id_el0(i).eq.10)Read(11,*)air_len0(i)
              Enddo
              Close(11)
              Return
              End Subroutine

              Subroutine readInMaterial
              Implicit None
              Integer(kind=4)::i,num
              Open(12,file='mater.dat')
              Read(12,*)
              Read(12,*)num
              Allocate(con_sic(num),con_sict(num))
              Do i=1,num
              Read(12,*)con_sic(i),con_sict(i)
              Enddo
              Read(12,*)
              Read(12,*)num
              Allocate(con_gra(num),con_grat(num))
              Do i=1,num
              Read(12,*)con_gra(i),con_grat(i)
              Enddo
              Read(12,*)
              Read(12,*)num
              Allocate(conduc_d(num))
              Do i=1,num
              Read(12,*)conduc_d(i)
              Enddo
              Close(12)
              Return
              End Subroutine

              Integer(kind=4) Function insd(x,y,uu1,uu2,vv1,vv2)
              Real(kind=8)::x,y,uu1,uu2,vv1,vv2
              Real(kind=8)::x1,x2,y1,y2
              insd=0
              x1=uu1
              x2=uu2
              y1=vv1
              y2=vv2
              If(x1.gt.x2)then
              x2=uu1
              x1=uu2
              y2=vv1
              y1=vv2
              Endif
              If(x.ge.x1.and.x.le.x2.and.y.ge.y1.and.y.le.y2)insd=1
              Return
              End Function

              Integer(kind=4) Function id_ele(i1,j1)
              Implicit None
              Integer(kind=4)::i,i1,j1
              Real(kind=8)::xx,yy
              id_ele=0
              xx=x(i1,j1)
              yy=y(i1,j1)
              Do i=1,id_total
              If(insd(xx,yy,xx_0(i),xx_1(i),yy_0(i),yy_1(i)).eq.1)then
              id_ele=id_el0(i)
              If(id_ele.eq.10)air_len=air_len0(i)
              Endif
              Enddo
              Return
              End Function

              Real(kind=8) Function rmu_c(i1,j1)
              Implicit None
              Integer(kind=4)::i1,j1,id
              rmu_c=1
              id=id_ele(i1,j1)
              If(id.eq.2)then
              rmu_c=1
              Endif
              If(id.eq.5)then
              rmu_c=1
              Endif
              Return
              End Function

              Real(kind=8) Function conduc(i,j)
              Implicit None
              Integer(kind=4)::i,j,id,iee
              Real(kind=8)::tcel,temp,ee
              conduc=0
              id=id_ele(i,j)
              If(id.eq.2.or.(id.ge.5.and.id.le.12))then
              conduc=conduc_d(id)
              Endif
              tcel=t(i,j)-273
              ee=(t(i,j)-273)/250.0+1
              iee=int(ee)
              If(iee.gt.9)iee=9
              If(iee.lt.1)iee=1
              If(id.eq.1.or.id.eq.8)then
              conduc=rinterp(con_grat,con_gra,10,tcel)
              conduc=conduc*1e4
              Endif
              If(id.eq.3)then
              temp=t(i,j)-273
              If(temp.le.1000)then
              conduc=con_sic(1)
              Else If(temp.ge.2500)then
              conduc=con_sic(7)
              Else
              ee=(temp-1000)/250.0+1
              iee=int(ee)
              conduc=con_sic(iee)*(1+iee-ee)+con_sic(iee+1)*(ee-iee)
              Endif
              Endif
              
              If(id.eq.13)then
              conduc=rinterp(con_sict,con_sic,7,tcel)
              Endif
              Return
              End Function
              
              Real(kind=8) Function rinterp(x,y,n,x0)
              Implicit None
              Real(kind=8),Allocatable::x(:),y(:)
              Integer(kind=4)::n,i,ii
              Real(kind=8)::x0
              Allocate(x(n),y(n))
              If(x(1).gt.x(n))Stop 'error in rinterp'
              If(x0.lt.x(1))ii=1
              If(x0.gt.x(n))ii=n-1
              Do i=1,n-1
              If(x0.ge.x(i).and.x0.le.x(i+1))then
              ii=i
              Endif
              Enddo
              rinterp=y(ii)+(y(ii+1)-y(ii))*(x0-x(ii))/(x(ii+1)-x(ii))
              Return
              End Function

      End Module materialElectroMagnetic
      
      Module solverCalculate
      Use meshElectroMagnetic
      Use materialElectroMagnetic
      Implicit none
c magnetic potential A f: far Real;fai Image
c      Real(kind=8),Allocatable::far(:,:),fai(:,:)
      Real(kind=8),Allocatable::ap(:,:),con(:,:),aip(:,:),aim(:,:)
      Real(kind=8),Allocatable::ajp(:,:),ajm(:,:),apc(:,:)
      Contains
              Subroutine elctroMagneticCofficience
              Implicit None
              Integer(kind=4)::i,j
              Real(kind=8)rmu,conduc0
              Allocate(ap(num_xtem,num_ytem),con(num_xtem,num_ytem))
              Allocate(aip(num_xtem,num_ytem),aim(num_xtem,num_ytem))
              Allocate(ajp(num_xtem,num_ytem),ajm(num_xtem,num_ytem))
              Allocate(apc(num_xtem,num_ytem))
              Do i=1,num_xtem
              Do j=1,num_ytem
              ap(i,j)=0
              con(i,j)=0
              aip(i,j)=0
              aim(i,j)=0
              ajp(i,j)=0
              ajm(i,j)=0
              apc(i,j)=0
              Enddo
              Enddo

              Do j=2,num_ytem-1
              aim(2,j)=ak1(2,j)/(0.5*hksi(2,j))
              aip(num_xtem-1,j)=ak1(num_xtem,j)/(0.5*hksi(num_xtem-1,j))
              Do i=2,num_xtem-2
              aip(i,j)=ak1(i+1,j)/(0.5*hksi(i,j)+0.5*hksi(i+1,j))
              aim(i+1,j)=aip(i,j)
              Enddo
              Enddo

              Do i=2,num_xtem-1
              ajm(i,2)=ae1(i,2)/(0.5*heta(i,2))
              ajp(i,num_ytem-1)=ae1(i,num_ytem)/(0.5*heta(i,num_ytem-1))
              Do j=2,num_ytem-2
              ajp(i,j)=ae1(i,j+1)/(0.5*heta(i,j)+0.5*heta(i,j+1))
              ajm(i,j+1)=ajp(i,j)
              Enddo
              Enddo

              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
              rmu=rmu_c(i,j)
              conduc0=conduc(i,j)
              If(id_ele(i,j).eq.11)conduc0=0
      apc(i,j)=-xmiu*epsil*omega*omega*vol(i,j)+
     1          aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)+
     1          1/r(i,j)/r(i,j)*vol(i,j)+
     1          rmu*xmiu*conduc0*omega*vol(i,j)
              Enddo
              Enddo
              Return
              End Subroutine
      End Module solverCalculate
      
      Program main
      Use meshElectroMagnetic
      Implicit None
      Integer(kind=4)::i,j
      Call readInMeshData()
      Call creatTemMesh()
      Write(*,*)nxx,nyy
      Write(*,'(2f6.1)')long_x,long_y
      Stop
      End Program
