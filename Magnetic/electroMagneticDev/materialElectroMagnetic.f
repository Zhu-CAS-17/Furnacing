      Module materialElectroMagnetic
      Use meshElectroMagnetic
      Use fieldPhysical
      Use externalControl
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
c      Real(kind=8),Allocatable::con(:,:)
c con(:,:),current densit J.
      Contains
              Subroutine readInCoil
              Implicit None
              Open(10,file='coil.dat')
              Read(10,*)freq
              Read(10,*)coil_xx0,coil_yy0,coil_dxx,width_coil,xx_r,yy_r
              Close(10)
              omega=freq*2*3.1415926
              Write(*,*)'readInCoil from coil.dat'
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
              Write(*,*)'Call readInFurnace from geom.dat' 
              Return
              End Subroutine

              Subroutine readInMaterial
c conductivity: /ohm/m
c only consider conduction of sic and graphite heater
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
              Write(*,*)'readInMaterial from mater.dat'
              Return
              End Subroutine
c current density based on input.dat and coil.dat
c              Subroutine add_con(x0,y0)
c              Implicit None
c              Real(kind=8)::x0,y0,xxo,yyo,xxi,yyi,area
c              Real(kind=8)::dxdksi,dydksi,dxdeta,dydeta,xjacb
c              Integer(kind=4)::i,j,id,i1,j1
c              Integer(kind=4),Allocatable::iid(:),jjd(:)
c              Allocate(iid(num_xtem),jjd(num_xtem))
c              Allocate(con(num_xtem,num_ytem))
c              id=0
c              area=0.0
c              Do i=2,num_xtem-1
c              Do j=2,num_ytem-1
c              xxo=0.75*0.0254
c              xxi=(0.75-0.048)*0.0254
c              yyo=0.1875*0.0254
c              yyi=(0.1875-0.048)*0.0254
c              If(dabs(x(i,j)-x0).le.xxo.and.
c     1           dabs(y(i,j)-y0).le.yyo.and.
c     1           (dabs(x(i,j)-x0).ge.xxi.and.
c     1            dabs(y(i,j)-y0).ge.yyi))then
c              id=id+1
c              dxdeta=xc(i,j+1)-xc(i,j)
c              dydeta=yc(i,j+1)-yc(i,j)
c              dxdksi=x(i,j)-x(i-1,j)
c              dydksi=y(i,j)-y(i-1,j)
c              If(i.eq.2.or.i.eq.num_xtem)dxdksi=dxdksi*2.0
c              If(i.eq.2.or.i.eq.num_xtem)dydksi=dydksi*2.0
c              xjacb=dxdksi*dydeta-dxdeta*dydksi
c              area=area+xjacb
c              If(id.gt.6000)then
c              Write(*,*)'need iid with elements of',id
c              Stop 'error in add_con, check iid'
c              Endif
c              iid(id)=i
c              jjd(id)=j
c              Endif
c              Enddo
c              Enddo
c              Write(*,*)'cells in coil',id
c              Write(*,*)'testlast'
c
c              Do i=1,id
c              i1=iid(i)
c              j1=jjd(i)
c              con(i1,j1)=con(i1,j1)+xmiu*cden/area*vol(i,j)
c              Enddo
c
c              Return
c              End Subroutine

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
              
c material conductivity
              Real(kind=8) Function conduc(i,j)
c 1:graph 2:soft felt 3: sic powder 4:argon 5:hard feld 6: graphe foam
c 7:felt 8:temperature constant 9:crucible 10:air with conduction
c 11:copper 12:bottom graphite 13: crystal
c 14:gar with apparent conductivity 50w/m/k
              Implicit None
              Integer(kind=4)::i,j,id,iee
              Real(kind=8)::tcel,temp,ee
              conduc=0.0
              id=id_ele(i,j)
c              Write(*,*)'id',id
              If(id.eq.2.or.(id.ge.5.and.id.le.12))then
              conduc=conduc_d(id)
              Endif
              tcel=t(i,j)-273
              ee=(t(i,j)-273)/250.0+1
              iee=int(ee)
c              Write(*,*)'iee',iee
              If(iee.gt.9)iee=9
              If(iee.lt.1)iee=1
              If(id.eq.1.or.id.eq.8)then
              conduc=rinterp(con_grat,con_gra,size(con_gra),tcel)
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
              conduc=rinterp(con_sict,con_sic,size(con_sic),tcel)
              Endif
              Return
              End Function
              
              Real(kind=8) Function rinterp(x,y,n,x0)
              Implicit None
c              Real(kind=8),Allocatable::x(:),y(:)
              Integer(kind=4)::n,i,ii
              Real(kind=8)::x(n),y(n)
              Real(kind=8)::x0
c              Allocate(x(n),y(n))
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
