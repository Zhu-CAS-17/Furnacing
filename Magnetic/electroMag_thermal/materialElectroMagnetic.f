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
