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
