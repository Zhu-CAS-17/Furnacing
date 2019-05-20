      Module solverCalculate
c      Use meshElectroMagnetic
      Use materialElectroMagnetic
c      Use fieldPhysical
      Implicit None
c magnetic potential A f: far Real;fai Image
c      Real(kind=8),Allocatable::far(:,:),fai(:,:)
      Real(kind=8),Allocatable::ap(:,:),aip(:,:),aim(:,:)
      Real(kind=8),Allocatable::ajp(:,:),ajm(:,:),apc(:,:)
c      Complex(kind=8),Allocatable::apc(:,:)
c old code apc is complex data
      Real(kind=8),Allocatable::con(:,:)
      Real(kind=8),Allocatable::br(:,:),bi(:,:)
      Real(kind=8)::tempfar,tempfai,deltfar,deltfai
      Integer(kind=4)::loop=1,convergency=0
      Real(kind=8)::residual=1.0e-6,relaerr=1.0e-6,loopmax=10e6
      Contains
              Subroutine electroCofficience_const
              Implicit None
              Integer(kind=4)::i,j
              Real(kind=8)rmu,conduc0,xx0,yy0
              Allocate(ap(num_xtem,num_ytem))
              Allocate(aip(num_xtem,num_ytem),aim(num_xtem,num_ytem))
              Allocate(ajp(num_xtem,num_ytem),ajm(num_xtem,num_ytem))
              Allocate(apc(num_xtem,num_ytem),con(num_xtem,num_ytem))
              Allocate(br(num_xtem,num_ytem),bi(num_xtem,num_ytem))
              Do i=1,num_xtem
              Do j=1,num_ytem
              ap(i,j)=0
              con(i,j)=0
              aip(i,j)=0
              aim(i,j)=0
              ajp(i,j)=0
              ajm(i,j)=0
              apc(i,j)=0
              br(i,j)=0
              bi(i,j)=0
              Enddo
              Enddo
c domain control volume
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
c coil position based on coil.dat
              Do i=1,5
              xx0=coil_xx0+coil_dxx*(i-1)
              yy0=coil_yy0
              Call add_con(xx0,yy0)
              Enddo

              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
c              rmu=rmu_c(i,j)
              conduc0=conduc(i,j)
c id_ele=11 copper
              If(id_ele(i,j).eq.11)conduc0=0
c if using complex formate, use the following method
c      apc(i,j)=-xmiu*epsil*omega*omega*vol(i,j)+
c     1          aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)+
c     1          1/r(i,j)/r(i,j)*vol(i,j)+
c     1          rmu*xmiu*conduc0*omega*vol(i,j)
      apc(i,j)=ajp(i,j)+ajm(i,j)+aip(i,j)+aim(i,j)
     1           +vol(i,j)/r(i,j)/r(i,j)
     1           -xmiu*epsil*omega*omega*vol(i,j)
c assumption: rmu=1.0
c      br(i,j)=1.0*xmiu*omega*conduc0*vol(i,j)*fai(i,j)+con(i,j)
c      bi(i,j)=1.0*xmiu*omega*conduc0*vol(i,j)*far(i,j)
c examine conduc0
              Enddo
              Enddo
              Write(*,*)'Call electroCofficience_const'
c chect the fai(i,j) and far(i,j)
              Open(10,file='coefmesh1.dat')
              Write(10,90)
              Write(10,91)num_xtem,num_ytem
              Do j=1,num_ytem
              Do i=1,num_xtem
              Write(10,'(4es12.4)')y(i,j),x(i,j),apc(i,j)
              Enddo
              Enddo
c              Write(10,92)((y(i,j),i=1,num_xtem),j=1,num_ytem)
c              Write(10,92)((x(i,j),i=1,num_xtem),j=1,num_ytem)
c              Write(10,92)((apc(i,j),i=1,num_xtem),j=1,num_ytem)
   90 format(1x,'variables="r","z","apcr","apci"')
   91 format(1x,'zone t= "zone 1",i=',i3,' , j=',i3,' , f= block')
   92 format(6(1pe12.4,1x))
              Close(10)
c              Stop
              Return
              End Subroutine electroCofficience_const
              
              Subroutine electroSourceTerm
              Implicit None
              Integer(kind=4)::i,j
              Real(kind=8)::rmu,conduc0
c              Allocate(br(num_xtem,num_ytem),bi(num_xtem,num_ytem))
              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
              rmu=rmu_c(i,j)
              conduc0=conduc(i,j)
c id_ele=11 copper,ignoring the eddy current in copper
              If(id_ele(i,j).eq.11)conduc0=0
c rmu=1.0
      br(i,j)=1.0*xmiu*omega*conduc0*vol(i,j)*fai(i,j)+con(i,j)
      bi(i,j)=1.0*xmiu*omega*conduc0*vol(i,j)*far(i,j)
              Enddo
              Enddo
              Return
              End Subroutine
c current density based on input.dat and coil.dat
c current density J equals I/A
              Subroutine add_con(x0,y0)
              Implicit None
              Real(kind=8)::x0,y0,xxo,yyo,xxi,yyi,sumarea
              Real(kind=8)::dxdksi,dydksi,dxdeta,dydeta,xjacb
              Integer(kind=4)::i,j,id,i1,j1
              Integer(kind=4),Allocatable::iid(:),jjd(:)
c              Real(kind=8),Allocatable::area(:)
              Character(len=20) :: nfile
              Integer(kind=4)::counts=0
              Allocate(iid(num_xtem),jjd(num_xtem))
c              Allocate(area(num_xtem*num_ytem/20))
              id=0
              sumarea=0.0
              xxo=0.75*0.0254
              xxi=(0.75-0.048)*0.0254
              yyo=0.1875*0.0254
              yyi=(0.1875-0.048)*0.0254
c              xxo=xx_r
c              yyo=yy_r
c              xxi=xx_r-width_coil
c              yyi=yy_r-width_coil
              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
              If((dabs(x(i,j)-x0).le.xxo.and.
     1           dabs(y(i,j)-y0).le.yyo).and.
     1           (dabs(x(i,j)-x0).ge.xxi.and.
     1            dabs(y(i,j)-y0).ge.yyi))then
c the number of grids involved in ccoper coil
              id=id+1
              dxdeta=xc(i,j+1)-xc(i,j)
              dydeta=yc(i,j+1)-yc(i,j)
              dxdksi=x(i,j)-x(i-1,j)
              dydksi=y(i,j)-y(i-1,j)
              If(i.eq.2.or.i.eq.num_xtem)dxdksi=dxdksi*2.0
              If(i.eq.2.or.i.eq.num_xtem)dydksi=dydksi*2.0
c area of grid
              xjacb=dxdksi*dydeta-dxdeta*dydksi
              sumarea=sumarea+xjacb
              If(id.gt.6000)then
              Write(*,*)'need iid with elements of',id
              Stop 'error in add_con, check iid'
              Endif
c              area(id)=xjacb
              iid(id)=i
              jjd(id)=j
              Endif
              Enddo
              Enddo
              Write(*,*)'cells in coil',id
              Write(*,*)sumarea

c              counts=counts+1
c              Write(nfile,*)counts
c              nfile='nodes'//trim(adjustl(nfile))//'.dat'
c              Open(10,file=nfile)
c              Write(10,1001)
c              Write(10,1002)num_xtem,num_ytem
              Do i=1,id
              i1=iid(i)
              j1=jjd(i)
c multiply the vol(i1,j1) denotes the control volume
              con(i1,j1)=con(i1,j1)+cden/sumarea*vol(i1,j1)*xmiu
              Enddo
c              Write(10,1003)((y(i,j),i=1,num_xtem),j=1,num_ytem)
c              Write(10,1003)((x(i,j),i=1,num_xtem),j=1,num_ytem)
c              Write(10,1003)((con(i,j),i=1,num_xtem),j=1,num_ytem)
c1001  format(1x,'variables="x","y","con"')
c1002  format(1x,'zone t="zone 1", i=',i3,' , j=',i3,' , f=block')
c1003  format(6(1pe12.5,1x))
c              Close(10) 
cc             stop
              Return
              End Subroutine

              Subroutine ifConvergency
              Implicit None
c              Real(kind=8)::residual=1.0e-10
      If((deltfar.le.residual.and.deltfai.le.residual).or.
     1   (loop.ge.loopmax)) convergency=1
              Write(*,*)convergency
              Return
              End Subroutine  
              
              Subroutine gaussSeidelSolver
              Implicit None
              Integer(kind=4)::i,j
              Real(kind=8)::reladeltr,reladelti
              Do i=2,num_xtem-1
              Do j=2,num_ytem-1
      tempfar=1.0/apc(i,j)*(aip(i,j)*far(i+1,j)+aim(i,j)*far(i-1,j)+
     1                      ajp(i,j)*far(i,j+1)+ajm(i,j)*far(i,j-1)+
     1                      br(i,j))
              deltfar=tempfar-far(i,j)
              reladeltr=abs(deltfar)/(far(i,j)+1.0e-10)
              
      tempfai=1.0/apc(i,j)*(aip(i,j)*fai(i+1,j)+aim(i,j)*fai(i-1,j)+
     1                      ajp(i,j)*fai(i,j+1)+ajm(i,j)*fai(i,j-1)+
     1                      bi(i,j))
              deltfai=tempfai-fai(i,j)
              reladelti=abs(deltfai)/(fai(i,j)+1.0-10)
              far(i,j)=tempfar
              fai(i,j)=tempfai
              relaerr=max(reladeltr,reladelti)
              If(deltfar.gt.1.0e-10.or.deltfai.gt.1.0e-10)then
              Write(*,'(3es16.5)')deltfar,deltfai,relaerr
              else
              EndIf
c              If(deltfar.gt.1.0e-4.or.deltfai.gt.1.0e-4)then
c              Write(*,*)'deltafar=', deltfar, 'deltafi=', deltfai
c              Endif
              Enddo
              Enddo
              Call fieldBoundaryCondition 
c              Call ifConvergency
              Call electroSourceTerm
              Return
              End Subroutine
      End Module solverCalculate
