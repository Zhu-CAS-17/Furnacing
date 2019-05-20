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
