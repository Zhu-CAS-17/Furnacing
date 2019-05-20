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
