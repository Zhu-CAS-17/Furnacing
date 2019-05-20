      Module fieldPhysical
      Use meshElectroMagnetic
      Implicit None
c magnetic potential A f: far Real;fai Image
      Real(kind=8),Allocatable::t(:,:),far(:,:),fai(:,:)
      Contains
              Subroutine fieldInit
              Implicit None
              Integer(kind=4)::i,j
              Allocate(far(num_xtem,num_ytem),fai(num_xtem,num_ytem))
              Allocate(t(num_xtem,num_ytem))
              Do i=1,num_xtem
              Do j=1,num_ytem
              far(i,j)=0
              fai(i,j)=0
              t(i,j)=293
              Enddo
              Enddo
              Write(*,*)'Call fieldInit'
              Return
              End Subroutine

              Subroutine fieldBoundaryCondition
              Implicit None
              Integer(kind=4)::i,j
              Do i=1,num_xtem
c axisymmetric
              far(i,1)=0
              fai(i,1)=0
c Robin condition//side wall
c there are two corner points, adapting Robin condition first.
           far(i,num_ytem)=far(i,num_ytem-1)/(1+0.5*heta(i,num_ytem-1))
           fai(i,num_ytem)=fai(i,num_ytem-1)/(1+0.5*heta(i,num_ytem-1))
              Enddo
c axial direction gradient equals 0
              Do j=1,num_ytem
              far(1,j)=far(2,j)
              fai(1,j)=fai(2,j)
              far(num_xtem,j)=far(num_xtem-1,j)
              fai(num_xtem,j)=fai(num_xtem-1,j)
              Enddo
              Return
              End Subroutine
      End Module fieldPhysical
