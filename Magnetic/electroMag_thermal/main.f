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
