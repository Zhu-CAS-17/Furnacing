      Program main
c      Use externalControl
      Use solverCalculate
      Use postElectroMagnetic
      Implicit None
      Integer(kind=4)::i,j
      Call readInControl
      Call readInFurnace
      Call readInCoil
      Call readInMaterial
      Call readInMeshData
      Call creatTemMesh
      Call electroCalculationMeshInfo
      Call fieldInit
      Call electroCofficience_const
      Call electroSourceTerm
      Do loop=1,10
      Write(*,100)loop
100   format(1x,'loop= ',i4)
c      Call fieldBoundaryCondition 
      Call gaussSeidelSolver
      Call electroSourceTerm
c      Call ifConvergency
c      If(convergency.eq.1)exit
      Enddo
      Write(*,*)'Program completed'
      Stop
      End Program
