      Module postElectroMagnetic
      Use materialElectroMagnetic
      Implicit None
      Contains
              Subroutine checkMagneticVectorPotential 
              Implicit None
              Integer(kind=4)::i,j
              Open(10,file='fa.dat')
              Write(10,50)
              Write(10,51)num_xtem,num_ytem
              Write(10,52)((y(i,j),i=1,num_xtem),j=1,num_ytem)
              Write(10,52)((x(i,j),i=1,num_xtem),j=1,num_ytem)
              Write(10,52)((far(i,j),i=1,num_xtem),j=1,num_ytem)
              Write(10,52)((fai(i,j),i=1,num_xtem),j=1,num_ytem)
   50 format(1x,'variables="x","y","far","fai"')
   51 format(1x,'zone t= "zone 1" , i=',i3,' , j=',i3,' , f= block')
   52 format(6(1pe12.5,1x)) 
              Close(10)
              Return
              End Subroutine
      End Module postElectroMagnetic 
