      Module externalControl
      Implicit None
      Real(kind=8)::init,dtm,lastt
      Integer(kind=4)::idesign,icontrol
      Real(kind=8)::pres_sys,cden0,crys_len,coil_move,felt_incr
      Real(kind=8)::t_ref,cden
      Real(kind=8),Allocatable::relax(:)
      Integer(kind=4),Allocatable::nsolve(:),ntimes(:)
      Contains
             Subroutine readInControl
             Implicit None
             Integer(kind=4)::i
cnsolve for u v p pc t ccen
             Allocate(relax(6),nsolve(6),ntimes(6))
             Open(1,file='input.dat')
             Read(1,*)init,dtm,lastt
             Read(1,*)idesign,pres_sys,cden0,crys_len
             Read(1,*)coil_move,felt_incr
             Read(1,*)icontrol
             Read(1,*)t_ref

             Read(1,*)(relax(i),i=1,6)
             Read(1,*)(nsolve(i),i=1,6)
             Read(1,*)(ntimes(i),i=1,6)
             Close(1)
             cden=cden0

             Write(*,*)'readInExternal input.dat'
             Return
             End Subroutine
      End Module externalControl 
