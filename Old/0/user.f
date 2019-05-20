      subroutine output_user
c///////////////////////////////////////////////////////////////////////////////
c   export the temperature in the crystal for thermal stress calculation
c       write(10,*)'variables ="x","y","u","v","T(K)",
c     1    "concentration","partial pressure","parti","p","pc","st"')
      include 'sc.h'
      ibot=1
	itop1=1
	itop2=1
	jleft=1
	jright=1
      
	do i=l1,1,-1
	if(x(i,1) .lt. (0.175-crys_len)) then
	ibot=i+1
      goto 800
	endif
	enddo
800   continue

	do i=1,l1
	if(x(i,1).GT.0.175) then
	itop1=i-1
	goto 900
	endif
	enddo
900   continue
	do i=1,l1
	if(x(i,1).GT.0.178) then
	itop2=i-1
	goto 1010
	endif
	enddo
1010  continue
      do j=1,m1
!modify parr	if(y(1,j).GT.0.025) then
      if(y(1,j).GT.0.03) then
	jright=j-1
	goto 1011
	endif
	enddo
1011  continue
      open(unit=10, file='cryst1.dat')
      rewind 10
      write(10,*)'title=','"zbzhang"'
      write(10,*)'variables ="i","r","z","T(K)"'
	write(10,*) itop1,ibot
      write(10,1012) jright-jleft+1,itop1-ibot+1
1012  format(1x,'zone t= "crystal",i=',i3,' , j=',i3,' , f= point')
c      write(10,*) x(ibot,1),x(itop,1),y(1,jleft),y(1,jright)
c      write(10,*) ibot,itop,jleft,jright
      ii=0
	do i=itop1,ibot,-1
	  do j=jleft,jright
	  ii=ii+1
	  write(10,1013) ii,y(i,j),x(i,j),t(i,j)
	  enddo
	enddo
	close(10)
1013   format(I4,2f10.4,f10.2)

      open(unit=10, file='cryst2.dat')
      rewind 10
      write(10,*)'title=','"zbzhang"'
      write(10,*)'variables ="i","r","z","T(K)"'
      write(10,1012) jright-jleft+1,itop2-ibot+1
c      write(10,*) x(ibot,1),x(itop,1),y(1,jleft),y(1,jright)
c      write(10,*) ibot,itop,jleft,jright
      ii=0
	do i=itop2,ibot,-1
	  do j=jleft,jright
	  ii=ii+1
	  write(10,1013) ii,y(i,j),x(i,j),t(i,j)
	  enddo
	enddo
	close(10)    

c z_pres is set in radi_ini
      i_rb1=interid(z_pres,x,l1)
      i_rb2=interid(z_pres2,x,l1)


50    continue

      tmax1=0
      i1=0
      j1=0
      do j1=1,j_rr-1
      do i1=i_rb1+1,i_rb2-1
      if(t(i1,j1).gt.tmax1.and.i_powder(i1,j1).eq.0)then
      tmax1=t(i1,j1)
      imaxp=i1
      jmaxp=j1
      endif
      enddo
      enddo

c temperature distribution in the growth chamber 
      ip1=1
	ip2=1
	ip3=1
      do i=1,l1
	if (x(i,2).gt.0.12 .and. id_ele(i,2).eq.4) then
	ip1=i
	goto 1014
	endif
	enddo
1014  continue
	
	do i=ip1,l1
	if (id_ele(i,2).ne.4) then
	ip2=i
	goto 1015
	endif
	enddo
1015  continue


	open(unit=10, file='chamber_temp.dat')
      rewind 10
      write(10,*)'title=','"zbzh"'
      write(10,*)'variables ="r","z","T(K)"'
      write(10,1018) jright,i_rt-imaxp+1
1018  format(1x,'zone t= "chamber",i=',i3,' , j=',i3,' , f= point')
      do i=imaxp,i_rt
	  do j=1,jright
	   write(10,1019) y(i,j),x(i,j),t(i,j)
	  enddo
	enddo
      close(10)
1019  format(2f10.4,f10.2)
      end
