c by Qi-Sheng Chen Nov. 4, 1999
c Dept. of Mech. Eng., SUNY at Stony Brook, Stony Brook, NY 11794-2300
c email: qschen@pml.eng.sunysb.edu
c 
	subroutine spline(x,y,n,yp1,ypn,y2)
	implicit real*8 (a-h,o-z)
	integer n, nmax
	real*8 yp1,ypn,x(n),y(n),y2(n)
	parameter (nmax=500)
	integer i,k
	real*8 p,qn,sig,un,u(nmax)
	if(yp1.gt..99e30)then
	y2(1)=0
	u(1)=0
	else
	y2(1)=-0.5
	u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif
	do i=2,n-1
	sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	p=sig*y2(i-1)+2
	y2(i)=(sig-1.)/p
	u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1   /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
	enddo
	if(ypn.gt..99e30)then
	qn=0
	un=0
	else
	qn=0.5
	un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	do k=n-1,1,-1
	y2(k)=y2(k)*y2(k+1)+u(k)
	enddo
	return
	end

	subroutine splint(xa,ya,y2a,n,x,y)
	implicit real*8 (a-h,o-z)
	integer n
	real*8 x,y,xa(n),y2a(n),ya(n)
	integer k,khi,klo
	real*8 a,b,h
	klo=1
	khi=n
1	if(khi-klo.gt.1)then
	k=(khi+klo)/2
	if(xa(k).gt.x)then
	khi=k
	else
	klo=k
	endif
	goto 1
	endif
	h=xa(khi)-xa(klo)
      if(h.eq.0)then
      write(*,*)'bad xa input in splint'
      stop
      endif

	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	y=a*ya(klo)+b*ya(khi)+
     1    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6
	return
	end

