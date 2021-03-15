module INparameter
    implicit none
    
    real*8,parameter ::pi = 3.1415926d0
    integer,parameter ::N = 100000 !ѭ������, loop times
   
    integer,parameter ::DN = 61 !�ײ㵥����ԭ�����������������1����֤�Գ�, number of atoms in the bottom layer in one direction (should be an odd number)
    integer,parameter ::ddN = 2 !��ʽԭ����, number of atoms belonging to a unit cell of the bottom layer
    integer,parameter ::dnumtot =  ddN*DN*DN    
    integer,parameter ::UN = 31 !���㵥����ԭ����,number of atoms in the top layer in one direction 
    integer,parameter ::uuN = 1 !��ʽԭ����, number of atoms belonging to a unit cell of the top layer
    integer,parameter ::unumtot =  uuN*uN*uN
    real*8,parameter ::Rvdw1 = 2.1d0  !�ײ�ԭ�ӵķ����߶�˹�뾶, van der Walls radius of atoms in the bottom layer
    real*8,parameter ::Rvdw2 = 2.1d0  !����ԭ�ӵķ����߶�˹�뾶, van der Walls radius of atoms in the top layer
    real*8,parameter ::AA = 1.01d0*(Rvdw1+Rvdw2) !ʵ���жϾ���, criterion for bonding
    real*8,parameter ::DD =((AA)*(AA)-1.d0*(Rvdw1+Rvdw2)*1.d0*(Rvdw1+Rvdw2))   !Լ���жϾ����ƽ��
    
    real*8,parameter ::R11x= 4.04958d0  !����1 �����ʸx����, basis of the bottom layer (a1, x component)
    real*8,parameter ::R11y = 0.d0  !����1 �����ʸy����, basis of the bottom layer (a1, y component)
    real*8,parameter ::R12x = 0.d0  !����1 �����ʸx����, basis of the bottom layer (b1, x component)
    real*8,parameter ::R12y = 4.04958d0  !����1 �����ʸy����, basis of the bottom layer (b1, y component)
    real*8,parameter ::R21x = 4.76020d0  !����2 �����ʸx����, basis of the top layer (a2, x component)
    real*8,parameter ::R21y = 0.d0  !����2 �����ʸy����, basis of the top layer (a2, y component)
    real*8,parameter ::R22x = 2.3801d0  !����2 �����ʸx����, basis of the top layer (b2, x component)
    real*8,parameter ::R22y = 4.12245413d0  !����2 �����ʸy����,, basis of the top layer (b2, y component)
    
    real*8,parameter ::Rax = 2.02479d0  !����1 ��ʽ���, displacement between two atoms in a unit cell of bottom layer (x component)  
    real*8,parameter ::Ray = 2.02479d0  !����1,   displacement between two atoms in a unit cell of bottom layer (y component)  
    real*8,parameter ::Rbx = 0.d0  !����2 ��ʽ���,  displacement between two atoms in a unit cell of top layer (x component)  
    real*8,parameter ::Rby = 0.d0  !����2, displacement between two atoms in a unit cell of top layer (x component)  
    
end module INparameter

program InterfaceN
	use INparameter
    integer :: i, j, k, m, km, flagd, flagu, flag, sitep
    integer :: num, idum
    real*8 :: xx, yy, dX, dY, RR, cphi, sphi, NB, NBtot, SX ,SY
    real*8 :: PN1(2*dnumtot), PN2(2*unumtot), PN3(2*unumtot)
    real*8 :: R11(2),R12(2),R21(2),R22(2),Ra(2),Rb(2),Rc(2)
    real*8 ::SP(2), thetaR(unumtot)
    
    OPEN (21,file="translation.txt") !ƽ����, translation displacement
    OPEN (24,file="rotate.txt") !��ת��, rotation displacement
    OPEN (27,file="matching.txt") !ƥ������ƥ���ʣ�ƽ��ƥ����  , matching
    !OPEN (3, file="other.txt") !��֤ԭ��
    
    idum =1
    NBtot=0
	R11(1) = R11x
	R11(2) = R11y
	R12(1) = R12x
	R12(2) = R12y
	R21(1) = R21x
	R21(2) = R21y
	R22(1) = R22x
	R22(2) = R22y
	Ra(1) = Rax
	Ra(2) = Ray
	Rb(1) = Rbx
	Rb(2) = Rby
	
	!�ײ�ԭ��λ��, bottom layer
	call create(PN1,ddN,DN,R11,R12)
	!��ʽ���
	flagd = 1
	do while (ddN-flagd)
		flagd = flagd +1
		call multi(PN1,ddN,DN,Ra,flagd)
	end do
	
	!��ʼ����ԭ��λ��, top layer
	call create(PN2,uuN,UN,R21,R22)
	!��ʽ���
	flagu = 1
	do while (uuN-flagu)
		flagu = flagu +1
		call multi(PN2,uuN,UN,Rb,flagu)
	end do
	
	!��ʼ�������, random displacement
	do m =1, N
		
		!���ö���ԭ��λ��
		do km = 1, 2*unumtot 
			PN3(km) = PN2(km)
		end do
		
		!����任		
		Rc=0.d0
		thetaR=0.d0
		xx=0.d0
		yy=0.d0
		cphi=0.d0
		sphi=0.d0
		SX=0.d0
		SY=0.d0
		!ƽ��, translation
		xx = ran1(idum) !X
		yy = ran1(idum) !Y
		Rc(1) = R11(1)*(xx-0.5d0)+R12(1)*(xx-0.5d0)
		Rc(2) = R11(2)*(yy-0.5d0)+R12(2)*(yy-0.5d0)
		!��ת���ɿ���, rotation
		cphi=1.d0-2.d0*ran1(idum) !1-2*rand(0,1) 
		if (ran1(idum) .GT. 0.5) then
			sphi=sqrt(1.d0-cphi*cphi)
		else 
			sphi=-sqrt(1.d0-cphi*cphi)
		end if
		!��ʼ�任
		do k = 1, (unumtot)
			!SP(1) = PN3(2*(uuN*(UN+1)*(UN+1))-1) !��ת���ģ���ʵ����ԭ��
			!SP(2) = PN3(2*(uuN*(UN+1)*(UN+1))) !Ϊ��֤���μ���
			SX = PN3(2*k-1)
			SY = PN3(2*k)
			PN3(2*k-1) = SX*cphi - SY*sphi + Rc(1)
			PN3(2*k) = SX*sphi + SY*cphi + Rc(2) 
		end do
		
		!�ж�ƥ��, matching		
		dX=0.d0
		dY=0.d0
		RR=0.d0
		num=0		
		do i=1,(dnumtot)
			do j=1,(unumtot)
				dX = PN1(2*(i-1)+1) - PN3(2*(j-1) +1)
				dY = PN1(2*(i-1)+2) - PN3(2*(j-1) +2)
					if ((dX*dX) .LE. DD .AND. (dY*dY) .LE. DD) then
						RR = dX*dX + dY*dY
						if (RR .LE. DD) then
							num = num +1
						end if 
					end if
			end do
		end do
		
		!���, output
		NB = (1.d0*num)/(1.d0*unumtot)
		NBtot = NBtot + NB
		write (21,*) Rc(1), Rc(2)  
		write (24,*) cphi, sphi 
		write (27,*) num, NB, (NBtot/(1.d0*m)) 
		!write (3,*) SP(1), SP(2) !��֤�õ��Ƿ�Ϊԭ��
	end do
	
	
end program InterfaceN

!���ɳ�ʼԭ��λ��, create atoms
subroutine create (PN,ccN,CN,R1,R2)
	use INparameter
	integer :: CN, ccN, mm, mn , flagc, ii, jj
	real*8 :: PN(2*ccN*CN*CN)
	real*8 :: R1(2),R2(2)

    flagc =0
    do mm = 0,CN-1
		do mn = 0,CN-1
			flagc = flagc+1
			ii = mm-(CN-1)/2
			jj = mn-(CN-1)/2
			PN(2*flagc-1) = ii*R1(1) + jj*R2(1)
			PN(2*flagc) = ii*R1(2) + jj*R2(2)
		end do 
	end do
	
end subroutine create

!��ʽ���
subroutine multi (PN,ccN,CN,Rxy,site)
	use INparameter
	integer :: CN, mm, flagm, ii, jj, ccN
	real*8 :: PN(2*ccN*CN*CN)
	real*8 :: Rxy(2*ccN-2)
	integer :: site
	
	flagm = (site-1)*CN*CN
	do mm = 1, ((site-1)*CN*CN)
		flagm =flagm +1
		PN(2*flagm-1) = PN(2*mm-1) + Rxy(2*(site-1)-1)
		PN(2*flagm) = PN(2*mm) + Rxy(2*(site-1))
	end do 
	
end subroutine multi

!��������� rand(0,1) ����idum Ϊ����ֻ�����(0,1)
real*8 function ran1(idum)
	implicit none
	integer,intent(inout) ::idum
	integer ::IA,IM,IQ,IR,NTAB,NDIV
	real*8 ::AM,EPS,RNMX
	parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1.0+(IM-1.0)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	integer ::j,k,iv(NTAB),iy 
	data iv /NTAB*0/, iy /0/
	if (idum .le. 0 .or. iy .eq. 0) then					
		idum=max(-idum,1)						
		do j=NTAB+8,1,-1							
            k=idum/IQ
			idum=IA*(idum-k*IQ)-IR*k
			if(idum.lt.0) idum=idum+IM
			if(j.le.NTAB) iv(j)=idum
		enddo 
		iy=iv(1)
	endif
	k=idum/IQ										
	idum=IA*(idum-k*IQ)-IR*k						
	if (idum.lt.0) idum=idum+IM						
	j=1+iy/NDIV										
	iy=iv(j)										
	iv(j)=idum										
	ran1=min(AM*iy,RNMX)							
end function ran1
