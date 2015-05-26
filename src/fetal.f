
	parameter (nn=15,n1=3,n2=3,npar=9)
	real*8 param(npar),cov(nn,nn)
c workspace
	real Ag(nn,nn),Am(nn,nn),Af(nn,nn),As(nn,nn),Afam(nn,nn)
	real Am2(nn,nn), Af2(nn,nn),Afmat(nn,nn)

	nvar =  n1 + n2

c setup covariation matrix: 
	  do i=1,nvar
	    do j=1,nvar
              Am(i,j) = 1.
	      Af(i,j) = 0.
	      Afam(i,j)=1.
	      As(i,j) = 1.
	      Am2(i,j) = 0.
	      Af2(i,j) = 0.
	  enddo
	  enddo
c
c genetic-maternal, fetal-mother and family env
c
	  do i=n1+1,nvar
	    do j=1,n1
		Am(i,j)=0.5
		Afam(j,i)=0.
		Afmat(i,j)=0.125
	    enddo
	  enddo
	  do i=1,n1
	    do j=n1+1,nvar
		Am(i,j)=0.5
		Afam(j,i)=0.
		Afmat(i,j)=0.125
	    enddo
	  enddo
c genetic fetal-paternal, fetal-maternal
	do i=1,nvar
	    Af(i,i) = 1.
	    Afmat(i,i)=1.
	enddo
        do i=2,n1
	  do j=1,i-1
	    Af(i,j)=0.25
	    Af(j,i)=0.25
	    Afmat(i,j)=0.25
	    Afmat(j,i)=0.25
	  enddo
	enddo
        do i=n1+2,nvar
	  do j=n1+1,i-1
	     Af(i,j)=0.25
	     Afmat(i,j)=0.25
	     Af(j,i)=0.25
	     Afmat(j,i)=0.25
	  enddo
	enddo

	do i=1,nvar
	  write(6,*)(Afam(i,j),j=1,nvar)
	enddo
	print*
	do i=1,nvar
	  write(6,*)(Am(i,j),j=1,nvar)
	enddo
	end
