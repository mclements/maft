c June 13 2002:
c replacing father effect with fetal effect
c
c idata has m2=9 columns: 
c n11=      number of first preg from first sister
c pe11=     number among n11 with pe
c n12=      number of other pregnancies from first sister
c pe12=     number among n12 with pe
c n21, pe21, n22, pe22
c nfam =    numer of families with associated characteristics
c
 	subroutine ssfetal(param,abseps,npar,idata,m1,m2,loglik,prob)
	parameter (nn=15)
	integer idata(m1,m2), INFIN(nn)
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar),prob(m1)
	integer pe11,pe12,pe21,pe22


c default values for mvnorm call
	maxpts = 25000
	releps = 0.
	do i=1,nn
	  delta(i)=0.
	enddo

c	nvar=2
c	do i=1,nvar
c	  lower(i)=-1.5
c	  upper(i)=-1.5
c	  infin(i)=1
c	enddo
c	correl(1) = .2
c	call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
c    &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
c	print*,'first call to mvt=', value
c	return
c
c parameters: 1= intercept for pe1
c             2= intercept for pe2, 
c             3=genetic-maternal,
c 	      4=genetic father ==> fetal
c             5=common sibling
c
c 
	t1 = param(1)
	t2 = param(2)
	sigm2 = param(3)
	sigf2 = param(4)
	sigs2 = param(5)
c
	loglik = 0.
c
c process data: 
c
	do idat=1,m1
          prob(idat)=0.
	  n11  = idata(idat,1)
	  pe11  = idata(idat,2)
	  n12 = idata(idat,3)
	  pe12 = idata(idat,4)
	  n21 = idata(idat,5)
	  pe21= idata(idat,6)
	  n22= idata(idat,7)
	  pe22= idata(idat,8)
	  nfam= idata(idat,9)
	  n1= n11+n12
	  n2= n21+n22
          nvar = n1 + n2

	  do j=1,nvar
	    lower(j) = t2
	    upper(j) = t2
	  enddo
c if the first pregnancy exists
          if (n11 .eq. 1) then
	    lower(1) = t1
	    upper(1) = t1
	  endif
	  if (n21 .eq. 1) then
	    lower(n1+1)= t1
	    upper(n1+1)= t1
	  endif

	  call ssfetcov(param,npar,n11,n21,n1,n2,correl)
c          print*,(correl(i),i=1,nvar*(nvar-1)/2)
c
c start do loop
          do j=1,nvar
            infin(j) = 1
          enddo
c
c 
c from first sister
	if (n11 .eq. 0) then
	  do k=1,pe12
	    infin(k) = 0
	  enddo
	  mult1 = combi(n12,pe12)
	endif
	if (n11 .eq. 1) then
	  do k=1,pe11
   	       infin(k) = 0
	  enddo
	  do k=n11+1,n11+pe12
	    infin(k) = 0
	  enddo
	  mult1 = combi(n12,pe12)
	endif
c second sister
	if (n21 .eq. 0) then
	  do k=n1+1,n1+pe22
	    infin(k) = 0
	  enddo
	  mult2 = combi(n22,pe22)
	endif
	if (n21 .eq. 1) then
	  do k=n1+1,n1+pe21
   	       infin(k) = 0
	  enddo
	  do k=n1+n21+1,n1+n21+pe22
	    infin(k) = 0
	  enddo
	  mult2 = combi(n22,pe22)
	endif
c
c        print*,nvar,(lower(ii),ii=1,nvar),(upper(ii),ii=1,nvar)
        call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = mult1*mult2
	    prob(idat) = cmult * value
c	    print*,idat,((1-infin(ii)),ii=1,nvar),prob(idat),inform

        loglik = loglik + nfam * log(prob(idat))

c end idat
	enddo
	
	loglik = -loglik
	return
	end 


c computing covariance matrices, including interactions
c  between maternal and later deliveries 
c
c n1 = number of daughters from first father
c n2 = number of daughters from second father
c
	subroutine ssfetcov(param,npar,n11,n21,n1,n2,correl)
	parameter (nn=15)
	real*8 param(npar),cov(nn,nn),correl(*)
c workspace
	real Ag(nn,nn),Am(nn,nn),Af(nn,nn),As(nn,nn)
	real Am2(nn,nn), Af2(nn,nn),Afam(nn,nn),Afmat(nn,nn)

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
	    Af(i,j)=0.5
	    Af(j,i)=0.5
	    Afmat(i,j)=0.5
	    Afmat(j,i)=0.5
	  enddo
	enddo
        do i=n1+2,nvar
	  do j=n1+1,i-1
	     Af(i,j)=0.5
	     Afmat(i,j)=0.5
	     Af(j,i)=0.5
	     Afmat(j,i)=0.5
	  enddo
	enddo

c ......................          interaction terms
	do i=1+n11,n1
	  do j=1+n11,n1
	    Am2(i,j) = 1.
	    Af2(i,j) = 1.
	  enddo
	enddo
	do i=(n1+1+n21),nvar
	  do j=(n1+1+n21),nvar
	    Am2(i,j) = 1.
	    Af2(i,j) = 1.
	  enddo
	enddo
	do i=1+n11,n1
	  do j=(n1+1+n21),nvar
	    Am2(i,j) = 0.5
	  enddo
	enddo
	do i=(n1+1+n21),nvar
	  do j=1+n11,n1
	    Am2(i,j) = 0.5
	  enddo
	enddo
c parameters: 1,2=intercept, 3=genetic-maternal,
c 	      4=genetic-fetal, 5=common sibling,
c             6=family environment

	sigm2 = param(3)
	sfet2 = param(4)
	sigs2 = param(5)
        sigfam2= param(6)

c total variance
	do i =1,nvar
	  do j=1,nvar
	  cov(i,j)=sigm2*Am(i,j) + sfet2*Afmat(i,j) + 
     &             sigs2*As(i,j) + sigfam2*Afam(i,j) + 0.
c     &             sigm22*Am2(i,j) + sigf22*Af2(i,j)
	  enddo
	enddo

c	do i=1,nvar
c	  print*,(cov(i,j),j=1,nvar)
c	enddo
c	print*,''

c extract correlation matrix: lower triangular part
	  do i=2,nvar
	    do j=1,i-1
	      vari = cov(i,i) + 1.
	      varj = cov(j,j) + 1.
	      correl(j+((i-2)*(i-1))/2)=cov(i,j)/sqrt(vari)/sqrt(varj)
	    enddo
	  enddo

	return
	end

c computing combinations:
c	function combi(n,k)
c	combi=1.
c	if (k .eq. 0) return
c	if (k .eq. n) return
c	do i=1,k
c	  combi = combi*float(n-i+1)/float(i)
c	enddo
c	return
c	end
c
