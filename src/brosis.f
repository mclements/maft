cc data from brother-sister: including effects of PE on first pregnancy
c  feb14: adding the interaction effects
c R program: inter.r
c
c	parameter(npar=8,nn=15,m1=277,m2=9)
c	integer idata(m1,m2)
c	real*8 param(npar),correl(nn*(nn-1)/2)
c	real*8 loglik,prob(m1)
c	data param/-1.8,-1.5,1.,1.,.0,1.,2.,1./
c twosisf.dat is generated in twosisf.r
c	  open(1,file='c:/papers/twins/twosisf.dat')
c	  do i=1,m1
c	    read(1,*)(idata(i,j),j=1,m2)
c	  enddo
c	  
c	  n1 = 3
c	  n2 = 3
c	  n11 = 1
c	  n21= 1
c	  call bscov(param,npar,n11,n21,n1,n2,correl)
c	  print*,(correl(i),i=1,15)
c
c	  call bsll(param,npar,idata,m1,m2,loglik,prob)
c	  print*,'log likelihood = ',loglik
cc	open(1,file='c:/papers/twins/intprob.out')  
c	write(1,*)prob
c	end
c
c
c idata has m2=9 columns: 
c n11=      number of first preg from first sister
c pe11=     number among n11 with pe
c n12=      number of other pregnancies from first sister
c pe12=     number among n12 with pe
c n21, pe21, n22, pe22
c nfam =    numer of families with associated characteristics
c
 	subroutine bsll(param,npar,idata,m1,m2,loglik,prob)
	parameter (nn=15)
	integer idata(m1,m2), INFIN(nn)
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar),prob(m1)
	integer pe11,pe12,pe21,pe22


c default values for mvnorm call
	maxpts = 25000
        abseps = 0.001 
	releps = 0.

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
c 
	t1 = param(1)
	t2 = param(2)
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
c
c adding interaction terms
	  call bscov(param,npar,n11,n21,n1,n2,correl)
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
	subroutine bscov(param,npar,n11,n21,n1,n2,correl)
	parameter (nn=15)
	real*8 param(npar),cov(nn,nn),correl(*)
c workspace
	real Ag(nn,nn),Am(nn,nn),Af(nn,nn),As(nn,nn)
	real Am2(nn,nn), Af2(nn,nn),Abros(nn,nn),Abros2(nn,nn)

c parameters: 1=intercept, 3=genetic-maternal,
c 	      4=genetic-paternal, 5=common sibling
c             6=mother X late deliv 7= father X late delv
c             8= common husband-wife environment
	sigm2 = param(3)
	sigf2 = param(4)
	sigs2 = param(5)
	sigm22= param(6)
	sigf22= param(7)
        sigfam= param(8)
	sigfam2= param(9)
	nvar =  n1 + n2

c setup covariation matrix: 
	  do i=1,nvar
	    do j=1,nvar
              Am(i,j) = 1.
	      Af(i,j) = 1.
	      As(i,j) = 1.
	      Am2(i,j) = 0.
	      Af2(i,j) = 0.
	      Abros(i,j) = 0.
              Abros2(i,j) = 0.
	  enddo
	  enddo
c
c genetic-maternal and paternal, bro-sis correlation
c
	  do i=n1+1,nvar
	    do j=1,n1
		Am(i,j)=0.
		Af(j,i)=0.
                Abros(i,j) = 0.5
	    enddo
	  enddo
	  do i=1,n1
	    do j=n1+1,nvar
		Am(i,j)=0.
		Af(j,i)=0.
		Abros(i,j) = 0.5
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
	    Abros2(i,j) = 0.5
	  enddo
	enddo
	do i=(n1+1+n21),nvar
	  do j=1+n11,n1
	    Abros2(i,j) = 0.5
	  enddo
	enddo

c total variance
	do i =1,nvar
	  do j=1,nvar
	  cov(i,j)=sigm2*Am(i,j) + sigf2*Af(i,j) + sigs2*As(i,j) +
     &             sigm22*Am2(i,j) + sigf22*Af2(i,j)+
     &             sqrt(sigm2)*sqrt(sigf2)*Abros(i,j)+
     &             sqrt(sigm22)*sqrt(sigf22)*Abros2(i,j)+
     &             sigfam*Am(i,j) + sigfam2*Am2(i,j)
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
