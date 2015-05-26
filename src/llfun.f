
c minimum main
c	parameter (nvar=2,npar=3,m1=4,m2=4)
c	real*8 correl(nvar*(nvar-1)/2),param(npar)
c input output
c	integer idata(m1,m2)
c	real*8 loglik,prob(m1)
c workspace for global
c	data (idata(1,j),j=1,4)/1,2,0,603/
c	data (idata(2,j),j=1,4)/99,2,0,603/
c	data (idata(3,j),j=1,4)/99,2,1,267/
c	data (idata(4,j),j=1,4)/99,2,2,130/
c
c	param(1) = -1.
c	param(2) = 1.
c	param(3) = 0.5
c	nparent=0

c	call covar(param,npar,nvar,nparent,correl)
c	print*,correl

c	call fun1(param,npar,idata,m1,m2,loglik,prob)
c	  print*,'loglik=',loglik
c	  print*,'prob=',prob
c	  print*,'sumprob= ', (prob(1)+prob(2)+prob(3))
c	end
c
c
c computing likelihood: long data structure: allowed unknown mothers
c   using covar2() with 4 variance components
c
	subroutine fun2(param,npar,idata,m1,m2,loglik,prob)
 	parameter (nn=15)
	integer idata(m1,m2),m1,m2,npar
	INTEGER ND, INFIN(nn), MAXPTS, INFORM, IVLS
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar),prob(m1)

c default values for mvnorm call
	maxpts = 25000
        abseps = 0.001 
	releps = 0.
c
	beta0 = param(1)
	sigg2 = param(2)
	sigm2 = param(3)
	sigf2 = param(4)
	sigs2 = param(5)
	tlim = beta0/sqrt(sigg2 + sigm2 + sigf2 + sigs2 + 1.)

	loglik = 0.
c
c process data
c
	do idat=1,m1
	  nparent=1
	  if (idata(idat,1) .eq. 99) nparent = 0
	  numdau  = idata(idat,2)
          nvar = nparent + numdau
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covar2(param,npar,nparent,numdau,correl)
c
 	  meverpe = idata(idat,1)
	  numdaupe = idata(idat,3)
	    if (meverpe .eq. 1) infin(1) = 0
	    do k=1+nparent,nparent+numdaupe
   	       infin(k) = 0
	    enddo
c	    print*,(infin(ii),ii=1,nvar)
	    call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = combi(numdau,numdaupe)
	    prob(idat) = cmult * value
	    loglik = loglik + idata(idat,4) * log(prob(idat))
	  

c end idat
	enddo
	
	loglik = -loglik
	return
	end 


c
c computing likelihood: long data structure: allowed unknown mothers
c
	subroutine fun1(param,npar,idata,m1,m2,loglik,prob)
 	parameter (nn=15)
	integer idata(m1,m2),m1,m2,npar
	INTEGER ND, INFIN(nn), MAXPTS, INFORM, IVLS
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar),prob(m1)

c default values for mvnorm call
	maxpts = 25000
        abseps = 0.001 
	releps = 0.
c
	beta0 = param(1)
	sigb2 = param(2)
	sigc2 = param(3)
	tlim = beta0/sqrt(sigb2 + sigc2 + 1.)

	loglik = 0.
c
c process data
c
	do idat=1,m1
	  nparent=1
	  if (idata(idat,1) .eq. 99) nparent = 0
	  numdau  = idata(idat,2)
          nvar = nparent + numdau
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covar(param,npar,nvar,nparent,correl)
c
 	  meverpe = idata(idat,1)
	  numdaupe = idata(idat,3)
	    if (meverpe .eq. 1) infin(1) = 0
	    do k=1+nparent,nparent+numdaupe
   	       infin(k) = 0
	    enddo
c	    print*,(infin(ii),ii=1,nvar)
	    call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = combi(numdau,numdaupe)
	    prob(idat) = cmult * value
	    loglik = loglik + idata(idat,4) * log(prob(idat))
	  

c end idat
	enddo
	
	loglik = -loglik
	return
	end 

c computing likelihood: long data structure:
c
	subroutine fun0(param,loglik,npar,idata,m1,m2)
 	parameter (nn=15)
	integer idata(m1,m2),m1,m2,npar
	INTEGER ND, INFIN(nn), MAXPTS, INFORM, IVLS
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar)

c default values for mvnorm call
	maxpts = 25000
        abseps = 0.00001 
	releps = 0.
c
	beta0 = param(1)
	sigb2 = param(2)
	sigc2 = param(3)
	tlim = beta0/sqrt(sigb2 + sigc2 + 1.)

	loglik = 0.
c
c process data
c
	do idat=1,m1
	  if (idata(idat,1) .lt. 99) nparent = 1
	  if (idata(idat,1) .eq. 99) nparent = 0
	  numdau  = idata(idat,2)
          nvar = nparent + numdau
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covar(param,npar,nvar,nparent,correl)
c
 	  meverpe = idata(idat,1)
	  numdaupe = idata(idat,3)
	  if (nparent .eq. 1) then
	    if (meverpe .eq. 1) infin(1) = 0
	    do k=1+nparent,1+numdaupe
   	       infin(k) = 0
	    enddo
	    call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = combi(numdau,numdaupe)
	    prob = cmult * value
	    loglik = loglik + idata(idat,4) * log(prob)
	  endif

c end idat
	enddo
	
	loglik = -loglik
	return
	end 




	subroutine fun(param,loglik,npar,idata,m1,m2)
 	parameter (nn=15)
	integer idata(m1,m2),m1,m2,npar
	INTEGER ND, INFIN(nn), MAXPTS, INFORM, IVLS
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar)

c default values for mvnorm call
	maxpts = 25000
        abseps = 0.001 
	releps = 0.
c
	beta0 = param(1)
	sigb2 = param(2)
	sigc2 = param(3)
	tlim = beta0/sqrt(sigb2 + sigc2 + 1.)

	loglik = 0.
c process data
	do idat=1,m1
	  nparent = idata(idat,1)
	  nchild  = idata(idat,2)
 	  np = idata(idat,3)
      nvar = nparent + nchild
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covar(param,npar,nvar,nparent,correl)

	  if (nparent .eq. 1) then
	    if (np .eq. 1) infin(1) = 0
	    do j=0,nchild
		  do k=2,1+j
		    infin(k) = 0
		  enddo
		  call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                   MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
		  cmult = combi(nchild,j)
		  prob = cmult * value
		  loglik = loglik + idata(idat,4+j) * log(prob)
	     enddo
	  endif

c end idat
	enddo
	
	loglik = -loglik
	return
	end



c
c computing covariance matrices and extracting correlation
c
	subroutine covar(param,npar,nvar,nparent,correl)
	parameter (nn=15)
	real*8 param(npar),cov(nn,nn),correl(*)
c workspace
	real Ac(nn,nn),Ag(nn,nn)

	beta0 = param(1)
	sigb2 = param(2)
	sigc2 = param(3)

c setup covariation matrix: one mother and several fullsibs
	  do i=1,nvar
 	    do j=1,nvar
	      Ag(i,j)=0.
	      Ac(i,j)=0.
	    enddo
	  enddo

	  do i=1,nvar
	    Ag(i,i) = 1.
	    Ac(i,i) = 1.
	    do j=i+1,nvar
	      Ag(i,j)=0.5
	      Ag(j,i)=0.5
	    enddo
	  enddo
	if (nparent .eq. 1) then
	  do i=2,nvar
	    do j=i+1,nvar
	      Ac(i,j)= 1.
	      Ac(j,i)= 1.
	    enddo
	   enddo
	endif

	if (nparent .eq. 0) then
	  do i=1,nvar
	  do j=1,nvar
	    Ac(i,j) = 1.
	  enddo
	  enddo
	endif

c total variance
	do i =1,nvar
	  do j=1,nvar
	    cov(i,j) = sigb2*Ag(i,j) + sigc2*Ac(i,j)
	  enddo
	enddo

c	do i=1,nvar
c	  print*,(cov(i,j),j=1,nvar)
c	enddo
c	print*,''

c extract correlation matrix: lower triangular part
	var = sigb2 + sigc2 + 1.
	  do i=2,nvar
	    do j=1,i
	     correl(j + ((i-2)*(i-1))/2) = cov(i,j)/var
	    enddo
	  enddo

	return
	end
c
c using 4 variance components:
c
	subroutine covar2(param,npar,nparent,n1,correl)
	parameter (nn=15)
	real*8 param(npar),cov(nn,nn),correl(*)
c workspace
	real Ag(nn,nn),Am(nn,nn),Af(nn,nn),As(nn,nn)

c parameters: 1=intercept, 2=genetic-fetal, 3=genetic-maternal,
c 	      4=common family, 5=common sibling
	beta0 = param(1)
	sigg2 = param(2)
	sigm2 = param(3)
	sigf2 = param(4)
	sigs2 = param(5)
	nvar = nparent + n1 

c setup covariation matrix: 
	  do i=1,nvar
 	    do j=1,nvar
	      Ag(i,j)=0.5
	      Am(i,j)=0.25
	      Af(i,j)=1.
	      As(i,j)=1.
	    enddo
	    Ag(i,i) = 1.
	    Am(i,i) = 1.
	  enddo
c
	if (nparent .eq. 1) then
	  do i=1+nparent,nvar
c genetic-maternal
	     Am(1,i) = 0.5
	     Am(i,1) = 0.5
c common siblings
	     As(i,1) = 0.
	     As(1,i) = 0.
	  enddo
	endif


c total variance
	do i =1,nvar
	  do j=1,nvar
	  cov(i,j)=sigg2*Ag(i,j) + sigm2*Am(i,j) + 
     &             sigf2*Af(i,j) + sigs2*As(i,j)
	  enddo
	enddo

c	do i=1,nvar
c	  print*,(cov(i,j),j=1,nvar)
c	enddo
c	print*,''

c extract correlation matrix: lower triangular part
	var = sigg2 + sigm2 + sigf2 + sigs2 + 1.
	  do i=2,nvar
	    do j=1,i
	     correl(j + ((i-2)*(i-1))/2) = cov(i,j)/var
	    enddo
	  enddo
	return
	end

c computing combinations:
	function combi(n,k)
	combi=1.
	if (k .eq. 0) return
	if (k .eq. n) return
	do i=1,k
	  combi = combi*float(n-i+1)/float(i)
	enddo
	return
	end

  
