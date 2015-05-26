cc data from 2 sisters:
c
c	parameter(npar=4,nn=15,m1=9,m2=5)
c	integer idata(m1,m2)
c	real*8 param(npar),correl(nn*(nn-1)/2)
c	real*8 loglik,prob(m1)
c	data param/-1.2,0.25,0.25,.0/
c	  open(1,file='c:/papers/twins/sister2.dat')
c	  do i=1,m1
c	    read(1,*)(idata(i,j),j=1,m2)
c	  enddo
c
c	  nparent= 0
c	  n1 = 4
c	  n2 = 2
c	  call covsis2(param,npar,nparent,n1,n2,correl)
c	  print*,(correl(i),i=1,10)
cc
c	  call llsis2(param,npar,idata,m1,m2,loglik,prob)
c	  print*,'log likelihood = ',loglik
c	end
c
c
c idata has m2=5 columns: 
c n1=       number of pregnancies from first sister
c npe1=     number among n1 with pe
c n2=       number of pregnancies from second sister, zero is allowed
c npe2=     number among n1 with pe
c nfam =    numer of families with associated characteristics
c
 	subroutine llsis2(param,npar,idata,m1,m2,loglik,prob)
	parameter (nn=15)
	integer idata(m1,m2),m1,m2,npar
	INTEGER ND, INFIN(nn), MAXPTS, INFORM, IVLS
	DOUBLE PRECISION CORREL(nn*(nn-1)/2), LOWER(nn), UPPER(nn), 
     &    DELTA(nn), RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
	double precision loglik,param(npar),prob(m1)

c default values for mvnorm call
	maxpts = 25000
        abseps = 0.0001 
	releps = 0.
c
c parameters: 1=intercept, 2=genetic-maternal,
c 	      3=genetic father, 4=common sibling
	tlim = param(1)
	sigm2 = param(2)
	sigf2 = param(3)
	sigs2 = param(4)
c	tlim = param(1)/sqrt(sigm2 + sigf2 + sigs2 + 1.)
c
	loglik = 0.
c
c process data: nparent = 0 for no parents
c
	do idat=1,m1
          nparent=0
	  n1  = idata(idat,1)
	  n2  = idata(idat,3)
          nvar = nparent + n1 + n2
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covsis2(param,npar,nparent,n1,n2,correl)
c
	  npe1 = idata(idat,2)
	  npe2 = idata(idat,4)
	    do k=1+nparent,nparent+npe1
   	       infin(k) = 0
	    enddo
	    do k=1+nparent+n1,nparent+n1+npe2
		infin(k)=0
	    enddo

	    call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = combi(n1,npe1)* combi(n2,npe2)
	    prob(idat) = cmult * value
c	    print*,(infin(ii),ii=1,nvar),prob(idat)
	    loglik = loglik + idata(idat,5) * log(prob(idat))
c	print*,'nvar=',nvar,(infin(i),i=1,3),'prob=',prob(idat)	  

c end idat
	enddo
	
	loglik = -loglik
	return
	end 

c computing covariance matrices and extracting correlation
c
c n1 = number of daughters from first father
c n2 = number of daughters from second father
c
	subroutine covsis2(param,npar,nparent,n1,n2,correl)
	parameter (nn=15)
	real*8 param(npar),cov(nn,nn),correl(*)
c workspace
	real Ag(nn,nn),Am(nn,nn),Af(nn,nn),As(nn,nn)

c parameters: 1=intercept, 2=genetic-maternal,
c 	      3=genetic-paternal, 4=common sibling
	tlim = param(1)
	sigm2 = param(2)
	sigf2 = param(3)
	sigs2 = param(4)
	nvar = nparent + n1 + n2

c setup covariation matrix: 
	  do i=1,nvar
	    do j=1,nvar
              Am(i,j) = 1.
	      Af(i,j) = 1.
	      As(i,j) = 1.
	  enddo
	  enddo
c
c genetic-maternal and paternal
c
	  do i=nparent+n1+1,nvar
	    do j=1,nparent+n1
		Am(i,j)=0.5
		Af(j,i)=0.
	    enddo
	  enddo
	  do i=1,nparent+n1
	    do j=nparent+n1+1,nvar
		Am(i,j)=0.5
		Af(j,i)=0.
	    enddo
	  enddo

c total variance
	do i =1,nvar
	  do j=1,nvar
	  cov(i,j)=sigm2*Am(i,j) + sigf2*Af(i,j) + sigs2*As(i,j)
	  enddo
	enddo

c	do i=1,nvar
c	  print*,(cov(i,j),j=1,nvar)
c	enddo
c	print*,''

c extract correlation matrix: lower triangular part
	var = sigm2 + sigf2 + sigs2 + 1.
	  do i=2,nvar
	    do j=1,i
	     correl(j + ((i-2)*(i-1))/2) = cov(i,j)/var
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
