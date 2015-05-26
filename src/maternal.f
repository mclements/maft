cc general maternal halfsibs, with 2 fathers
c	parameter(npar=5,nn=15,m1=18,m2=6)
c	integer idata(m1,m2)
c	real*8 param(npar),correl(nn*(nn-1)/2)
c	real*8 loglik,prob(m1)
c	data param/-1.,1.,0.25,0.25,.5/
c	  open(1,file='c:/papers/twins/maternal.dat')
c	  do i=1,m1
c	    read(1,*)(idata(i,j),j=1,m2)
c	  enddo
cc
cc	  nparent= 1
cc	  n1 = 2
cc	  n2 = 2
cc	  call covmat(param,npar,nparent,n1,n2,correl)
cc	  print*,(correl(i),i=1,10)
cc
c	  call llmat(param,npar,idata,m1,m2,loglik,prob)
c	  print*,'log likelihood = ',loglik
c	end
c
c
c idata has m2=6 columns: 
c mpe = mother pe status 0,1,99
c n1=       number of daughters from first father
c npe1=     number among n1 with pe
c n2=       number of daughters from second father, zero is allowed
c npe2=     number among n1 with pe
c nfam =    numer of families with associated characteristics
c
 	subroutine llmat(param,npar,idata,m1,m2,loglik,prob)
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
c parameters: 1=intercept, 2=genetic-fetal, 3=genetic-maternal,
c 	      4=common family, 5=common sibling
	beta0 = param(1)
	sigg2 = param(2)
	sigm2 = param(3)
	sigf2 = param(4)
	sigs2 = param(5)
	tlim = beta0/sqrt(sigg2 + sigm2 + sigf2 + sigs2 + 1.)
c
c reparameterize:
	tlim = param(1)	
	loglik = 0.
c
c process data
c
	do idat=1,m1
	  nparent=1
	  if (idata(idat,1) .eq. 99) nparent = 0
	  n1  = idata(idat,2)
	  n2  = idata(idat,4)
          nvar = nparent + n1 + n2
	  do j=1,nvar
	    infin(j) = 1
	    lower(j) = tlim
	    upper(j) = tlim
	  enddo
	  call covmat(param,npar,nparent,n1,n2,correl)
c
 	  mpe = idata(idat,1)
	  npe1 = idata(idat,3)
	  npe2 = idata(idat,5)
	    if (mpe .eq. 1) infin(1) = 0
	    do k=1+nparent,nparent+npe1
   	       infin(k) = 0
	    enddo
	    do k=1+nparent+n1,nparent+n1+npe2
		infin(k)=0
	    enddo
c	    print*,(infin(ii),ii=1,nvar)
	    call MVTDST(nvar, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 
	    cmult = combi(n1,npe1)* combi(n2,npe2)
	    prob(idat) = cmult * value
	    loglik = loglik + idata(idat,6) * log(prob(idat))
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
	subroutine covmat(param,npar,nparent,n1,n2,correl)
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
	nvar = nparent + n1 + n2

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
c genetic-fetal
c
	  do i=nparent+1,nparent+n1
	    do j=nparent+n1+1,nparent+n1+n2
		Ag(i,j)=0.25
		Ag(j,i)=0.25
	    enddo
	  enddo

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
