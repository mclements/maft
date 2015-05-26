c 15 May 2008
c assumes Sigma is the same across families, so correl is only a vector
c
	subroutine proba0(tnr,nr,ncomp,np,prob,varcomp,beta,Y,correl,
     &             X,vars,ortho,relerr,abserr)
	parameter (nn=15)
	double precision varcomp(*),vars(*),beta(*),prob(*),correl(*)
        double precision sumvar,relerr,abserr
        integer tnr,nr,ncomp,np, Y(*), X(*), ortho
c local workspace
	integer yi(nn), infin(nn)
	DOUBLE PRECISION  LOWER(nn), UPPER(nn), xb(nn), DELTA(nn), 
     &      RELEPS, ABSEPS, ERROR, VALUE, MVINIT, MVFUNC



c default values for mvnorm call
	maxpts = 25000
	releps = relerr
        abseps = abserr
	do i=1,nn
	  delta(i)=0.
	enddo

	do idat = 1,nr
c  define matrices
	  do i =1,tnr
            yi(i) = Y((idat-1)*tnr +i)
            xb(i) = 0.
            do j = 1,np
               Xij = X((idat-1)*tnr*np + (j-1)*tnr +i)
               xb(i) = xb(i) + Xij*beta(j)
            enddo
          enddo

c
	sumvar=1.
   	do i =1,ncomp
	  sumvar = sumvar + varcomp(i)
	enddo

c standardize
c if (ortho .eq. "yes") then    ## assumed yes
            do i = 1,tnr
              xb(i) = xb(i)*dsqrt(sumvar)/sqrt(vars(i))
            enddo
c      endif

        do i = 1,tnr
           if(yi(i)>0.5) then
              lower(i) = -1000.
              upper(i) = xb(i)
              infin(i) = 0
           endif
           if (yi(i)<0.5) then
              lower(i) = xb(i)
              upper(i) = 1000.
              infin(i) = 1
           endif
        enddo
       call MVTDST(tnr, 0, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM ) 

        prob(idat) = value
c        print*, prob(idat)
c end idat loop
	enddo        

 	end

