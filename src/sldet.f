c from det-work.f in papers/ascertainment/survival/fortran
c
c subroutine to compute the penalty = sum(log(det(..)))
c
	subroutine sldet(nr,nc,W,invD,vare,logdet)
c        parameter (LDA=2, pi=3.141593d0)
c	parameter (LDA=nr, pi=3.141593d0)
	parameter (pi=3.141593d0)
	double precision W(nr,nc), invD(nr,nr), vare

        double precision A(nr,nr), logdet
        CHARACTER*1  UPLO

	LDA = nr

	logdet = 0.d0
	do ic = 1,nc
c for every col of W:
c create nrxnr matrix A = (diag(W) + vare*invD)/(2*pi*vare)
	  do i=1,nr
	    do j=1,nr
               A(i,j) = invD(i,j)/(2.d0*pi)
            enddo
            A(i,i) = (W(i,ic) + vare*invD(i,i))/(2.d0*pi*vare)
           enddo
c          det  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
c         print*,A,det
        UPLO = 'U'
        call DPOTF2( UPLO, nr, A, LDA, INFO)
	
	do i=1,nr
          logdet = logdet + 2.d0*dlog(A(i,i))
        enddo
c end ic loop
        enddo

	return
	end

