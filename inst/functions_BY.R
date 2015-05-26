#Original code by Ben Yip
#
#Modified by Robert Szulkin, 2015-05-15
#
#Description: Remove all unnecessary functions and code regarding twins 


## DESRCIPTIONS OF ALL FUNCITONS
##
## THE FUNCTIONS CAN BE DIVIDED IN THREE CATEGORY
## 1) LIKELIHOOD & Std err CALCULATION
## 2) DATA MANAGEMENT/Ascertainment
## 3) SIMULATION


## 1) LIKELIHOOD & Std err CALCULATION
##    fixreg                 obtain the likelihood of survival data, and ignore the familial correlation
##    Q.likelihod.v2         Using Gauss-Seidel optimization procedure to calculate the adjusted profile h-likelihood
##    |
##    -> Stratify            Stratify the data, see (2)
##    -> UP.beta.like        Given var comp, estimate/update fix parameter, beta, and random effects
##    |  |
##    |  -> UPDATE.ystar     
##    |  -> UP.u
##    |     |
##    |     -> Covar.ls.3    Create covariance/correlation matrices for arbritary family structure
##    -> HP.like.nuclear     Dependending the argument of "ret", and will return either the h-likelihood or adjusted h-likelihood 
##       |
##       -> L1fun            Obtain the "non-random" part of h-likelihood 
##       -> L2fun            Obtain the random part of h-likelihoo/adj h-likelihood
##          |
##          -> sldet         obtain the log(detD) part by calling fortran dode
##    sandwich               sandwich variance formula for variance components
##    |   |
##    |   -> L1fun            
##    |   -> L2fun            
##    |
##    -> deriv1              obtain partial derivatives of fn(x)  
##    Vbeta.nuclear          obtain std error of the fixed parameters
##    OBTAIN.wi              obtain the w-elements, see HA paper for more details 

# 
# YP:
# 24 June 2008
#    give natural names to all variables 
#
# 16 June 2008 version 
# removing the surv.strat list structure    
# allows 
#    - twins or nuclear family
#    - mixture of family types
# ascer.prob = 1 or 'p' (name of ascer.prob in the survdat)
Gauseidel <- function(comp,betas,model,fdat,stratum,subjects='famid',famtype='nuclear',
                      entry,event,ascerp,uloop=7,bloop=30,iter){
  if(missing(ascerp)) return('arg ascerp is missing: ascerp=1 no ascertainment, otherwise a vector of ascerp prob')
  if(missing(event)) return('must specify event')
  
  starting <- c(comp,betas)
  iter.est <- matrix(0,nrow=length(starting)+5,ncol=iter+1)
  iter.est[,1] <- c(starting,rep(0,5))
  rownames(iter.est) <- c(names(comp),names(betas),c('u1','u2','u3','hlike','hplike'))
  betas.up <- betas
  U.up <- 0
  for(i in 1:iter){
    res <- Q.likelihood.v2(comp,model=model,betas=betas.up,
                           fdat=fdat,stratum=stratum,OPT=FALSE,Uvect=U.up,
                           subjects=subjects,famtype=famtype,entry=entry,event=event,
                           ascerp=ascerp,limit=c(0,10),uloop=uloop,bloop=bloop)
    betas.up <- res$Fix
    U.up <- res$U
    fdatb <- cbind(fdat,U=res$U,MU=c(res$MU), YC=c(res$YC), BC=c(res$BC))
    fdatb <- Stratify(fdatb,stratum)
    res2 <- try(optim(comp,HP.like.nuclear,xdat=fdatb,subjects=subjects,ret='adj.hlike',
                 entry=entry, famtype=famtype, event=event, ascer.prob=ascerp,hessian=TRUE),TRUE)
    if(!inherits(res2,'try-error')){
      iter.est[,i] <- c(res2$par,res$Fix,res$U[c(1,11,23)],res$hlike,res2$val)
      comp <- res2$par
    }
    if(inherits(res2,'try-error')) iter <- i
    print(i)
  }
  return(list(fix.obj=res,var.obj=res2,summary=iter.est))
}






UP.beta.nu <- function(model,survdat,startval,ascer.prob=1,
                       family.type=c('nuclear','twin'),
                       entry='entry', event='death',
                       subject='famid',
                       varcomp,
                       rel.err=0.01,max.iter=10,
                       SE=FALSE){
  
  ## define the fixed, random effects  
  famid <- unique(survdat[,subject])
  random.effect=c('g','c','a','f','A','C')
  tf <- terms(model,specials='r')
  fix.effect <- rownames(attr(terms(model), "factors"))[-attr(tf, "specials")$r]
  Y <- as.character(as.formula(model)[2])
  fix.effect <- fix.effect[!fix.effect == Y]
  nrPar <- length(fix.effect)
  tf <- unlist(strsplit(rownames(attr(tf, "factors"))[attr(tf, "specials")$r],''))
  r.effect <- random.effect[random.effect %in% tf]

    
  ## some fix constant
  fixed  <- startval[fix.effect]
  rand.var <- varcomp[r.effect]
  var.e <- varcomp['e']

    if(sum(colnames(survdat)%in%c('tnr','knr2'))!=2){
      return('tnr=total nr of member, knr2=total nr of child indicator is missing for nuclear data') 
    }
    if(sum(colnames(survdat)%in%c('tnr','knr2'))==2){
      famtype <- survdat[,c('tnr','knr2')]
      famtype <- famtype[,1]*10 + famtype[,2] #apply(famtype,1,paste,collapse=' ')
      uftype = unique(famtype)
    }

  
  ## iterate the estimation back and forward
  ## dar <- 1:nrow(survdat)
  fixed.up <- fixed*0
  iter <- 0
  U = rep(0, nrow(survdat))

  while(!(sum((abs(fixed.up-fixed)/fixed) <rel.err)==length(fix.effect) || iter==max.iter)){
    if(iter!=0){fixed <- fixed.up}
    
    for(itype in uftype){
      pick = famtype==itype
      tempdat <- survdat[pick,]
      XB <- tempdat[,fix.effect]%*%fixed
      MU <- XB+ U[pick]
      Ystar <- UPDATE.ystar(tempdat[,Y],tempdat[,event],MU,var.e,tempdat[,entry])
      yc <- Ystar-XB
      
       
        nrchild <- tempdat[1,'knr2'];Tnr <- tempdat[1,'tnr']; type='fam.nuc'
      ui <- UP.u(rand.var,var.e,r.effect,Tnr,nrchild,type,yc)     
      U[pick] <- ui
    }
    MU <- survdat[,fix.effect]%*%fixed +  U
    Ystar <- UPDATE.ystar(survdat[,Y],survdat[,event],MU,var.e,survdat[,entry])
    YC <- Ystar-U

# YP: simpler lines
    ifelse(ascer.prob==1, weight <- ascer.prob, weight <- 1/survdat[,ascer.prob])
    XYC <- crossprod(survdat[,fix.effect],weight*YC)
    XX <- crossprod(survdat[,fix.effect]*sqrt(weight))
    ## update beta
    fixed.up <- solve(XX,XYC)
    iter <- iter+1; ##cat(iter, fixed.up, var(U), '\n')
  }
  
  mu <- survdat[,fix.effect]%*%fixed.up + U
  YC <- survdat[,Y] - mu
  ifelse(!missing(entry),  BC <- survdat[,entry] - mu, BC <- -Inf )## so that pnorm(BC)~0
  return(list(U=U, YC=YC, BC=BC, 
              MU=mu, event=survdat[,event],
              famtype=famtype,
              est.beta=fixed.up,iteration=iter))
}


OBTAIN.wi <- function(y,d,mu,b,var.e){
  sd.e <- sqrt(var.e)
  m <- (y-mu)/sd.e
  V <- dnorm(m)/(1-pnorm(m))
  if(sum(b)!=0){ 
    m.star <- (b-mu)/sd.e
    V.star <- dnorm(m.star)/(1-pnorm(m.star))
  }
  if(sum(b)==0){V.star <- 0; m.star <- 0}
  V[V %in% c(NaN,Inf,-Inf)] <- 0
  V.star[V.star %in% c(NaN,Inf,-Inf)] <- 0
  w <- d + (1-d)*(V*(V-m)) - V.star*(V.star-m.star)
  return(w)
}



Q.likelihood.v2 <- function(comp,model,betas,fdat,subjects,famtype,ascerp=1,
                            entry,event,stratum=c('tnr','knr2'),Uvect=0,
                            limit,OPT=TRUE,uloop=5,bloop=30,rel.err=0.01,abs.err=0.001,
                            out.iter=10){
  
  if(any(comp<0.000001)) return(exp(20))
  if(sum(comp)>limit[2] | sum(comp)<limit[1]) return(exp(20))
  Hlike <- exp(20)
  iter <- 0
  step2 <- 1
  while( (abs(step2 - Hlike)) > 3 || iter<out.iter){
    if(iter==0){betas.up <- betas
                if(length(Uvect)==1){U <- rep(0,nrow(fdat))}
                if(length(Uvect)!=1){U <- Uvect}}
    if(iter>0){U <- step1$U; Hlike <- step2}
    step1 <- UP.beta.like(model,survdat=fdat,startval=betas.up,ascer.prob=ascerp,varcomp=comp,U=U,
                          family.type=famtype,entry=entry,event=event,subject=subjects,
                          rel.err=rel.err,bloop=bloop,abs.err=abs.err,uloop=uloop)
    betas.up <- c(step1$est)
    names(betas.up) <- names(betas)
    fstrat <- cbind(fdat,U=c(step1$U),MU=c(step1$MU),YC=c(step1$YC), BC=c(step1$BC))
    fstrat <- Stratify(fstrat, strata=stratum)
    step2 <- HP.like.nuclear(comp,xdat=fstrat,subject=subjects,famtype=famtype,
                             entry=entry,event=event,ascer.prob=ascerp,
                             ret='hlike')

    iter <- iter+1
  }
  if(Hlike >= exp(20)) return(Hlike)
  hplike <- HP.like.nuclear(comp,xdat=fstrat,subject=subjects,famtype=famtype,
                            entry=entry,event=event,ascer.prob=ascerp,
                            ret='adj.hlike')
  if(OPT==TRUE){  return(hplike) }
  if(OPT==FALSE){
    list(U             = step1$U,
         MU            = step1$MU,
         YC            = step1$YC,
         BC            = step1$BC,
         hlike         = step2,
         Adj.prof.like = hplike,
         Varcomop      = comp,
         iter          = iter,
         Fixest        = betas.up)
  }
}

UP.beta.like <- function(model,survdat,startval,ascer.prob=1,
                         family.type=c('nuclear','twin'),
                         entry,event,subject,varcomp,U=0,
                         rel.err=0.01,bloop=50,abs.err=0.001,uloop=10){

  ## define the fixed, random effects  
  famid <- unique(survdat[,subject])
  random.effect=c('g','a','A','c','f')
  tf <- terms(model,specials='r')
  fix.effect <- rownames(attr(terms(model), "factors"))[-attr(tf, "specials")$r]
  Y <- as.character(as.formula(model)[2])
  fix.effect <- fix.effect[!fix.effect == Y]
  nrPar <- length(fix.effect)
  tf <- unlist(strsplit(rownames(attr(tf, "factors"))[attr(tf, "specials")$r],''))
  r.effect <- random.effect[random.effect %in% tf]
    
  ## some fix constant
  fixed  <- startval[fix.effect]
  rand.var <- varcomp[r.effect]
  var.e <- varcomp['e']
  survdat <- cbind(survdat,U=U)


  ## define the data type
  
  if(family.type=='nuclear'){
    if(sum(colnames(survdat)%in%c('tnr','knr2'))!=2){
      return('tnr=total nr of member, knr2=total nr of child indicator is missing for nuclear data') 
    }
    if(sum(colnames(survdat)%in%c('tnr','knr2'))==2){
       famtype <- survdat[,c('tnr','knr2')]
       famtype <- famtype[,1]*10 + famtype[,2] 
       uftype = unique(famtype)
    }
  }
  
  ## iterate the estimation back and forward
  ## dar <- 1:nrow(survdat)
  fixed.up <- fixed*0
  beta.iter <- 0
  while(!(sum(abs(fixed.up-fixed)<rel.err)==length(fix.effect) || beta.iter==bloop)){
    if(beta.iter!=0){fixed <- fixed.up}    

    for(itype in uftype){
      tempdat <- survdat[famtype==itype,]
      XB <- tempdat[,fix.effect]%*%fixed
      ui <- tempdat[,'U']
      ui.up <- ui+0.1
      
       
        nrchild <- tempdat[1,'knr2'];Tnr <- tempdat[1,'tnr']; type='fam.nuc'
      u.iter <- 0
      while(!(sum((abs(ui-ui.up))<abs.err)==length(ui) | u.iter==uloop)){
        if(u.iter>0){ui <- ui.up}
        MU <- XB+ui
        Ystar <- UPDATE.ystar(tempdat[,Y],tempdat[,event],MU,var.e,tempdat[,entry])
        yc <- Ystar-XB
        ui.up <- UP.u(rand.var,var.e,r.effect,Tnr,nrchild,type,yc)
        u.iter <- u.iter+1
      }
      survdat[famtype==itype,'U'] <- ui.up
    }
    U <- survdat[,'U']
    MU <- survdat[,fix.effect]%*%fixed +  U
    Ystar <- UPDATE.ystar(survdat[,Y],survdat[,event],MU,var.e,survdat[,entry])
    YC <- Ystar-U
    
    ifelse(ascer.prob==1, weight <- ascer.prob, weight <- 1/survdat[,ascer.prob])
    XYC <- crossprod(survdat[,fix.effect],weight*YC)
    XX <- crossprod(survdat[,fix.effect]*sqrt(weight))
    ## update beta
    fixed.up <- solve(XX,XYC)
    beta.iter <- beta.iter+1; ## cat(iter, fixed.up, '\n')
  }
  U <- survdat[,'U']
  MU <- survdat[,fix.effect]%*%fixed +  U
  YC <- survdat[,'Y']-MU
  BC <- survdat[,entry]-MU
  return(list(U=U,MU=MU,YC=YC,BC=BC,beta.iter=beta.iter,u.iter=u.iter,est.beta=fixed.up))
}  

UPDATE.ystar <- function(y,d,mu,var.e,b){
  sd.e <- sqrt(var.e)
  m <- (y-mu)/sd.e
  V <- dnorm(m)/(1-pnorm(m))
  if(sum(b)!=0){ 
    m.star <- (b-mu)/sd.e
    V.star <- dnorm(m.star)/(1-pnorm(m.star))
  }
  if(sum(b)==0){V.star <- 0}
  V[V %in% c(NaN,Inf,-Inf)] <- 0
  V.star[V.star %in% c(NaN,Inf,-Inf)] <- 0
  ystar <- y*d + (1-d)*(mu + sd.e*V)-sd.e*V.star
  return(ystar)
}

UP.u <- function(rand.var,var.e,r.effect,Tnr,nrchild,ftype,yc){
  
  
    Sigma <- Covar.ls.3(rand.var,r.effect,tnr=Tnr,knr2=nrchild,fam.type=ftype)
    A <- solve(diag(Tnr) + var.e*solve(Sigma))
    yc <- matrix(yc,nr=Tnr)
    ## u <-  crossprod(A,yc) #apply(yc,2,FUN=function(x){solve(diag(Tnr) + var.e*solve(Sigma), x) })
    u <-  A%*%yc # no need to to crossprod
  
  return(c(u))
}

Covar.ls.3 <-   function(comp,effect=c('g','f','c'),tnr,knr2=2,fam.type='fam.nuc',corr.mat='n'){
  
  
    r.g  <- matrix(0.5,ncol=tnr,nrow=tnr)
    pnr <- tnr-knr2
    if(pnr>1){
      r.g[1:pnr,1:pnr] <- 0 # assume that the data is sorted with (parent,child)
    }
    diag(r.g) <- 1
    z.g <- diag(tnr)

    ## r.f <- 1
    ## z.f <- matrix(1,ncol=1,nrow=tnr)
    r.f <- matrix(1,tnr,tnr)
    z.f <- diag(tnr)

##     r.a <- 1
##     z.a <- matrix(c(rep(1,pnr),rep(0,knr2)),ncol=1,nrow=tnr)
    r.a1 <- matrix(1,pnr,pnr)
    r.a4 <- matrix(0,knr2,knr2)
    r.a2 <- matrix(0,nc=knr2,nr=pnr)
    r.a3 <- matrix(0,nr=knr2,nc=pnr)
    r.a <- rbind(cbind(r.a1,r.a2),
                 cbind(r.a3,r.a4))
    z.a <- diag(tnr)

    r.c1 <- matrix(0,pnr,pnr)
    r.c4 <- matrix(1,knr2,knr2)
    r.c2 <- matrix(0,nc=knr2,nr=pnr)
    r.c3 <- matrix(0,nr=knr2,nc=pnr)
    r.c <- rbind(cbind(r.c1,r.c2),
                 cbind(r.c3,r.c4))
    z.c <- diag(tnr)

    r.A <- r.a
    diag(r.A) <- 1
    z.A <- diag(tnr)

    r.C <- r.c
    diag(r.C) <- 1
    z.C <- diag(tnr)

  
    ## check <- !effect %in% c('g','a','A','c','C','u1','u2','f','PC')
    check <- !effect %in% c('g','a','c','A','C','f')
    if(sum(check)>0){      return(paste(effect[check],'not allow'))     }
        
    R <- NULL
    Z <- NULL
    R$g <- r.g
    R$f <- r.f
    R$a <- r.a
    R$A <- r.A
    R$c <- r.c
    R$C <- r.C

    Z$g <- z.g
    Z$a <- z.a
    Z$c <- z.c
    Z$A <- z.A
    Z$C <- z.C
    Z$f <- z.f
  

  names(comp) <- effect 
  r <- R[effect]
  z <- Z[effect]
  if(corr.mat %in% c('y','Y','yes','Yes')){    return( list(R=r,Z=z) )    }

  nr <- length(comp)
  sigma <- matrix(0,tnr,tnr)
  for(i in 1:nr){
    sigma <- sigma+r[[i]]*comp[i]
  }  
  return( sigma )
}






HP.like.nuclear <- function(varcomp,xdat,subjects,entry,event,
                            famtype,fix.obj,ascer.prob=1,ret=c('adj.hlike'))
{
  ## sum the adj (pseudo) prof likelihood 
  if(any(varcomp<0.000001)) return(exp(20))
  if(ret=='hlike'){
    a <- 0
    for(i in 1:length(xdat)){
      idat = xdat[[i]]
      prob <- ifelse(ascer.prob==1, 1, idat[1,ascer.prob])
      L1 <- L1fun(varcomp=varcomp,xdat=xdat[[i]],event=event,error='e') 
      L2 <- L2fun(varcomp=varcomp,xdat=idat,subject=subjects,
                  entry=entry,event=event,ret='hlike')
      a <- a+ (L1 + L2)/prob
    }
    hlike <- a
    return(-hlike)
  }
  if(ret=='adj.hlike'){
    a <- 0
    for(i in 1:length(xdat)){
      idat = xdat[[i]]
      prob <- ifelse(ascer.prob==1, 1, idat[1,ascer.prob])
      L1 <- L1fun(varcomp=varcomp,xdat=xdat[[i]],event=event,error='e') 
      L2PEN <- L2fun(varcomp=varcomp,xdat=idat,subject=subjects,entry=entry,
                     event=event,ret='adj.hlike')
      a <- a+ (L1 + L2PEN)/prob
    }
    adj.hlike <- a
    return(-adj.hlike)
  }
}

 
L1fun <- function(varcomp, xdat, event='death', error='e', individual=FALSE){
  if(sum(c('YC','BC') %in% colnames(xdat))!=2) {
    return('YC or BC are missing in xdat')
  } 
  ## individual=TRUE gives the individual family contribution to log-likelihood 
  var.e <- varcomp[error]
  YC <- xdat[,'YC']
  BC <- xdat[,'BC']
  d <- xdat[,event]
  tnr = xdat[1,'tnr']
  
  if(!individual){
    LIKE <- sum(d)*log(2*pi*var.e)+crossprod(d*YC)/var.e
    standardize <- pnorm(YC/sqrt(var.e))
    ##  standardize[standardize==1] <- 0.999999999
    CENS <- sum((1-d)*log(1-standardize))
    standardize <- pnorm(BC/sqrt(var.e))
    ##  standardize[standardize==1] <- 0.999999999
    TRUN <- sum(log(1-standardize))
    l1 <- -0.5*LIKE + CENS - TRUN
    return(l1)
  }
  if(individual){
    LIKE <- d*(log(2*pi*var.e)+YC^2/var.e)
    CENS <- (1-d)*log(1-pnorm(YC/sqrt(var.e)))
    TRUN <- log(1-pnorm(BC/sqrt(var.e)))
    l1.ind <- -0.5*LIKE + CENS - TRUN
    l1.ind <- colSums(matrix(l1.ind, nrow=tnr))   # combine tnr rows per family
    return(l1.ind)
  }
}


L2fun <- function(varcomp, xdat, subject='famid',
            entry='entry', event='death', ret='adj.hlike',individual=FALSE){

  var.e <- varcomp['e']
  r.effect <- names(varcomp)[names(varcomp)!='e']
  rand.var <- varcomp[r.effect]
  
  Y <- xdat[,1]
  d <- xdat[,event]
  b <- xdat[,entry]
  MU <- xdat[,'MU']
  tnr <- xdat[1,'tnr']
  Knr2 <- xdat[1,'knr2']
  
  COV <- Covar.ls.3(rand.var,r.effect,tnr=tnr,Knr2,'fam.nuc')
  invD <- solve(COV)
  NUCu <- matrix(xdat[,'U'],nr=tnr)
  invD.NUCu <- crossprod(invD, NUCu)   
  

  if (!individual){
    l2nuc <- sum(NUCu * invD.NUCu)
    logdet <- (nrow(xdat)/tnr)*log(det(2*pi*COV))
    L2 <- -0.5*(l2nuc+logdet)
    if(ret=='hlike'){return(L2)}
    W.element <- OBTAIN.wi(Y,d,MU,b,var.e)
    nuc.W <- matrix(W.element,nr=tnr)
    PEN <-  sldet(nuc.W,invD,var.e)           

    L2PEN <- (L2-0.5*PEN)
  }
  
  if (individual){
    l2nuc.ind <- colSums(NUCu * invD.NUCu)
    logdet <- log(det(2*pi*COV))
    L2.ind <- -0.5*(l2nuc.ind + logdet)
    if(ret=='hlike'){return(L2.ind)}
    W.element <- OBTAIN.wi(Y,d,MU,b,var.e)
    nuc.W <- matrix(W.element,nr=tnr)
    PEN.ind <-  sldet(nuc.W, invD, var.e, individual=TRUE)           
    
    L2PEN <- (L2.ind -0.5*PEN.ind)
  }
  return(L2PEN)
}



sldet <- function(W,invD,var.e,individual=FALSE){
  nr = nrow(W)
  nc = ncol(W)
  logdet = 0
  ilogdet = rep(1,nc)
  ret = .Fortran('sldet', 
           nr = as.integer(nr), nc= as.integer(nc),
           W = as.double(W),
           invD = as.double(invD),
           vare = as.double(var.e),
           logdet = as.double(logdet),
           ilogdet = as.double(ilogdet))

  if(!individual) return(ret$logdet)
  if (individual) return(ret$ilogdet)
}





Stratify <- function(fdat,strata){ ## BY 24/6
  famtype <- fdat[,strata]
  #if(sum(strata%in%c('tnr','knr2'))!=2) return('strata need tnr and knr2')
  if(length(strata)==1){
    strat <- famtype
    k=1
  }
  
  if(length(strata)==2){
    if(sum(c('tnr','knr2')%in%strata)==2){
      strat <- famtype[,1]*10 + famtype[,2]
      k=2
    }
    if(sum(c('tnr','knr2')%in%strata)!=2){
      strat <- apply(fdat[,strata],1,paste,collapse=' ')
      k=3	
    }
  }
  
  if(length(strata)==3){ 
    if(sum(c('tnr','knr2','p')%in%strata)==3){
      strat <- famtype[,1]*10 + famtype[,2] + famtype[,3]/10
      k=2
    }
    if(sum(c('tnr','knr2','p')%in%strata)!=3){
      strat <- apply(fdat[,strata],1,paste,collapse=' ')
      k=3
    }
  }
  
  if(length(strata)>3){
    strat <- apply(fdat[,strata],1,paste,collapse=' ')
    k=3
  }
  
  
  uftype <- sort(unique(strat))
  if(k==1) namn <- as.character(uftype)
  if(k==2) namn <- apply(sapply(strsplit(as.character(uftype),''),function(x){paste(strata,'=',x,sep='')}),2,paste,collapse=' ')
  if(k==3) namn <- apply(sapply(strsplit(uftype,' '),function(x){paste(strata,'=',x,sep='')}),2,paste,collapse=' ')
  dat.strat <- NULL
  for(i in 1:length(uftype)){
    dat.strat[[i]] <- fdat[strat == uftype[i],]
  }
  names(dat.strat) <- namn 
  return(dat.strat)
}







