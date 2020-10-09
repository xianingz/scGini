#Calculate Gini index for Gamma distribution
gammaGini <- function(alpha){
  nu <- gamma((2*alpha+1)/2)
  de <- alpha*gamma(alpha)*sqrt(pi)
  return(nu/de)
}

#ADAM gradiente descent
est.nb <- function(poi, lr=1e-2){
  alpha=0.1
  ita=0.1
  v.alpha=0
  s.alpha=0
  beta1=0.9
  beta2=0.999
  epsilon = 1e-6
  for(i in c(1:10000)){
    dalpha = sum(digamma(poi+alpha)-digamma(alpha)+log(1-ita))

    v.alpha = beta1*v.alpha + (1-beta1)*dalpha
    s.alpha = beta2*s.alpha + (1-beta2)*dalpha^2
    v.alpha.norm = v.alpha/(1-beta1^i)
    s.alpha.norm = s.alpha/(1-beta2^i)

    alpha = alpha + lr*v.alpha.norm/(sqrt(s.alpha.norm)+epsilon)

    ita = sum(poi)/(length(poi)*alpha+sum(poi))

    if(abs(dalpha) < 1e-5){
      break
    }
  }
  return(c(ita,alpha,i))
}

#Gamma-Possion pdf
gampoi <- function(y, alpha, eta){
  if(y<=100){
    ga1 <- gamma(y+alpha)/(gamma(alpha)*factorial(y))
  }else{
    ga1 <- y^(alpha-1)/(gamma(alpha))
  }
  t2 <- eta^y * (1-eta)^alpha
  return(ga1*t2)
}

#Gama-Poisson cdf
gampoi.cul <- function(y, alpha, est){
  res = 0
  if(y<0){
    return(0)
  }
  for(i in c(0:y)){
    res = res + gampoi(i, alpha, est)
  }
  return(res)
}

##Estimator for each single cell
gnf.nb <- function(data, er=1e-2){
  pvs=c()
  ds=c()
  pb <- txtProgressBar(min = 0, max = dim(data)[2], style = 3)
  alfs <- c()
  ct = 0
  cels = c()
  ginis = c()
  sginis = c()
  cginis = c()
  for(i in colnames(data)){
    #print(i)
    ct = ct + 1
    setTxtProgressBar(pb, ct)
    cont = 0
    sc <- data[,i]
    sc.fre <- as.data.frame(table(sc))
    sc.fre$prob <- sc.fre$Freq/sum(sc.fre$Freq)
    er2 = er
    for(j in c(1:5)){
      sc.est <- tryCatch(
        est.nb(sc, er2),
        error = function(e) {
          NULL
        })
      if(!is.null(sc.est)){
        break
      }
      er2 = er2/(j+1)
    }

    ginis <- c(ginis, ineq::Gini(sc))
    if(is.null(sc.est)){
      alfs <- rbind(alfs, c(NA,NA,NA))
      cginis <- c(cginis, NA)
    }else{
      alfs <- rbind(alfs, sc.est)
      cginis <- c(cginis, gammaGini(sc.est[2]))
    }
  }
  close(pb)
  alfs <- cbind(alfs, ginis, cginis)
  alfs <- data.frame(alfs)
  rownames(alfs) <- colnames(data)
  colnames(alfs) <- c("eta", "alpha", "Iter", "Gini", "cGini")
  return(alfs)
}
