"mctest" <- function(x, y, paired = TRUE, nsim = 1000, plot = TRUE) {
  if(missing(x))
    x <- eval(parse(prompt="enter a vector with data values \n(an object, an expression or use the format c() to enter a numerical vector),\n  x = "))
  if(missing(y))
    y <- eval(parse(prompt="enter a vector with data values \n(an object name, an expression or use the format c() to enter a numerical vector),\n  y = "))
  nx <- length(x)
  ny <- length(y)
  if(paired){
    if(nx != ny)
      stop("this function requires x and y with same length")
    df <- nx - 1
    dxy <- x-y
    d <- mean(dxy)/sqrt(var(dxy)/nx)
    for (i in 1:nsim) {
      r <- runif(nx)
      xx <- x*(r<0.5)+y*(r>=0.5)
      yy <- x*(r>=0.5)+y*(r<0.5)
      dxy <- xx - yy
      d <- c(d, mean(dxy)/sqrt(var(dxy)/nx))
    }
  }
  else{
    xy <- c(x, y)
    df <- nx + ny - 2
    s2 <- ((nx-1) * var(x) + (ny-1) * var(y))/df
    d <- (mean(x)-mean(y))/sqrt(s2)
    for (i in 1:nsim) {
      ind <- sample(1:(nx+ny))
      xx <- xy[ind <= nx]
      yy <- xy[ind > nx]
      s2 <- ((nx-1) * var(xx) + (ny-1) * var(yy))/df
      d <- c(d, (mean(xx)-mean(yy))/sqrt(s2*((1/nx)+(1/ny))))
    }
  }
  Pmc <- sum(d>=d[1])/(nsim+1)
  ##  p2 <- sum(d<=d[1])/(nsim+1)
  Pt <- pt(d[1], df = df, lower=FALSE)
  res <- list(p=c(upper.tail = Pmc, lower.tail = 1-Pmc),
              pt = c(upper.tail = Pt, lower.tail = 1-Pt),
              data.statistic = d[1], sim.statistic = d[-1])
  attr(res, "paired") <- paired
  attr(res, "df") <- df
  class(res) <- "mctest"
  if(plot) plot.mctest(res)
  return(res)
}

"plot.mctest" <- function(x, tcurve = TRUE, ...){
  df <- eval(attr(x, "df"))
  xyhist <- hist(x$sim, prob=TRUE, plot=FALSE)
  ymax <- max(max(xyhist$dens, na.rm=TRUE), dt(0, df = df)) 
  xmax <- max(max(xyhist$breaks, na.rm=TRUE), x$data)
  xmin <- min(min(xyhist$breaks, na.rm=TRUE), x$data)
  ldots <- list(...)
  MCstatistics <- x$sim
  if(is.null(ldots$xlim) & is.null(ldots$ylim))
    hist(MCstatistics, prob=TRUE, xlim = c(xmin, xmax), ylim=c(0, ymax),...)
  else
    hist(x$sim, freq=FALSE, ...)
  arrows(x$data, ymax/5, x$data, 0, lwd=2)
  text(x$data, ymax/5, round(x$data, dig=4), pos=3, col="red")
  st <- round(x$p[1], dig=4)
  text(.75*xmax, .9*ymax,
       substitute(P[upper] == Pval, list(Pval=st)), col="red")
  if(tcurve){
    ap <- function(x) {return(dt(x, df=df))}
    curve(ap(x), from = xmin, to = xmax, add=TRUE)
    prob.t <- pt(x$data.st, df = df, lower=FALSE)
    stt <- round(x$p[1], dig=4)
    text(.75*xmax, .75*ymax,
         substitute(P(t)[upper] == Pval, list(Pval=stt)), col="red")
  }
  return(invisible())
}

"print.mctest" <- function(x, ...){
  paired <- eval(attr(x, "paired"))
  if(paired) cat("Paired data\n")
  else cat("Two-sample data\n")
  cat(paste("data statistics = ", x$data,"\n"))
  cat("\n")
  cat("probabilities based on Monte Carlo simulations:\n")
  print(round(x$p, dig=6))
  cat("\n")
  cat("probabilities based on the \"t\" distribution:\n")
  print(round(x$pt, dig=6))
  return(invisible())
}

