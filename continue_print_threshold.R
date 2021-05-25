  ## this is a modified function of the GRINDR continue function, that not only plots the bifurcation diagram but also
  #returns the threshold value.
  continue_print_threshold <- function(state=s, parms=p, odes=model, x=1, step=0.01, xmin=0, xmax=1, y=2, ymin=0, ymax=1.1, log="", time=0, positive=FALSE, add=FALSE, ...) 
  {  
    # continue a steady state
    dots <- list(...)
    if (!is.null(dots)) {
      unknown <- names(dots[!names(dots) %in% c(args_steady,args_plot)])
      if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
      dots_steady <- dots[names(dots) %in% args_steady]
    }else dots_steady <- NULL
    if (!is.numeric(x)) x <- index(x,names(parms))
    if (!is.numeric(y)) y <- index(y,names(state))
    logx <- ifelse(grepl('x',log), TRUE, FALSE)
    clrs <- c("red","black","blue")
    lwds <- c(2,1,1)
    if (missing(xmax) & parms[x] >= 1) xmax <- 2*parms[x]
    if (missing(xmin) & parms[x] < 0) xmin <- 2*parms[x]
    if (!missing(xmin) & xmin >= parms[x]) stop("xmin should be smaller than parameter")
    if (!missing(xmax) & xmax <= parms[x]) stop("xmax should be larger than parameter")
    
    FUN <- function(lastState,lastDom,step) {
      lastP <- p0
      preLastState <- lastState
      nok <- 0
      while (xmin < lastP & lastP < xmax & ymin < lastState[y] & lastState[y] < ymax) {
        if (logx) parms[x] <- lastP*(1+step)
        else parms[x] <- lastP + step
        #predState <- lastState + (lastState-preLastState)
        q <- do.call('steady',c(list(y=lastState,fun=odes,parms=parms,time=time,positive=positive),dots_steady))
        newState <- q$y  # should be steady state and closeby
        if (attr(q,"steady") & abs(sum(newState-lastState))/(1e-9+abs(sum(lastState))) < 0.1) {
          jac <- jacobian.full(y=newState,fun=odes,parms=parms)
          dom <- sign(max(sort(Re(eigen(jac)$values))))
          if (dom != lastDom) cat("Bifurcation at",names(parms[x]),"=",parms[x],"\n")
          if (logx) lines(c(parms[x]/(1+step),parms[x]),c(lastState[y],newState[y]), col=clrs[dom+2],lwd=lwds[dom+2])
          else lines(c(parms[x]-step,parms[x]),c(lastState[y],newState[y]), col=clrs[dom+2],lwd=lwds[dom+2])
          preLastState <- lastState
          lastState <- newState
          lastDom <- dom
          lastP <- parms[x]
         # lalala <- NA
          if (nok > 10 & abs(step) < actualStep) step <- sign(step)*min(2*abs(step),actualStep)
          nok <- nok + 1
        }else{
          nok <- 0
          if (abs(step) > actualStep/100) step <- step/2
          else{ # Go back one step, overpredict, call steady, and turn
            parms[x] <- lastP
            predState <- lastState + 5*(lastState-preLastState)
            q <- do.call('steady',c(list(y=predState,fun=odes,parms=parms,time=time,positive=positive),dots_steady))
            newState <- q$y  # should be steady state and not the same
            #print(c(lastState,predState,newState,parms[x]))
            if (attr(q,"steady") & abs(sum(newState-lastState))/(1e-9+abs(sum(lastState))) > 0.001) {
              cat("Turning point point at",names(parms[x]),"=",parms[x],"\n")
              lalala<-parms[x]
              print("lalala")
              jac <- jacobian.full(y=newState,fun=odes,parms=parms)
              dom <- sign(max(sort(Re(eigen(jac)$values))))
              middle <- (lastState[y]+newState[y])/2
              lines(c(parms[x],parms[x]),c(lastState[y],middle), col=clrs[lastDom+2],lwd=lwds[lastDom+2])
              lines(c(parms[x],parms[x]),c(middle,newState[y]), col=clrs[dom+2],lwd=lwds[dom+2])
              step <- -step
              preLastState <- lastState
              lastState <- newState
              lastDom <- dom
              lastP <- parms[x]
              return(lalala)
            }else{
              cat("Final point at",names(parms[x]),"=",parms[x],"\n")
              #lalala <- NA 
              cat("If this looks wrong try changing the step size\n")
              break
            }
          }
        }
      }
     return(lalala)}
    p0 <- parms[x]
    q0 <- do.call('steady',c(list(y=state,fun=odes,parms=parms,time=time,positive=positive),dots_steady))
    
    if (attr(q0,"steady")) {
      cat("Starting at",names(parms[x]),"=",parms[x],"with:\n")
      print(q0$y)
      bary <- q0$y[y]
      if (missing(ymax) & bary >= 1.1) ymax <- 2*bary
      if (missing(ymin) & bary < 0) ymin <- 2*bary
      if (!missing(ymin) & ymin >= bary) stop("ymin should be smaller than y-variable")
      if (!missing(ymax) & ymax <= bary) stop("ymax should be larger than y-variable")
      if (!add)
        do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=names(p0),ylab=names(state)[y],log=log,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
      orgWarn <- getOption("warn")
      options(warn = -1)
      jac <- jacobian.full(y=q0$y,fun=odes,parms=parms)
      dom <- sign(max(sort(Re(eigen(jac)$values))))
      if (logx) actualStep <- step
      else actualStep <- step*xmax
      lalala<-FUN(lastState=q0$y,lastDom=dom,actualStep)
      lalala<-FUN(lastState=q0$y,lastDom=dom,-actualStep)
      options(warn = orgWarn)
    } else cat("No convergence: start closer to a steady state")
  return(lalala)
    }
    