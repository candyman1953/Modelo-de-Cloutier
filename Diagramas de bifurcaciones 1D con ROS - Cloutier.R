setwd("~/Documentos/TESIS/Github 1")
source("Grind.R")
        
#Definimos nuestro sistema de ecuaciones (Cloutier 2012)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dROS <- K1*(1+S1+d_aSYN*((aSYN/k_aSYN)^4)/(1+(aSYN/k_aSYN)^4))-K2*ROS*S2; 
    daSYN <- K3*ROS*S3-K4*aSYN*S4
    return(list(c(dROS, daSYN)))  
  }) 
}  
        
#Declaramos los parámetos
p <- c(S1=0.5, S2=1, S3=1, S4=1, K1=720, K2=720, K3=0.007, K4=0.007, d_aSYN=15, k_aSYN=8.5) # p is a named vector of parameters
#Condiciones iniciales
s <- c(ROS=0.1, aSYN=0.1)
        

#Debemos dar un valor de S1 en el rango de biestabilidad
p["S1"] <- .1
#Calculamos raíz cercana al punto (0,0) con Newton-Raphson y su estabilidad
lowSS <- newton(c(ROS=0,aSYN=0))
#Calculamos raíz cercana al punto (7,7) con Newton-Raphson y su estabilidad
midSS <- newton(c(ROS=7,aSYN=7))
#Calculamos raíz cercana al punto (20,20) con Newton-Raphson y su estabilidad
highSS <- newton(c(ROS=20,aSYN=20))

#Graficamos el diagrama de bifurcaciones en S1 con la variable ROS
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
        
        
#Hacemos el análisis ahora variando a "S2"

p["S1"] <- .1
p["S2"] <- 1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)
        
#Declaramos S2=0.8
p["S2"] <- 0.8
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)

#Declaramos S2=0.6
p["S1"] <- 0.1
p["S3"] <- 1
p["S2"] <- 0.6

lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))

continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=35, log="", time=0, positive=TRUE, add=TRUE)

#############################################################################################
####################################### VARIAMOS A S3 #######################################
#############################################################################################
        
#Damos a S1, S2 y S3 sus valores nominales 
p["S1"] <- 0.1
p["S2"] <- 1
p["S3"] <- 1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
        
#Diagrama con S3=1.3
p["S3"] <- 1.3
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.0001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
        
#Diagrama con S3=1.3
p["S3"] <- 1.8
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
        

############################################################################################
########################## VARIAMOS A S2 Y S3 SIMULTÁNEAMENTE #####################################
############################################################################################

#Regresamos a valores nominales nuevamente
p["S2"] <- 1
p["S3"] <- 1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
        
#Damos S2=0.8 y S3=1.3
p["S2"] <- 0.8
p["S3"] <- 1.3
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
        
#Damos S2=0.6 y S3=1.8 -> Observamos en este caso un sólo estado estable con altos niveles de ROS
p["S2"] <- 0.6
p["S3"] <- 1.8
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, add=TRUE) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=30, log="", time=0, positive=TRUE, add=TRUE)
  