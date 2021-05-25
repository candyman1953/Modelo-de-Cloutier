setwd("~/Documentos/TESIS/Github 1")
source("Grind.R")
source("continue_print_threshold.R")

#Definimos nuestro sistema de ecuaciones (Clouthier 2012)
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dROS <- K1*(1+S1+d_aSYN*((aSYN/k_aSYN)^4)/(1+(aSYN/k_aSYN)^4))-K2*ROS*S2; 
    daSYN <- K3*ROS*S3-K4*aSYN*S4
    return(list(c(dROS, daSYN)))  
  }) 
}  

#Declaramos los parámetros
p <- c(S1=0.5, S2=1, S3=1, S4=1, K1=720, K2=720, K3=0.007, K4=0.007, d_aSYN=15, k_aSYN=8.5) # p is a named vector of parameters
#Condiciones iniciales
s <- c(ROS=0.1, aSYN=0.1)

#Graficamos el diagrama de bifurcación con valores nominales
#Este paso es necesario para poder guardar valores más adelante
p["S1"] <- .1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=24) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=24, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=0, xmax=4,y="ROS", ymin=0, ymax=24, log="", time=0, positive=TRUE, add=TRUE)

#Jugamos con los índice para ver qué valores de Si's nos da valores S1 mayores a 0
#A partir de ello, partimos la forma de llenar la matriz

#Partimos los valores de los intervalos a analizar para S2 y S3
#S2 de 0.1 a 1 cada 0.1 unidades
S2_3D=seq(0.1,1,0.1)
#S3 de de 1 a 2.5 cada 0.1 unidades
S3_3D=seq(1,2.5,0.1)

#Definimos nuestra matriz con 16 renglones (S3) y 10 columnas (S2)
S1_crit_values <- matrix(nrow=16, ncol=10, NA)


p["S1"] <- 0
S3_values=seq(1,1.2,0.1)
S2_values=seq(.5,1,0.1)
#Definimos los índices
ii=1
#jj=5 para recorrer los valores de la matriz y comenzar con S2=0.5 
jj=5
for (S2 in S2_values ){
  # a?adimos un for para S3
  for (S3 in S3_values) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    # xmin debe tener un valor menorr a 0
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    #Guardamos el valor de S1 crit
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=1
  jj=jj+1
}

##Una vez siguiendo los comentarios y recomendaciones de dónde hay valores
##bien definidos para S1, pude ver las casillas a llenar

p["S1"] <- 0
#Llenaremos los valores de S3 de 1.3 a 2.1 y de S2 de 0.9 a 1
S3_values2=seq(1.3,2.1,0.1)
S2_values2=seq(0.9,1,0.1)
ii=4
#Para comenzar con S3=1.3
jj=9
#Para comenzar con S2=0.9
for (S2 in S2_values2 ){
  for (S3 in S3_values2) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE) 
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=4
  jj=jj+1
}

#Va de nuevo

p["S1"] <- 0

S3_values3=seq(1.3,1.7,0.1)
S2_values3=seq(0.7,0.8,0.1)


ii=4
#Para S3 desde 1.3
jj=7
#Para S2 desde 0.7

for (S2 in S2_values3 ){
  for (S3 in S3_values3) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=4
  jj=jj+1
}

#Again...

p["S1"] <- 0
S3_values4=seq(2.2,2.3,0.1)
S2_values4=seq(1,1,0.1)

ii=13
#Comenzamos con S3=2.2
jj=10
#S2=1

for (S2 in S2_values4 ){
  for (S3 in S3_values4) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=13
  #Ya s?lo es un valor para S3, por ello ii=13
  jj=jj
}

#De nuevo...

p["S1"] <- 0
S3_values5=seq(1.8,1.9,0.1)
S2_values5=seq(0.8,0.8,0.1)

ii=9
#S3 desde 1.8
jj=8
#S2 = 0.8

for (S2 in S2_values5 ){
  for (S3 in S3_values5) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=9
  jj=jj
}

#?ltimo para valores de S2 y S3 con los que tiene sentido calcular el valor de S1

p["S1"] <- 0
S3_values6=seq(1.3,1.4,0.1)
S2_values6=seq(0.6,0.6,0.1)

ii=4
#TOmamos S3 de 1.3 a 1.4
jj=6
#S2 = 0.6
for (S2 in S2_values6 ){
  for (S3 in S3_values6) {
    p["S2"] = S2 
    p["S3"] = S3 
    S2
    S3
    
    lowSS <- newton(c(ROS=0,aSYN=0))
    midSS <- newton(c(ROS=7,aSYN=7))
    highSS <- newton(c(ROS=100,aSYN=100))
    
    if (!is.null(lowSS)){
      S1Thres=continue_print_threshold(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }else{ S1Thres=continue_print_threshold(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, add=TRUE) 
    }
    continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-1, xmax=5,y="ROS", ymin=-1, ymax=100, log="", time=0, positive=TRUE, add=TRUE)
    
    S1_crit_values[ii,jj]=S1Thres
    S1Thres
    S1_crit_values
    ii=ii+1
  }
  ii=9
  jj=jj
}


##### AQUÍ TERMINAN LOS VALORES DE S1 QUE ESTÁN BIEN DEFINIDOS (S1>0)#####

#Imprimimos la matriz
S1_crit_values

###Ahora le ponemos cero a las entradas vacías (donde aparece "NA")###
for(d in 1:nrow(S1_crit_values)) {
  for(f in 1:ncol(S1_crit_values)) {
    if(is.na(S1_crit_values[d,f] == TRUE)){
      S1_crit_values[d,f] <- 0
    }
  }
}

#Verificamos que todo esté en orden
S1_crit_values

### ?Yeeeeeeiiih! :D ###

#Guardamos la matriz en un documento .csv:
write.csv(S1_crit_values, file="ValoresdeS1.csv")

#Llamamos a la paquetería para poder imprimir con  la matriz
library(plot3D)
library(scatterplot3d)
#Imprimimos la matriz con colores (2D)
image2D(S1_crit_values, x = S3_3D, y = S2_3D,  clab="S1 crit", xlab = "S3", ylab = "S2")

#"Levantamos" la matriz de acuerdo al valor de S1_crit (3D)
persp3D(x = S3_3D, y = S2_3D, z = S1_crit_values, xlab = " ", ylab = " ", zlab = " ", xlim = c(1,2.5), ylim = c(0.1,1), zlim = c(0,2.5), 
        ticktype = "detailed", inttype = 2, main = "inttype = 2")

#Le ponemos título que indica que la gráfica de valores umbrales y cuadrícula para
#mejor visualización
persp3D(x = S3_3D, y = S2_3D, z = S1_crit_values, xlab = " ", ylab = " ", zlab = " ", xlim = c(1,2.5), 
        ylim = c(0.1,1), zlim = c(0,2.5), ticktype = "detailed", bty = "b2", main = "Valores umbrales de S1 variando S2 y S3")