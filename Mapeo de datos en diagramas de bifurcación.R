library("ggplot2")
setwd("~/Documentos/TESIS/Datos gráficas 3")

#Datos de niveles de expresión de a-sinucleína publicados en #Virdi et al
aSYN_expr <- read.csv(file="aSYn_Prod_datos.csv")
View(aSYN_expr)
#Graficamso los datos de expresión
barplot(aSYN_expr[,2], main = "Niveles de expresión del gen SNCA", names.arg = (aSYN_expr[,1]), 
        col = c("#49B6FF", "#A480CF", "#FF499E"),
        ylab="", ylim = c(0,5), las=1)

#Datos de niveles de ROS celular publicados en #Virdi et al
ROS_data <- read.csv(file="ROS_data.csv")
View(ROS_data)
#Graficamos los datos de ROS
barplot(ROS_data[,2], main = "Porcentaje de incremento ROS al día 28", names.arg = (ROS_data[,1]), 
        col = c("#49B6FF", "#A480CF", "#FF499E"),
        ylab="Porcentaje de incremento ROS", ylim = c(0,200), las=1)

#Niveles de GSH en neuronas dopaminérficas publicados en Virdi et al.
GSH_data <- read.csv(file="GSH_data.csv")
#Normalizamos los datos
#WT
WT_GSH=GSH_data[1,2]/100
WT_GSH
#A53T
A53T_GSH <- GSH_data[2,2]/100
A53T_GSH  #0.788
#Tripl
Tripl_GSH <- GSH_data[3,2]/100
Tripl_GSH #0.915

View(GSH_data)
#Graficamos los valores de GSH
barplot(GSH_data[,2], main = "Niveles de GSH al día 28", names.arg = (GSH_data[,1]), 
        col = c("#49B6FF", "#A480CF", "#FF499E"),
        ylab="Niveles de GSH", ylim = c(0,100), las=1)


############## VISUALIZACIÓN DE DATOS EN DIAGRAMAS DE BIFURCACIÓN #########################
###########################   LLAMAMOS A GRIND ################################################
source("Grind.R")
library(plot3D)
library(scatterplot3d)

#Definismo el sistema de Cloutier et al, 2012
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    # NO DEJAR ESPACIOS
    dROS <- K1*(1+S1+d_aSYN*((aSYN/k_aSYN)^4)/(1+(aSYN/k_aSYN)^4))-K2*ROS*S2; 
    daSYN <- K3*ROS*S3-K4*aSYN*S4
    #
    return(list(c(dROS, daSYN)))  
  }) 
}  

#Valor de los parámetos
p <- c(S1=0, S2=1, S3=1, S4=1, K1=720, K2=720, K3=0.007, K4=0.007, d_aSYN=15, k_aSYN=8.5) # p is a named vector of parameters
#Condiciones iniciales
s <- c(ROS=0.1, aSYN=0.1) 

############################################################################################
##################  MAPEO DE DATOS DE EXPRESIÓN DE A-SINUCLEÍNA  ##########################
############################################################################################
#Diagrama de bifurcación en S3 con valores nominales
p["S3"] <- 1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S3", step=0.001, xmin=0.99, xmax=5,y="ROS", ymin=0, ymax=20) 
continue(state=midSS, parms=p, odes=model, x="S3", step=0.001, xmin=0.99, xmax=5,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S3", step=0.001, xmin=0.99, xmax=5,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)

#Guardamos como ROSss al valor estacionario de ROS en S3=1 (valores nominales)
ROSss <- 1.002907
#Graficamos línea vertical en S1=0
abline(v=1, lty=2, col="#49B6FF")
#Fraficamos el valor de ROSss en y con x=0=S1
scatter2D(x=1, y=ROSss, add = TRUE, col = c("#49B6FF"))

#Calculamos los valores estacionarios de ROS de los mutantes con regla de tres:
#Y graficamos los niveles de expresión
#A53T
A53T_ROS <- ROSss*(ROS_data[2,2]/100)
A53T_ROS #1.731017
abline(v=aSYN_expr[2,2], lty=2, col="#A480CF")
scatter2D(x=aSYN_expr[2,2], y=A53T_ROS , col = c("#A480CF"), add = TRUE)
aSYN_expr[2,2]
#Error_GSH_A53T <- A53T_ROS -1.278782
#Error = 0.4522355

#Niveles de ROS y expresión en variante de la triplicación del gen SNCA 
Tripl_ROS <- ROSss*(ROS_data[3,2]/100)
Tripl_ROS
abline(v=aSYN_expr[3,2], lty=2, col="#FF499E")
scatter2D(x=aSYN_expr[3,2], y=Tripl_ROS , col = c("#FF499E"), add = TRUE)
v=aSYN_expr[3,2]
#Error_Tripl_A53T <- Tripl_ROS - 1.09745

############### Promediamos los decaímientos de funciones mitocondriales de Zambón #########
DM_A53T <- (0.35*3 + 0.3)/4
DM_A53T  #0.3375  
DM_Tripl <- (0.9+0.4)/4
DM_Tripl  # 0.0325

##############################################################################################
#######################  DIAGRAMAS DE BIFURCACIÓN MAPEANDO DATOS  #############################
#############################################################################################

#Diagrama de pacientes A53T
p["S1"] <- 0
p["S2"] <- A53T_GSH
p["S4"] <- 0.65
p["S3"] <- 1
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1, y="ROS", ymin=0, ymax=22) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1,y="ROS", ymin=0, ymax=22, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1,y="ROS", ymin=0, ymax=22, log="", time=0, positive=TRUE, add=TRUE)
abline(v=DM_A53T, lty=2, col="#A480CF")
scatter2D(x=DM_A53T, y=A53T_ROS , col = c("#A480CF"), add = TRUE)

#Diagrama de pacientes con triplicación del gen SNCA
p["S1"] <- 0
p["S2"] <-Tripl_GSH
p["S4"] <- 1
p["S3"] <- 1/0.7
lowSS <- newton(c(ROS=0,aSYN=0))
midSS <- newton(c(ROS=7,aSYN=7))
highSS <- newton(c(ROS=20,aSYN=20))
continue(state=lowSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1.5, y="ROS", ymin=0, ymax=20) 
continue(state=midSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1.5,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
continue(state=highSS, parms=p, odes=model, x="S1", step=0.001, xmin=-0.1, xmax=1.5,y="ROS", ymin=0, ymax=20, log="", time=0, positive=TRUE, add=TRUE)
#S1crit = 0.8616562  
abline(v=DM_Tripl, lty=2, col="#FF499E")
scatter2D(x=DM_Tripl, y=Tripl_ROS , col = c("#FF499E"), add = TRUE)

#### Umbral: S1= 2.025208 
###Error 0.102091
1.38759 - A53T_ROS
  ###0.3434275
