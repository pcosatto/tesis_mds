#Simulacion Main Script
source('aux_functions.R'); source('cmdscaling.R')
n_cores <- parallel::detectCores()
n_cores <- 1

#Primera simulacion
primera <- primera_simulacion(10)

#Simulacion Principal - n fijo
principal_metricas <- metricas_n_fijo(100)


#Simulacion Principal
principal_escenario1 <- simulacion_principal(100,1)
principal_escenario2 <- simulacion_principal(100,2)

#borrar todo lo que no sirve
remove(dim, dims, metodos, n_cores,
         Nrep, scenario, scenarios, size, sizes)
rm(cmdscaling_test,eigenvalues_calculation,simulacion_principal,
   preparar_tabla, eig, B, X, rep, t0, tf, tiempo,
   metricas_bloque2)

#guardar environment
rm(list=setdiff(ls(),setdiff(ls(),lsf.str())))
save.image("Simulacion Main Results.RData")






