library(CASdatasets) # Requiere instalación de Github
library(Rglpk)# Linear Programming
library(quadprog) # Quadratic Programming
library(Rsolnp) # e General Non-Linear Optimization
library(DEoptim) #Global Optimization by Differential Evolution
library(robustbase) # Basic Robust Statistics
library(readr)
library(readxl)
require(xts)
require(YieldCurve)
library(PerformanceAnalytics)
library(quantmod)
library(tidyverse)
##### Teoría de optimización ###################

## Programación no Lineal ## 
args(Rglpk_solve_LP)


LP_solver <- function(c, cstr = list(), trace = FALSE) {
       Aeq <- Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
       aeq <- Reduce(c, cstr[names(cstr) %in% "aeq"])
       A <- Reduce(rbind, cstr[names(cstr) %in% "A"])
       a <- Reduce(c, cstr[names(cstr) %in% "a"])
       sol <- Rglpk_solve_LP(obj = c,
      mat = rbind(Aeq, A),
      dir = c(rep("==", nrow(Aeq)),
               rep(">=", nrow(A))),
       rhs = c(aeq, a),
       verbose = trace)
           status <- sol$status
        solution <- if (status) rep(NA, length(c)) else sol$solution
        list(solution = solution, status = status)
       }





## Programación quadratica ##
QP_solver <- function(c, Q, cstr = list(), trace = FALSE) {

       Aeq <- Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
      aeq <- Reduce(c, cstr[names(cstr) %in% "aeq"])
      A <- Reduce(rbind, cstr[names(cstr) %in% "A"])
       a <- Reduce(c, cstr[names(cstr) %in% "a"])
       sol <- try(solve.QP(Dmat = Q,
                                dvec = -2 * c,
                                Amat = t(rbind(Aeq, A)),
                                bvec = c(aeq, a),
                                meq = nrow(Aeq)),
                       silent = TRUE)
      if (trace) cat(sol)
      if (inherits(sol, "try-error"))
       list(solution = rep(NA, length(c)), status = 1)
       else
          list(solution = sol$solution, status = 0)
       }

 

## Progamación no Lineal ##
NLP_solver <- function(par, f, cstr = list(), trace = FALSE) {Aeq <- Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
     aeq <- Reduce(c, cstr[names(cstr) %in% "aeq"])
     A <- Reduce(rbind, cstr[names(cstr) %in% "A"])
     a <- Reduce(c, cstr[names(cstr) %in% "a"])
     heq <- Reduce(c, cstr[names(cstr) %in% "heq"])
     h <- Reduce(c, cstr[names(cstr) %in% "h"])
     leqfun <- c(function(par) c(Aeq %*% par), heq)
     eqfun <- function(par)
       unlist(lapply(leqfun, do.call, args = list(par)))
     eqB <- c(aeq, rep(0, length(heq)))
    
       lineqfun <- c(function(par) c(A %*% par), h)
       ineqfun <- function(par)
         unlist(lapply(lineqfun, do.call, args = list(par)))
       ineqLB <- c(a, rep(0, length(h)))
       ineqUB <- rep(Inf, length(ineqLB))
      
         sol <- solnp(par = par,
                        fun = f,
                        eqfun = eqfun,
                        eqB = eqB,
                        ineqfun = ineqfun,
                        ineqLB = ineqLB,
                        ineqUB = ineqUB,
                        control = list(trace = trace))
        
           status <- sol$convergence
           solution <- if (status) rep(NA, length(par)) else sol$pars
           list(solution= solution, status = status)
            }





## Datos ## 
# Tickers
n_derivados <-c( "FB", "GOOG",  "TAOP", "BAC","AMD","GE","QCOM","AAL","UAL","MGM")   

# Descaragar Datos de YAHOO
for(i in  n_derivados) {
  getSymbols(i, src = "yahoo", from = "2010-01-01", to = "2014-07-30", periodicity = "daily")
  
}

# Obtener Valores Ajustados
list <- lapply(n_derivados, function(x) Cl(get(x)))
precio.cierre <- do.call(merge,list)


# Calculo de Retornos
retornos <- data.frame(apply(precio.cierre, 2, function(x) Delt(x, type = "log")),
                       fecha = index(precio.cierre)) %>%  na.omit() 

x<-retornos%>%select(-"fecha")
x <-x %>% rename(FB = FB.Close, GOOG = GOOG.Close, TAOP = TAOP.Close, BAC = BAC.Close,
                 AMD=AMD.Close, , GE = GE.Close, QCOM= QCOM.Close,  AAL= AAL.Close,MGM=MGM.Close, UAL=UAL.Close,  ) 
x<- as.matrix(x)

colnames(x)

## Restricciones ##

# Retorno Esperado
targetReturn <- function(x, target) {
   list(Aeq = rbind(colMeans(x)), aeq = target)
   }


# Full inversión 
fullInvest <- function(x) {
   list(Aeq = matrix(1, nrow = 1, ncol = ncol(x)), aeq = 1)
   }



# Long
longOnly <- function(x) {
   list(A = diag(1, ncol(x)), a = rep(0, ncol(x)))
  }




MV_QP <- function(x, target, Sigma = cov(x), ...,
                  cstr = c(fullInvest(x),
                              targetReturn(x, target),
                              longOnly(x), ...),
                  trace = FALSE) { size <- ncol(x)
                   c <- rep(0, size)
                   Q <- Sigma
                  
                     # optimization
                     sol <- QP_solver(c, Q, cstr, trace)
                    
                       # extract weights
                       weights <- sol$solution
                       names(weights) <- colnames(x)
                       weights
                     }


# Modelo Clásico
w <-MV_QP(x ,  mean(x) )    
barplot(w, ylim = c(0, 1), las = 2,border="red",
        col="blue",xlab = "Acciones", ylab = "Peso de la acción (w)",
        density=10,  main = "Optimización Clásica MV")




# Ahora con estimación robusta de la varianza
w2 <- MV_QP(x, mean(x), Sigma = covMcd(x)$cov)
barplot(w2, ylim = c(0, 1), las = 2,border="red",
        col="blue",xlab = "Acciones", ylab = "Peso de la acción (w)",
        density=10,  main = "Optimización RObusta MV")


# Portafolio de Minima varianza 
w3 <- MV_QP(x, cstr = c(fullInvest(x), longOnly(x)))
barplot(w3, ylim = c(0, 1), las = 2,border="red",
        col="blue",xlab = "Acciones", ylab = "Peso de la acción (w)",
        density=10,  main = "Modelo de Mínima Varianza")




# Conditional Value at Risk
CVaR_LP <- function(x, target, alpha = 0.95, ...,
                     cstr = c(fullInvest(x),
                                targetReturn(x, target),
                                longOnly(x)),
                     trace = FALSE) {# number of scenarios
   J <- nrow(x)
  
     # number of assets
     size = ncol(x)
    
       # objective coefficients
       c_weights <- rep(0, size)
       c_VaR <- 1
       c_Scenarios <- rep(1 / ((1 - alpha) * J), J)
       c <- c(c_weights, c_VaR, c_Scenarios)
      
         # extract values from constraint to extend them
         # with CVaR constraints
         Aeq <- Reduce(rbind, cstr[names(cstr) %in% "Aeq"])
         aeq <- Reduce(c, cstr[names(cstr) %in% "aeq"])
         A <- Reduce(rbind, cstr[names(cstr) %in% "A"])
         a <- Reduce(c, cstr[names(cstr) %in% "a"])
        
           # build first two blocks of the constraint matrix
           M1 <- cbind(Aeq, simple_triplet_zero_matrix(nrow(Aeq), J + 1))
           M2 <- cbind(A, simple_triplet_zero_matrix(nrow(A), J + 1))
          
             # identity matrix and vector of zeros
             I <- simple_triplet_diag_matrix(1, J)
            
               # block CVaR constraint (y x  alpha  z_j >= 0)
               M3 <- cbind(x, rep(1, J), I)
              
                 # block CVaR constraint (z_j >= 0)
                 M4 <- cbind(simple_triplet_zero_matrix(J, size+  1), I)
                
                   # vector of zeros used for the rhs of M3 and M4
                   zeros <- rep(0, J)
                  
                     # combine constraints
                     cstr <- list(Aeq = M1,
                                    aeq = aeq,
                                    A = rbind(M2, M3, M4),
                                    a = c(a, zeros, zeros))
                    
                       # optimization
                       sol <- LP_solver(c, cstr, trace = trace)
                      
                         # extract weights
                         weights <- sol$solution[1:size]
                         names(weights) <- colnames(x)
                        
                           # extract VaR and CVaR
                           VaR <- sol$solution[size + 1]
                           CVaR <- c(c %*% sol$solution)
                           attr(weights, "risk") <- c(VaR = VaR, CVaR = CVaR)
                          
                             weights
                           }


round(CVaR_LP(x, mean(x)), 3)







