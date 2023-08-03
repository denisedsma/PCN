# Script para gerar simulações
#
# - gera matriz baseado nas probabilidades de ising
# - inserir criterio para escolha de numero de iteraçoes
# - criar matriz com borda fixa, como na analise do pantanal
# - Novas simulaçoes para 2a ordem completa e 2a ordem incompleta
# - Escolha aleatoria do pixel avaliado, nao linear


#### Rodar funcoes para gerar amostras ####

setwd("C:/Users/debor/Documents/Mestrado/Dissertacao/pós dissertação/script")
source("pcn_algorithm.R")  


#### Gerar matrizes de ordem completa ####

#funcao auxiliar
rotate <- function(x) apply(x, 2, rev)

#Define a estrutura da arvore de dependencia
#Qual a ordem de dependencia
ord <- 2
#Quantos quadradinhos em cada ordem
nviz<-NULL
for(o in 1:ord) {
  config<-2*(2*o+1)+2*(2*o-1)
  nviz<-c(nviz,config)
}
#Valores possiveis de s_i - soma de pretos (+1) e brancos (-1) - dado a ordem
valor1 <- seq(from = (sum(nviz[1:ord])), to = -(sum(nviz[1:ord])), by=-2)
valor1 #de todos pretos ate 0 pretos

valor <- (sum(nviz[1:ord]):0) #numero de pretos na vizinhanca 
valor# de todos pretos ate 0 pretos

#Atencao: a ordem de "valor1" e "valor" deve ser a mesma!!! de max numero de pretos ate o min.

#vetor de probabilidade (de ising). prob(site i ser preto|vizinhanca)
beta <- 0.05
prob <- 1/(1+exp(-2*beta*(valor1))) #formula de ising = 1/(1+exp(-2*beta*s))



# Simula matriz de ordem COMPLETA!
# Gera lattice atraves da funcao genlattice
# n: lado da matriz 
# valor: possiveis valores de s_i (utilizar vetor gerado anteriormente)
# prob: vetor de probabilidades correspondente para cada valor de vizinhanca <valor>
# ord: ordem de dependencia
# burnin: numero de itera??es
# ini: ini = 0, se gostaria de comecar as itera?oes a partir de uma configura??o aleatoria;
#      ini = <amostra>, se gostaria de come?ar a partir de uma amostra j? existente;
# v1 e v2 : configuracao (numero de brancos) que sera monitarada entre iteracoes 
# iter : a funcao para de rodar quando a dif da proporcao de v1 e v2 na amostra fica abaixo de 1e-3 por <iter> iteracoes consecutivas
# jump : plota matriz e iter1 e iter2 a cada <jump> iteracoes
genlattice <- function(valor, prob, ord, n, v1 = 4, v2 = 8, iter = 10,
                       ini = 0, jump = 10, D = trunc(log(n/2))) {
  
  if (length(ini)==1){ #cria matriz aleatoria como amostra inicial
    m2 <- matrix(1*(0.5>runif((n+2*D)*(n+2*D))), n+2*D) #1 para valores menores que 0.5 e 0 para maiores que 0.5
    m2 <- m2 +1  #o que era 1 fica 2(preto), o que era 0 vira 1(branco)
    x <- m2[(D+1):(D+n),(D+1):(D+n)]
  }else{ #utiliza matriz inputada
    m2 = ini
    x = m2[(D+1):(D+n),(D+1):(D+n)]
  }
  #imagem da matriz inicial
  image(1:(n+2*D), 1:(n+2*D), t(rotate(m2)), col=gray(1:0), xlab="", ylab="")
  title(main = paste0("k =", 0))
  
  #freq da amostra original
  f0 <- arv.amostral(m2, n, 2)
  freq <- rbind(c(v1, ( f0[[1]][which(f0[[1]][,1]==v1),3] )/n^2, NA),
                c(v2, ( f0[[1]][which(f0[[1]][,1]==v2),3] )/n^2, NA) )
  iter1 <- 0
  iter2 <- 0
  
  #variavel para criar estimativa intervelar (salva 100 matrizes depois de convergir) 
  # *** convergiu <- FALSE 
  # *** c <- 0
  
  for (k in 1:1e3) { #loop das iteracoes
    for ( j in sample(1:n, n, replace = FALSE)) { #loop das colunas da matriz
      for ( i in sample(1:n, n, replace = FALSE)) { #loop das linhas da matriz
         
            s1 <- sum(m2[(D+i-ord):(D+i+ord),(D+j-ord):(D+j+ord)]-1) -1*(m2[(D+i),(D+j)]==2) #quantos pretos em D^j
            p1 <- prob[which(valor == s1)] 
            x[i,j] <- 1*(p1>runif(1)) ## como no algoritmo de Metropolis Hastings. Se p<runif, entao x[i,j] = 0 se p>runif, entao x[i,j]=1
            x[i,j] <- x[i,j] + 1 #transforma em escala de 2 e 1
            m2[(D+i),(D+j)] = x[i,j]

      }
    }
    #contagem de certas vizinhancas na amostra
    f1<-0
    f2<-0
    for (j in 1:n){
      for (i in 1:n){
        s <- ordem(m2, i, j, 1, 2, n) #numero de brancos na primeira ordem
        if(s == v1) f1 <- f1 + 1
        if(s == v2) f2 <- f2 + 1
      }
    }
    #tabela com proporcao das vizinhancas por iteracao
    freq <- rbind(freq, c(v1, f1/n^2, abs((f1/n^2)-freq[2*k-1,2])),
                  c(v2, f2/n^2, abs((f2/n^2)-freq[2*k,2])) )
    
    #numero de vezes q a diferenca entre proporcoes consecutivas é menor que 1e-3
    if (freq[2*k+1,3]<3e-3) iter1 <- iter1 + 1
    else iter1 <- 0
    if (freq[2*k+2,3]<3e-3) iter2 <- iter2 + 1
    else iter2 <- 0
    
    #plota matriz e iter1 e iter2 a cada <jump> iteracoes
    if (k %% jump == 0) {
      image(1:n, 1:n, t(rotate(x)), col=gray(1:0), xlab="", ylab="")
      title(main = paste0("k =", k))
      print(c(iter1, iter2))
      
      odd <- seq(1,2*k+2,2)
      even <- seq(2,2*k+2,2)
      #track difference between proportion of configurations from one sample to the other
      #ideally, we want to select when the values are bellow the red line, but bellow blue line will be good too.
      plot(freq[odd,3], type = "l", 
           main = paste0("difference of proportion of ",v1,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      plot(freq[even,3], type = "l", 
           main = paste0("difference of proportion of ",v2,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      
    }
    
    #se o numero de iteracoes consecutivas (iter1 e iter2)  for >= ao numero
    #pre-stabelecido <iter>, salva a matriz e termina a funcao.
    if( iter1 >= iter && iter2 >= iter){ # *** && convergiu == FALSE ){
      #plota grafico que acompanha "estabilizacao" da matriz
      image(1:n, 1:n, t(rotate(x)), col=gray(1:0), xlab="", ylab="")
      title(main = paste0("k =", k))
      odd <- seq(1,2*k+2,2)
      even <- seq(2,2*k+2,2)
      #track difference between proportion of configurations from one sample to the other
      #ideally, we want to select when the values are bellow the red line, but bellow blue line will be good too.
      plot(freq[odd,3], type = "l", 
           main = paste0("difference of proportion of ",v1,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      plot(freq[even,3], type = "l", 
           main = paste0("difference of proportion of ",v2,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      
      #print frequency table, and other values
      print(freq)
      print(c(k, iter1, iter2))
      
      saveRDS(freq, paste0("freq_",v1,"e",v2,"m_b005_2acompleta_iter10.rds"))
      
      #convergiu
      # *** convergiu <- TRUE
      break
    } 
    
    # *** if(convergiu == TRUE){
    # ***  write(m2, file = paste0("m",k,"_b005_2completa_iter10.txt"))
    # ***  c <- c + 1
    # ***}
    # ***if(c==100)  break
    
  }
  #amostra final
  m2
}


#### Simulações - 2a ordem completa ####
## beta = 0.05

start <- Sys.time() #19min
m2comp <- genlattice(valor = valor, prob = prob, ord = ord, ini = 0,
                     n = 500, v1 = 4, v2 = 8, iter = 10, jump = 25)

end <- Sys.time() 
end - start


start <- Sys.time() #11min
arv2comp <- arv.amostral(teste, 500, 2) 
end <- Sys.time() 
end - start

ind2comp <- Vdj(arvteste, 500)
poda2comp <- poda(indteste, 2)
poda2comp
#recupera 2 ordem completa (148 contextos)


# Gerar estimativa intervalar #

#gerar 100 matrizes (alterei genlattice para salvar 100 matrizes depois de convergir)
#comentarios da funcao genlattice com ***
start <- Sys.time() #19min
genlattice(valor = valor, prob = prob, ord = ord, ini = 0,
           n = 500, v1 = 4, v2 = 8, iter = 10, jump = 25)

end <- Sys.time() 
end - start

prob2 <- NULL
d <- 0
n <- 500
D <- trunc(log(n/2))

#gerar arvores e aplicar pcn pra amostra
start <-Sys.time()
for(i in 18:117 ){ #14 horas
  m <- scan(paste0("m",i,"_b005_2completa_iter10.txt"))
  m <- matrix(m, (n+2*D) , (n+2*D))
  arv <- arv.amostral(m, n, 2) 
  saveRDS(arv, file = paste0("arv",i,"_2acompleta_b005_iter10.rds"))
  ind <- Vdj(arv, n)
  pod <- poda(ind, 2)
  if( dim(pod[[1]])[1]==0 & dim(pod[[3]])[1]==0 ){
    d <- d + 1
    prob2 <- rbind(prob2, cbind(arv[[2]][,1:4],(arv[[2]][,4]-arv[[2]][,3])/arv[[2]][,4]))
  }
}
end <- Sys.time()
end-start

#quantas matrizes recuperam uma arvore de 2a ordem completa
d 

#probibilidade de site ser preto
#saveRDS(prob2, file = "prob_est_2aordem_completa.rds")


#se as 100 foram recuperadas, esperaria que prob1 teria 600 linhas
#calculate quantiles
est2<-NULL
config<-NULL
for (i in 0:8){
  for(j in 0:16){
    if ( !is.null(nrow(prob2[which(prob2[,1]==i & prob2[,2]==j),])) & nrow(prob2[which(prob2[,1]==i & prob2[,2]==j)!=0,]) ){
      est <- apply(prob2[which(prob2[,1]==i & prob2[,2]==j),c(1,2,5)], 2, quantile, probs=c(0.025,0.5,0.975))
      est2<-rbind(est2, c(i, j, est[7:9])) #prob de ser preto
    } else{
      config<-rbind(config, c(i,j))
    }
  }
}

#vizinhancas que nao apareceram nas 100 matrizes
config 

#write(est2, "IC_2acompleta_b005.txt") #151x5
est2 <- scan("IC_2acompleta_b005.txt")
est2 <- matrix(est2, 151, 5 )

#quantile type=8
est2_type8<-NULL
config<-NULL
for (i in 0:8){
  for(j in 0:16){
    if ( !is.null(nrow(prob2[which(prob2[,1]==i & prob2[,2]==j),])) & nrow(prob2[which(prob2[,1]==i & prob2[,2]==j)!=0,]) ){
      est <- apply(prob2[which(prob2[,1]==i & prob2[,2]==j),c(1,2,5)], 2, quantile, probs=c(0.025,0.5,0.975),type=8)
      est2_type8<-rbind(est2_type8, c(i, j, est[7:9])) #prob de ser preto
    } else{
      config<-rbind(config, c(i,j))
    }
  }
}

#vizinhancas que nao apareceram nas 100 matrizes
config 

#write(est2_type8, "IC_2acompleta_b005_type8.txt") #151x5
est2_type8 <- scan("IC_2acompleta_b005_type8.txt")
est2_type8 <- matrix(est2_type8, 151, 5 )





#### Gerar matrizes de ordem variavel ####

#funcao auxiliar
rotate <- function(x) apply(x, 2, rev)

#Define a estrutura da arvore de dependencia
#Qual a ordem de dependencia
ord <- 2
#Quantos quadradinhos em cada ordem
nviz<-NULL
for(o in 1:ord) {
  config<-2*(2*o+1)+2*(2*o-1)
  nviz<-c(nviz,config)
}
#Valores possiveis de s_i - soma de pretos (+1) e brancos (-1) - dado a ordem
valor1 <- seq(from = (sum(nviz[1:ord])), to = -(sum(nviz[1:ord])), by=-2)
valor1 #de todos pretos ate 0 pretos

#vetor de probabilidade (de ising). prob(site i ser preto|vizinhanca)
beta <- 0.05
prob1 <- 1/(1+exp(-2*beta*(valor1))) #formula de ising = 1/(1+exp(-2*beta*s))
prob1


genmix <- function(valor1 = valor1, prob1 = prob1, ord, n, v1 = 4, v2 = 8, iter = 10,
                   ini = 0, jump = 10, D = trunc(log(n/2))){ 
  
  if (length(ini)==1){ #cria matriz aleatoria como amostra inicial
    m2 <- matrix(1*(0.5>runif((n+2*D)*(n+2*D))), n+2*D) #1 para valores menores que 0.5 e 0 para maiores que 0.5
    m2 <- m2 +1  #o que era 1 fica 2(preto), o que era 0 vira 1(branco)
    x <- m2[(D+1):(D+n),(D+1):(D+n)]
  } else{ #utiliza matriz inputada
    m2 = ini
    x = m2[(D+1):(D+n),(D+1):(D+n)]
  }
  
  #imagem da matriz inicial
  image(1:(n+2*D), 1:(n+2*D), t(rotate(m2)), col=gray(1:0), xlab="", ylab="")
  title(main = paste0("k =", 0))
  
  #freq da amostra original
  f0 <- arv.amostral(m2, n, 2)
  freq <- rbind(c(v1, ( f0[[1]][which(f0[[1]][,1]==v1),3] )/n^2, NA),
                c(v2, ( f0[[1]][which(f0[[1]][,1]==v2),3] )/n^2, NA) )
  iter1 <- 0
  iter2 <- 0
  
  for (k in 1:1e3) {
    for ( j in 1:n) {
      for ( i in 1:n) {
        
        s1 <- sum(m2[(D+i-1):(D+i+1),(D+j-1):(D+j+1)]-1) -1*(m2[(D+i),(D+j)]==2) #quantos pretos na primeira frame
        p1 <- prob1[which(valor1 == (s1-(8-s1)))] 
        x[i,j] <- 1*(p1>runif(1)) ## como no algoritmo de Metropolis Hastings. Se p<runif, entao x[i,j] = 0 se p>runif, entao x[i,j]=1
        x[i,j] <- x[i,j] + 1 #transforma em escala de 2 e 1
        m2[(D+i),(D+j)] = x[i,j] 
        
        if( s1 %in% c(3,4,5)){
          s2 <- sum(m2[(D+i-ord):(D+i+ord),(D+j-ord):(D+j+ord)]-1) -1*(m2[(D+i),(D+j)]==2) #quantos pretos na segunda frame
          p2 <- prob1[which(valor1 == (s2-(24-s2)))] 
          x[i,j] <- 1*(p2>runif(1)) ## como no algoritmo de Metropolis Hastings. Se p<runif, entao x[i,j] = 0 se p>runif, entao x[i,j]=1
          x[i,j] <- x[i,j] + 1 #transforma em escala de 2 e 1
          m2[(D+i),(D+j)] = x[i,j]
        } 
        
      }
    }
    
    #contagem de certas vizinhancas na amostra
    f1<-0
    f2<-0
    for (j in 1:n){
      for (i in 1:n){
        s <- ordem(m2, i, j, 1, 2, n) #numero de brancos na primeira ordem
        if(s == v1) f1 <- f1 + 1
        if(s == v2) f2 <- f2 + 1
      }
    }
    #tabela com proporcao das vizinhancas por iteracao
    freq <- rbind(freq, c(v1, f1/n^2, abs((f1/n^2)-freq[2*k-1,2])),
                  c(v2, f2/n^2, abs((f2/n^2)-freq[2*k,2])) )
    
    #numero de vezes q a diferenca entre proporcoes consecutivas é menor que 3e-3
    if (freq[2*k+1,3]<3e-3) iter1 <- iter1 + 1
    else iter1 <- 0
    if (freq[2*k+2,3]<3e-3) iter2 <- iter2 + 1
    else iter2 <- 0
    
    #plota matriz e iter1 e iter2 a cada <jump> iteracoes
    if (k %% jump == 0) {
      image(1:n, 1:n, t(rotate(x)), col=gray(1:0), xlab="", ylab="")
      title(main = paste0("k =", k))
      print(c(iter1, iter2))
      
      odd <- seq(1,2*k+2,2)
      even <- seq(2,2*k+2,2)
      #track difference between proportion of configurations from one sample to the other
      #ideally, we want to select when the values are bellow the red line, but bellow blue line will be good too.
      plot(freq[odd,3], type = "l", 
           main = paste0("difference of proportion of ",v1,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      plot(freq[even,3], type = "l", 
           main = paste0("difference of proportion of ",v2,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      
    }
    
    #se o numero de iteracoes consecutivas (iter1 e iter2)  for >= ao numero
    #pre-stabelecido <iter>, salva a matriz e termina a funcao.
    if( iter1 >= iter && iter2 >= iter){ # *** && convergiu == FALSE ){
      #plota grafico que acompanha "estabilizacao" da matriz
      image(1:n, 1:n, t(rotate(x)), col=gray(1:0), xlab="", ylab="")
      title(main = paste0("k =", k))
      odd <- seq(1,2*k+2,2)
      even <- seq(2,2*k+2,2)
      #track difference between proportion of configurations from one sample to the other
      #ideally, we want to select when the values are bellow the red line, but bellow blue line will be good too.
      plot(freq[odd,3], type = "l", 
           main = paste0("difference of proportion of ",v1,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      plot(freq[even,3], type = "l", 
           main = paste0("difference of proportion of ",v2,"'s"), ylab = "diff", xlab = "iter")
      abline(h=c(3e-3,1e-3), col=c("blue","red"))
      
      #print frequency table, and other values
      print(freq)
      print(c(k, iter1, iter2))
      
      saveRDS(freq, paste0("freq_",v1,"e",v2,"m_b005_2incompleta_iter10.rds"))
      
      #convergiu
      # *** convergiu <- TRUE
      break
    } 
    
    # *** if(convergiu == TRUE){
    # ***  write(m2, file = paste0("m",k,"_b005_2completa_iter10.txt"))
    # ***  c <- c + 1
    # ***}
    # ***if(c==100)  break
    
  }
  #amostra final
  m2
}


#### Simulações - 2a ordem incompleta ####
## beta = 0.05
# se tem 3, 4 ou 5 pretos na 1a ordem, vai pra 2a ordem 

start <- Sys.time() #17.5min
m2incomp <- genmix(valor1 = valor1, prob1 = prob1, ord = ord, ini = 0,
                     n = 500, v1 = 4, v2 = 8, iter = 10, jump = 10)

end <- Sys.time() 
end - start


start <- Sys.time() #10min
arv2incomp <- arv.amostral(m2incomp, 500, 2) 
end <- Sys.time() 
end - start

ind2incomp <- Vdj(arv2incomp, 500)
poda2incomp <- poda(ind2incomp, 2)
poda2incomp
#recupera 2 ordem completa (148 contextos)

#### Testar novo algoritmo (novo pdj) nas arvores de 2a ordem incompleta da dissertacao ####

prob2 <- NULL
prob1 <- NULL
d <- 0
n <- 150

start <-Sys.time()
for(i in c(202:213, 215:218, 220:229, 231:242,
           244, 246:247, 249, 251:258, 260:261,
           263:271, 274:276, 279:280, 283:287,
           289:297, 299) ){ 
  arv <- readRDS(paste0("arv",i,"_2aordem_incomp.rds"))
  ind <- Vdj(arv, n)
  pod <- poda(ind, 2)
  if( dim(pod[[1]])[1]==6 & dim(pod[[3]])[1]==0 ){
    d <- d +1
    frame1<- arv[[1]][which(arv[[1]][,1]%in% c(0,1,2,6,7,8)),]
    prob1<-rbind(prob1, cbind(frame1[,1:3],(frame1[,3]-frame1[,2])/frame1[,3]))
    frame2<- arv[[2]][which(arv[[2]][,1]%in% c(3,4,5)),]
    prob2<-rbind(prob2, cbind(frame2[,1:4],(frame2[,4]-frame2[,3])/frame2[,4]))
  }
}
end <- Sys.time()
end-start

#Apenas i = 251 recupera a arvore correta, de 81 matrizes!!!
d 
