#  ALGORITMO PCN 
#
# - Utiliza matrizes que contem borda
# - consertando a formula para N = n^2
# - Ordem maxima da arvore vira argumento da função (facil testar ordens diferentes)

#### ordem ####
# Calcula o valor observado da vizinhaca de ordem o para o site ij
# com argumento novo mirror para escolher se quer ou nao espelhar matriz.
# m2 - matriz (de 2(preto) e 1(branco))
# i - linha
# j - coluna
# o - ordem
# a - alfabeto
# n - tamanho da matriz n x n
#ordem adaptada para lidar com agua e quando a=2!!
ordem <- function(m2, i, j, o, a, n, D = trunc(log(n/2))){
  l1<-m2[(D+i-o),(D+j-o):(D+j+o-1)]
  l2<-m2[(D+i+o),(D+j-o+1):(D+j+o)]
  l3<-m2[(D+i-o+1):(D+i+o),(D+j-o)]
  l4<-m2[(D+i-o):(D+i+o-1),(D+j+o)]
  tbl<-table(c(l1,l2,l3,l4,1:a))
  obs<-matrix(tbl)[,1]-1
if(0 %in% names(tbl)) obs <- obs[1] + 1 + obs[a]
else obs<- obs[-a]
obs
}


#### obs.ordem ####
## Modifica a matriz dados de acordo com a vizinhanha observada (atualiza numero de sites e vizinhancas encontradas)
## Atencao para o uso com alfabeto  pois ordenamos em 1,2,3....max(alfabeto)
obs.ordem<-function(m, dados, achou, obs, o, i, j, a, n, D = trunc(log(n/2))){
  # Se ja existe esse contexto
  if( length(achou)==1){
    dados[[o]][achou,(o*(a-1)+a)]<-dados[[o]][achou,(o*(a-1)+a)] + 1
      if(m[(i+D),(j+D)]<a){
        dados[[o]][achou,(o*(a-1)+m[(i+D),(j+D)])] <- dados[[o]][achou,(o*(a-1)+m[(i+D),(j+D)])]+1
    }
  }
  # Nao existe esse contexto
  if(length(achou)==0){
    dados[[o]]<-rbind(dados[[o]], c(obs[1:(o*(a-1))],rep(0,(a-1)),1))
    if(m[(i+D),(j+D)]<a){
      dados[[o]][dim(dados[[o]])[1],(o*(a-1)+m[(i+D),(j+D)])] <- 1
    }
  }
  dados
}


#### arv.amostral ####
#Cria arvore amostral a partir de uma determinada matriz
#adaptado para agua quando a=2
#D = ordem maxima da arvore
arv.amostral<-function(m, n, a, D = trunc(log(n/2))){
  ## Entre com uma lista que sera usada de apoio
  ## A variavel dados eh uma lista. ## Cada item da lista correspondo a uma determinada ordem de vizinhanca
  ## Se A=3 : Para ordem de vizinhanca 2(porexemplo) dados[[2]] sera uma matrix com 9 colunas pois temos 3 cores
  ## As 6 primeiras sao as vizinhancas de primeira, segunda ordem ## A setima e oitava correpondem ao
  ## numero de brancos e cinzas encontrados
  ## A nona quantas vezes observamos essa vizinhanca independente se ij eh branco cinza ou preto.
  ## Logo se a ordem eh o e o alfabeto(nuemros de cores)=a temos para a*o+(a)
  ##
  dados<-NULL
  ## Essa lista tera log(n/2) elementos= ordem maxima de vizinhanca
  dados<-vector("list", (D))
  obs<-NULL
  # Para encontrar os contexto:
  # Calcule para o site i=j=1 o valor de todas as vizinhas e aloque na lista
  for ( o in 1: D) {
    obs<-c(obs, ordem(m, 1, 1, o, a, n))
    dados[[o]]<-rbind(dados[[o]], c(obs[1:(o*(a-1))],rep(0,a)))
  }
  #Para todos os outros sites faca o mesmo
  for(i in 1:n){
    for ( j in 1: n){
      obs<-NULL
     #Veja se o pixel avaliado é zero (água), se for, pule para o próximo.
        if(m[(i+D),(j+D)]!= 0){
            #Calcule o valor da vizinhanca para cada ordem
            for ( o in 1: D){
            obs<-c(obs, ordem(m, i, j, o, a, n))
            # Tente encontrar se ja existe essa vizinhanca
            # Liste a quantidade de strings diferentes de ordem o (elemento o da lista)
            achou<-1: dim(dados[[o]])[1]
            #achou
            #A cada ordem fique somente com aqueles possuem a mesma vizinhanca
            # Para que encontre um contexto deve ser testado desda raiz
            for( p in 1:(o*(a-1))) achou<-achou[which(dados[[o]][achou,p]==obs[p])]
            #achou
            #Atualize os dados com o contexto novo ou ja existente
            dados<-obs.ordem(m, dados, achou, obs, o, i, j, a, n)
            }
          obs
          achou
        }
      }
    }
  #remove from final tree if total sites observed are 0!
  if(m[(D+1),(D+1)]==0){
    for (ordem in D:1){
      if(dados[[ordem]][1,(ordem+2)]==0){
        dados[[ordem]]<- dados[[ordem]][-1,]
      }
    }
  }
  dados
}



#### pdj ####
#Calcula Pdj
#as vezes (nsa[a]/ns)^nsa[a] vai pra zero e o log disso = -Inf
#fazer nsa[a]*log(nsa[a]/ns) evita esse tipo de coisa
#exemplo: log((844/2822)^844) = -inf, mas 844*log(844/2822) = -1018.749
pdj<-function(dados,o,s,n){
  A = 2 #alfabeto = 2
  aprod<-0
  
  ns<-dados[[o]][s,(o+2)]
  nsa<- c(dados[[o]][s,(o+1)], dados[[o]][s,(o+2)]-dados[[o]][s,(o+1)])
  
  if (any(nsa == 0)){
    for (a in 1:A) {
      aprod<- aprod + log((nsa[a]/ns)^nsa[a])
    }
  }else {
    for (a in 1:A) {
      aprod<- aprod + nsa[a]*log(nsa[a]/ns)
    }
  }
  aprod<- aprod - log(n^(A-1)) #n é o lado da matriz!!
}


#### vdj ####
#Calcula as indicadoras Vdj
Vdj<-function(dados, n){
  ### Calculo dos Vs e Xs ####
  # Acrescente duas colunas
  for ( o in 1: length(dados)) dados[[o]]<-cbind(dados[[o]],0,0)
  # Faca primeiro para ultimos filhos
  u<-length(dados)
  nfilhos<-dim(dados[[u]])[1]
  for( f in 1:nfilhos) dados[[u]][f,(o+3):(o+4) ]<-c(pdj(dados,u,f,n),0)
  # Faca para  o restante
  for( o in (u-1):1){
    # Quantos strings?
    nstring<-dim(dados[[o]])[1]
    for( s in 1:nstring){
      #s=1
      # Quais sao os filhos de s
      filhos<-1: dim(dados[[(o+1)]])[1]
      for (p in 1: o) filhos<- filhos[which(dados[[(o+1)]][filhos,p]==dados[[o]][s,p])]
      #Qual eh o produto dos v dos filhos
      #vprod<-prod(dados[[(o+1)]][filhos,(o+4)])
      vprod<-sum(dados[[(o+1)]][filhos,(o+4)])
      # Testa para encontrar V e X
      dados[[o]][s,(o+3):(o+4)]<-c(max(pdj(dados,o,s,n),vprod),sum(vprod > pdj(dados,o,s,n))) #se true 1, se false 0
    }
  }
  dados
}

#### poda ####
#Podar arvore amostral
poda<-function(arv.dados,a){
  # Para todas as ordens
  for ( o in 1:(length(arv.dados)-1)){
    i=1   # auxiliar
    #TRansformar tudo em matriz
    if( is.vector(arv.dados[[o]]))arv.dados[[o]]<-matrix(arv.dados[[o]],nrow=1)
    if(dim(arv.dados[[o]])[1]==0)i<-dim(arv.dados[[o]])[1]+1 ## caso teime em ser vetor
    # Enquanto tiver vizinhanca nao analisada
    while (i<=dim(arv.dados[[o]])[1]){
      aux<-NULL
      if( is.vector(arv.dados[[o+1]]))arv.dados[[o+1]]<-matrix(arv.dados[[o+1]],nrow=1)
      #olhe na coluna dos indicadores X = o(a-1)+(a-1)+ntotal+V+X
      if (arv.dados[[o]][i,(a*o- o +a +2)]==0 && dim(arv.dados[[(o+1)]])[1]>0){
        corte<-arv.dados[[o]][i,1:(o*(a-1))]
        for(oo in (o+1):(length(arv.dados))){
          if( is.vector(arv.dados[[oo]]) )arv.dados[[oo]]<-matrix(arv.dados[[oo]],nrow=1)
          aux<-matrix(rep(corte,dim(arv.dados[[oo]])[1]),ncol=length(corte),byrow=T)
          aux2<-matrix(arv.dados[[oo]][,1:o],ncol=o)
          corte2<-which(rowSums(aux2-aux)==0)
          arv.dados[[oo]]<-arv.dados[[oo]][-corte2,]
        }
      }
      if (arv.dados[[o]][i,(a*o- o +a +2)]==1){
        arv.dados[[o]]<-arv.dados[[o]][-i,]
        i<-i-1
        if( is.vector(arv.dados[[o]]))arv.dados[[o]]<-matrix(arv.dados[[o]],nrow=1)
      }
      i<-i+1
    }
  }
  arv.dados
}

