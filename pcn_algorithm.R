#    PCN Algorithm

#
# - Use matrices containing borders
# - Fix the formula for N = n^2
# - Maximum order of the tree becomes a function argument (easy to test different orders)

#### Order ####
# Calculate the observed value of the neighborhood of order o for the site ij
# with the new "mirror" argument to choose whether or not to mirror the matrix.
# m2 - matrix (of 2 (black) and 1 (white))
# i - row
# j - column
# o - order
# a - alphabet
# n - matrix size n x n
# Order adapted to deal with water and when a=2!!

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
## Modify the data matrix based on the observed neighborhood (update the number of sites and discovered neighborhoods).
## Pay attention to the use with the alphabet as we order it from 1 to max(alphabet).


obs.ordem<-function(m, dados, achou, obs, o, i, j, a, n, D = trunc(log(n/2))){
  # Se ja existe esse contexto
  if( length(achou)==1){
    dados[[o]][achou,(o*(a-1)+a)]<-dados[[o]][achou,(o*(a-1)+a)] + 1
      if(m[(i+D),(j+D)]<a){
        dados[[o]][achou,(o*(a-1)+m[(i+D),(j+D)])] <- dados[[o]][achou,(o*(a-1)+m[(i+D),(j+D)])]+1
    }
  }
  #"There is no such context."
  if(length(achou)==0){
    dados[[o]]<-rbind(dados[[o]], c(obs[1:(o*(a-1))],rep(0,(a-1)),1))
    if(m[(i+D),(j+D)]<a){
      dados[[o]][dim(dados[[o]])[1],(o*(a-1)+m[(i+D),(j+D)])] <- 1
    }
  }
  dados
}


#### arv.amostral ####
#Create a sample tree from a given matrix
#Adapted for water when a=2
#D = maximum order of the tree

arv.amostral<-function(m, n, a, D = trunc(log(n/2))){
 ## Enter a list that will be used as support
## The variable "dados" is a list. Each item in the list corresponds to a specific neighborhood order.
## If A=3: For neighborhood order 2 (for example), dados[[2]] will be a matrix with 9 columns, as we have 3 colors.
## The first 6 columns represent the first and second-order neighborhoods.
## The seventh and eighth columns correspond to the number of white and gray colors found.
## The ninth column represents how many times we observed this neighborhood, regardless of whether ij is white, gray, or black.
## Therefore, if the order is o and the alphabet (number of colors) is a, we have a*o+(a) columns.

  dados<-NULL
  ## This list will have log(n/2) elements = maximum order of neighborhood.

  dados<-vector("list", (D))
  obs<-NULL
# To find the context:
# Calculate the value of all neighbors for the site i=j=1 and allocate them in the list.

  for ( o in 1: D) {
    obs<-c(obs, ordem(m, 1, 1, o, a, n))
    dados[[o]]<-rbind(dados[[o]], c(obs[1:(o*(a-1))],rep(0,a)))
  }
  #For all other sites, do the same.

  for(i in 1:n){
    for ( j in 1: n){
      obs<-NULL
     #Check if the evaluated pixel is zero (water), if it is, skip to the next one.

        if(m[(i+D),(j+D)]!= 0){
           # Calculate the value of the neighborhood for each order.

            for ( o in 1: D){
            obs<-c(obs, ordem(m, i, j, o, a, n))
           # Try to find if this neighborhood already exists.
# List the number of different strings of order o (element "o" from the list).

            achou<-1: dim(dados[[o]])[1]
            #achou
            # At each order, keep only those that have the same neighborhood.
# To find a context, it must be tested from the root.

            for( p in 1:(o*(a-1))) achou<-achou[which(dados[[o]][achou,p]==obs[p])]
            #achou
           # Update the data with the new or already existing context.

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
# Calculate Pdj.
# Sometimes (nsa[a]/ns)^nsa[a] goes to zero, and the logarithm of that is -Inf.
# Using nsa[a]*log(nsa[a]/ns) avoids such issues.
# For example: log((844/2822)^844) = -Inf, but 844*log(844/2822) = -1018.749.

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
# Calculate the indicator Vdj.

Vdj<-function(dados, n){
### Calculation of Vs and Xs ###
# Add two columns.

  for ( o in 1: length(dados)) dados[[o]]<-cbind(dados[[o]],0,0)
 # First, do it for the last children.

  u<-length(dados)
  nfilhos<-dim(dados[[u]])[1]
  for( f in 1:nfilhos) dados[[u]][f,(o+3):(o+4) ]<-c(pdj(dados,u,f,n),0)
 # Now, do it for the rest.

  for( o in (u-1):1){
    # How many strings?
    nstring<-dim(dados[[o]])[1]
    for( s in 1:nstring){
      #s=1
      #What are the children of s?

      filhos<-1: dim(dados[[(o+1)]])[1]
      for (p in 1: o) filhos<- filhos[which(dados[[(o+1)]][filhos,p]==dados[[o]][s,p])]
     #What is the product of the V's of the children?

      #vprod<-prod(dados[[(o+1)]][filhos,(o+4)])
      vprod<-sum(dados[[(o+1)]][filhos,(o+4)])
     # Test to find V and X.

      dados[[o]][s,(o+3):(o+4)]<-c(max(pdj(dados,o,s,n),vprod),sum(vprod > pdj(dados,o,s,n))) #if true 1, if false 0
    }
  }
  dados
}

#### poda ####
# Prune the sample tree.

poda<-function(arv.dados,a){
# For all orders.

  for ( o in 1:(length(arv.dados)-1)){
    i=1   # auxiliar
 # Transform everything into a matrix.

    if( is.vector(arv.dados[[o]]))arv.dados[[o]]<-matrix(arv.dados[[o]],nrow=1)
    if(dim(arv.dados[[o]])[1]==0)i<-dim(arv.dados[[o]])[1]+1 ### In case it insists on being a vector
# While there is unanalyzed neighborhood.

    while (i<=dim(arv.dados[[o]])[1]){
      aux<-NULL
      if( is.vector(arv.dados[[o+1]]))arv.dados[[o+1]]<-matrix(arv.dados[[o+1]],nrow=1)
      # Look at the column of X indicators. = o(a-1)+(a-1)+ntotal+V+X
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

