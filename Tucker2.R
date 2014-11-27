
## Copyright 2014 Giménez, Gustavo y Marticorena, Marta

## Author: Giménez, Gustavo y Marticorena, Marta

## Tucker2.R

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.



################################################################################
################################################################################
tucker2R <- function(datos,amb=2,stand=TRUE,nc1=2,nc2=2,niter=1000){
    amb==nc3
    if(!is.data.frame(datos)){stop("datos must be a data frame")}
    if(ncol(datos)%%amb!=0){stop("The variables must be the same in amb")}
    if (any(is.na(datos))){stop("There is at least one NA  'datos' must be complete")}
    if(amb>4){print("when amb > 4 The Tucker 3 is recommended!")}
    colnume <- ncol(datos)
    colnamb <- colnume/amb
    ambi <- amb-1
    for (j in 1:ambi){
        t <- -1+j
        t1 <- t*colnamb
        for (i in 1:colnamb){
            a <- i+t1
            b <- colnamb +i + t1
            fff <- colnames(datos)[a] == colnames(datos)[b]
            if (fff==FALSE){stop("datos must be the same variables")}
        }
    }

    if (stand == TRUE){datos.stan <- scale(datos)
                       center <- attributes(datos.stan)$`scaled:center`
                       Escala <- attributes(datos.stan)$`scaled:scale`
                       datos <- as.data.frame(datos.stan)
                   }else{
                       datos <- datos
                   }
    if(!(nc1<=nc2*nc3 && nc3<=nc1*nc2 && nc2 <= nc1*nc3)){stop("in the combination of components for a solution must use a Diffit criteria")}
    I<-dim(datos)[1]   
    J <-dim(datos)[2]  
    J <- J/amb 
    K<-amb   
    col <- J
    n1<-nc1 
    n2<-nc2 
    n3<-nc3 
   
    etiq_var <- colnames(datos)[1:J]
    etiq_varXamb <- rep(etiq_var,K)
    eti <- rep("amb",K)
    neti <- seq(1:K)
    etiq_amb <- paste(eti,neti,sep=".")
    etiq_ambXJ <- rep(etiq_amb,each=J)
    etiq <- paste(etiq_var,etiq_ambXJ,sep="-")
    etiq_ind <- rownames(datos)

    m <- array(0, dim=c(I,J,K))
    for (k in 1:K){
        n <- k-1
        l <- n*col
        for (i in 1:I)  
            {
                for(j in 1:J)
                    {
                        m[i,j,k]<-m[i,j,k]+datos[i,j+l];
                    }
            }
    }

    rownames(m) <- etiq_ind
    colnames(m) <- etiq_var
    dimnames(m)[3] <-list(etiq_amb)

    X1 <- array(0, dim=c(I,J*K))
    for (ii in 1:I)
        {
            c<-1; 
            for (kk in 1:K)
                {
                    for (jj in 1:J)
                        {
                            X1[ii,c]<-X1[ii,c]+m[ii,jj,kk];
                            c<-c+1
                        }
                }
        }

    X2 <- array(0, dim=c(J,I*K))
    for (jj in 1:J)
        {
            c<-1;  
            for (ii in 1:I)
                {
                    for (kk in 1:K)
                        {
                            X2[jj,c]<-X2[jj,c]+m[ii,jj,kk];
                            c<-c+1
                        }
                }
        }


    X3 <- array(0, dim=c(K,I*J))
    for (kk in 1:K)
        {
            c<-1; 
            for (jj in 1:J)
                {
                    for (ii in 1:I)
                        {
                            X3[kk,c]<-X3[kk,c]+m[ii,jj,kk];
                            c<-c+1
                        }
                }
        }


 
    W1=0 
    x=X1%*%t(X1)

    p<-0
    for (u in 1:I)
        {
            p<-p+x[u,u]
        }
    y=X2%*%t(X2);
    z=X3%*%t(X3);
    tam1 =nrow(X1);
    tam2 =ncol(X1);
    j1=1;

    a=svd(x)$u
    d1=svd(x)$d
    v=svd(x)$v

    b=svd(y)$u
    d2=svd(y)$d
    v=svd(y)$v

    a=a[,1:n1]
    b=b[,1:n2] 
    c=diag(n3) 

    iter <-  0 
    while (abs(j1)>=0.05){
        k1=kronecker(c,b) 
        G1=t(a)%*%X1%*%k1 
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l1=sum(diag(t1))
        m=l1
        j1=l1-W1
        k1=kronecker(c,b)
        x=X1%*%k1
        x=x%*%t(x) 
        a1=svd(x)$u
        d2=svd(x)$d
        v2= svd(x)$v
        a=a1[ ,1:n1]
        k1=kronecker(c,b)
        G1=t(a)%*%X1%*%k1
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l2=sum(diag(t1))
        m=l2
        j1=l2-l1
        k2=kronecker(c,a)
        y=X2%*%k2
        y=y%*%t(y)
        b1=svd(y)$u
        d3=svd(y)$d
        v3= svd(y)$v
        b=b1[ ,1:n2]
        k1=kronecker(c,b)
        G1=t(a)%*%X1%*%k1
        r1=kronecker(t(c),t(b))
        S1=a%*%G1%*%r1
        t1=(X1-S1)%*%t(X1-S1)
        l3=sum(diag(t1))
        m=l3
        j1=l3-l2
        iter=iter + 1
        if (iter >= niter) print("Warnings: not converg")
        if (iter >= niter) break
        W1=l3}
    A=a
    B=b
    C=c

    K=kronecker(C,B);
    D=G1%*%t(K);
    D=t(D);


    SCE=(p-m)/p * 100

    sca=0; 
    scb=0;
    sca=sum(sum(A^2));
    scb=sum(sum(D^2));
    sca=sca/tam1;
    scb=scb/tam2;
    scf=sqrt(sqrt(scb/sca));
    IND=a*scf
    IND
    INDi <- IND[,1:2]
    colnames(INDi) <- c("Dim-1","Dim-2")
    D=D/scf
    rownames(D) <- etiq
    rownames(IND) <- etiq_ind
    Dvar <- D[,1:2]
    colnames(Dvar) <- c("Dim-1","Dim-2")

    if (stand == TRUE){
        Resultado <- list(Scales=Escala,MEANS=center,individuos=etiq_ind,variables=etiq,
                          iteraciones=iter,SCExplicada=SCE,vartot=p)
    }
    if (stand == FALSE){
        Resultado <- list(individuos=etiq_ind,variables=etiq,
                          iteraciones=iter,SCExplicada=SCE,vartot=p)
    }

    Proyeccion <- list(variables=Dvar,Individuos=INDi)
    saltuck <- list(IND=IND,D=D,Ambientes=etiq_amb,
                    Resultados=Resultado,Proyecciones=Proyeccion)
    class(saltuck) <- "marta"
    return(saltuck)
}


plotear_tucker2 <- function(saltuck){
    if(class(saltuck)!="marta") stop("Object must be of class ’marta’")
    X1 <- saltuck$IND
    X <- as.data.frame(X1)
    Z1 <- saltuck$D
    Z <- as.data.frame(Z1)
    numcol <- dim(X)[2]
    numfil <- dim(X)[1]
    etiq <- rownames(saltuck$IND)
    etiq2 <- rownames(saltuck$D)
    ambsa <- saltuck$Ambientes
    ambs <- length(ambsa)
 
    plot(X[,1],
         X[,2],
         ylim=c(min(X,Z),max(X,Z)), 
         xlim=c(min(X,Z),max(X,Z)), type="n",
         ylab="DIM 2", xlab="DIM 1"
         )

    title(main = list("Interactive Biplot Tucker-v2", cex = 1.5,
              col = "red", font = 3))

    abline(h = 0, lty = 2); abline(v = 0, lty = 2) 

    for (i in 1:(numfil)) {
        points(X[i,1],X[i,2],pch=20, col="blue")
        text(X[i,1],X[i,2],labels=etiq[i],
             cex=0.6,pos=2)
    }

    dim(Z)[2] -> numcolz
    dim(Z)[1] -> numfilz

    ga <- 0:ambs
    fa <- ambs:1
    filXamb <- numfilz/ambs
    colores <- c("red","green","cyan","brown","coral","darkmagenta","hotpink")
    for (tt in 1:ambs){
        f <- 0
        gh <- 0
        f <- fa[tt]
        gh <- (filXamb)*ga[tt]
        gi <- ga[tt]+1
        for (k in (1+gh):(filXamb*gi)){
            arrows(0,0,Z[k,1],Z[k,2],length=0.05,
                   lwd=1.7, col=colores[tt])
            text(Z[k,1],Z[k,2],labels=etiq2[k],
                 cex=0.6,pos=2)
        }
    }
    legend("topright",legend=c(ambsa),fill=colores[1:ambs])

}
