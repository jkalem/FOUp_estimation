

 ##This program estimate the parameters in a FOU(3) process, that is FOU(l1,l2,l3) (with the different lambda's).
  ##Input x: the observed time series. It is mandatory to remove any trend or seasonal component
  ##Input aa: parameter a in the weight function, according with Theorem 4.
  ##(aa must be greater than or equal to 2p=6.
  ##Input bb: here the parameter b of the weight function according to Theorem 4 is b= aa+bb 
  ##(bb must be greater tha or equal to 3)
a=c(0.4829629131445341,-0.8365163037378077,0.2241438680420134,0.1294095225512603)/sqrt(2)
#In the following lines,  the binomial filters of different orders can be used instead de Daubechies filter in previous line. 
#a=c(-1,2,-1)/4
#a=c(-1,3,-3,1)/8
#a=c(1,-4,6,-4,1)/16
#a=c(-1,5,-10,10,-5,1)/32
#a=c(-1,6,-15,20,-15,6,-1)/64
#a=c(-1,7,-21,35,-35,21,-7,1)/128
#  a=c(-1,8,-28,56,-70,56,-28,8,-1)/256

  H_o_l1l2l3=function(x,Te,aa,bb)
  {
  x=x-mean(x)
  n=length(x)
  


  
 
  p=length(a)
  del=Te/n
  
  v1=rep(0,(length(x)-p+1))
  v2=rep(0,(length(x)-2*p))
  #Ahora calculamos el filtro de longitud 2p
  a2=rep(0,p)
  a2=rbind(a,a2)
  a2=a2[1:(2*p-1)]
  
  
  for(i in 1:(length(x)-p+1))
  {
    v1[i]=(sum(a*x[i:(i+p-1)]))^2
  }
  va=mean(v1)
  
  for(i in 1:(length(x)-2*p+2))
  {
    v2[i]=(sum(a2*x[i:(i+2*p-2)]))^2
  }
  va2=mean(v2)
  
  Hg=log2(va2/va)/2# Hg=Hgorro
  b=rep(0,length(a))
  
  for(j in 1:length(a))
  {
    A=seq(1,p)
    A=abs(A-j)^(2*Hg)
    b[j]=sum(a*A)
  }
  
  og=(-2*va/(del^(2*Hg)*sum(a*b)))^0.5# sigma gorro (ogorro)
  

n=length(x)

  U_Tn=function(beta)
  {
  l1=beta[1]
  l2=beta[2]
  l3=beta[3]
    j=seq(Te/n,Te,Te/n)
    w=abs(j)^aa/(1+abs(j)^(aa+bb))
    esp_dens=og^2*gamma(2*Hg+1)*sin(Hg*pi)*abs(j)^(5-2*Hg)/(2*pi*(l1^2+j^2)*(l2^2+j^2)*(l3^2+j^2))    
    
    Idelta=rep(NA,n)
    for (i in 1:n)
    {
      Idelta[i]=Te/(2*pi)*(mean(x*cos(j*i*Te/n))^2+mean(x*sin(j*i*Te/n))^2)
    }
    (Te/(2*pi*n))*sum((log(esp_dens)+Idelta/esp_dens)*w)
    
  }
    lg=constrOptim(c(0.3,0.5,0.8),U_Tn,NULL, ui=rbind(c(1,0,0),c(-1,1,0),c(0,-1,1)),ci=c(0.01,0.01,0.01))$par
 #lg=optim(c(0.3,0.5,0.8),U_Tn, lower=c(0.01,.015,0.02),method ="L-BFGS-B")$par
 lg=sort(lg)
 
  l1g=lg[1]
  l2g=lg[2]
  l3g=lg[3]
  print("Estimation of parameters: H, sigma, lambda_1, lambda_2 and lambda_3")
  c(Hg,og,l1g,l2g,l3g) 
  }
 
  
  