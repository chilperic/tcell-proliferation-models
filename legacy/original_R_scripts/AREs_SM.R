
deltat=1.000000e-03 
mup0<-40
mud0<-60
mup<-8
mud<-60
deltat=1.000000e-03 
lambda_01<-1.870167e-02
lambda_02<-1.930920e-02
lambda_03<-2.008399e-02 

lambda1<-3.422586e-02
lambda2<-4.924785e-02
lambda3<-3.834159e-02

d01<-1.000000e-03
d02<-1.000000e-03
d03<-1.000000e-03

d1<-3.425349e-02
d2<-4.327687e-02
d3<-3.673752e-02

mud01<-1/d01
mud02<-1/d02
mud03<-1/d03
  
mup01<-(1/lambda_01)+deltat
mup02<-(1/lambda_02)+deltat
mup03<-(1/lambda_03)+deltat

mud1<-1/d1
mud2<-1/d2
mud3<-1/d3


mup1<-(1/lambda1)+deltat
mup2<-(1/lambda2)+deltat
mup3<-(1/lambda3)+deltat

AREmup0<-((abs(mup0-mup01)/abs(mup01)+abs(mup0-mup02)/abs(mup02)+abs(mup0-mup03)/abs(mup03))/3)*100

AREmud0<-((abs(mud0-mud01)/abs(mud01)+abs(mud0-mud02)/abs(mud02)+abs(mud0-mud03)/abs(mud03))/3)*100

AREmup<-((abs(mup-mup1)/abs(mup1)+abs(mup-mup2)/abs(mup2)+abs(mup-mup3)/abs(mup3))/3)*100

AREmud<-((abs(mud-mud1)/abs(mud1)+abs(mud-mud2)/abs(mud2)+abs(mud-mud3)/abs(mud3))/3)*100


AREmup0
AREmud0
AREmup
AREmud
