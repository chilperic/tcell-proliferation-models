
mup0<-40
mud0<-60
mup<-8
mud<-60


mud01<-4.300041e+02
mud02<-1.000000e+03
mud03<-1.000000e+03 

mup01<-5.316612e+01
mup02<-4.707435e+01
mup03<-4.722025e+01

mud1<-2.341927e+01
mud2<-1.996041e+01
mud3<-2.217201e+01

mup1<-1.722716e+01
mup2<-1.375061e+01
mup3<-1.679847e+01 




AREmup0<-((abs(mup0-mup01)/abs(mup01)+abs(mup0-mup02)/abs(mup02)+abs(mup0-mup03)/abs(mup03))/3)*100

AREmud0<-((abs(mud0-mud01)/abs(mud01)+abs(mud0-mud02)/abs(mud02)+abs(mud0-mud03)/abs(mud03))/3)*100

AREmup<-((abs(mup-mup1)/abs(mup1)+abs(mup-mup2)/abs(mup2)+abs(mup-mup3)/abs(mup3))/3)*100

AREmud<-((abs(mud-mud1)/abs(mud1)+abs(mud-mud2)/abs(mud2)+abs(mud-mud3)/abs(mud3))/3)*100


AREmup0
AREmud0
AREmup
AREmud
