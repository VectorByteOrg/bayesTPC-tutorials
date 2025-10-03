### 1. Functional forms and their temperature derivatives (d_...)

briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b
}

d_briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i]))
    )}
    else {b[i]<-0}
  }
  b
}

briere.trunc<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
    if(b[i]>1) b[i]<-1
  }
  b
}

d_briere.trunc<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i]))
    )}
    else {b[i]<-0}
    if(t[i]>T0 && t[i]<Tm && c*t[i]*(t[i]-T0)*sqrt(Tm-t[i])>1) b[i]<-0
  }
  b
}


BG<-function(t, c, E, x, s_squared, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*exp((-E/(8.617*10^-5*(t[i]+273.15)))-(((t[i]+273.15)-x)^2/(2*s_squared))))}
    else {b[i]<-0}
  }
  b
}

d_BG<-function(t, c, E, x, s_squared, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-((c*exp((-E/(8.617*10^-5*(t[i]+273.15)))-(((t[i]+273.15)-x)^2/(2*s_squared))))*(E/(8.617*10^-5*(t[i]+273.15)^2)-((t[i]+273.15)-x)/s_squared))}
    else {b[i]<-0}
  }
  b
}

boltz<-function(t, c, E, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*exp(-E/(8.617*10^-5*(t[i]+273.15))))}
    else {b[i]<-0}
  }
  b
}

d_boltz<-function(t, c, E, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-((c*exp(-E/(8.617*10^-5*(t[i]+273.15))))*(E/(8.617*10^-5*(t[i]+273.15)^2)))}
    else {b[i]<-0}
  }
  b
}

linear<-function(t, inter, slope){
  b=c()
  for (i in 1: length(t))
  {
    if  (inter+slope*t[i]>0) {b[i]<-inter+slope*t[i]}
    else {b[i]<-0}
  }
  b
}

d_linear<-function(t, inter, slope){
  b=c()
  for (i in 1: length(t))
  {
    if  (inter+slope*t[i]>0) {b[i]<-slope}
    else {b[i]<-0}
  }
  b
}


quad.alt<-function(t, inter, n.slope, qd){
  b=c()
  b.temp = c()
  for (i in 1:length(t)){
    b.temp[i] <- inter+n.slope*t[i]+qd*t[i]^2
    b[i] <- 0*(b.temp[i]<=0)+b.temp[i]*(b.temp[i]>0)
  }
  b
}

d_quad<-function(t, inter, slope, qd){
  b=c()
  for (i in 1:length(t)){
    if (inter+slope*t[i]+qd*t[i]^2>0) {b[i]<-slope+2*qd*t[i]}
    else {b[i]<-0}
  }
  b
}

quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<-qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
  }
  b
}

d_quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i]<-qd*(2*t[i]-T0-Tm)}
    else {b[i]<-0}
  }
  b
}

quad.2.trunc<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<-qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
    if(b[i]>1){ b[i] <-1 }
  }
  b
}

d_quad.2.trunc<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i]<-qd*(2*t[i]-T0-Tm)}
    else {b[i]<-0}
    if(qd*(t[i]-T0)*(t[i]-Tm)>1) b[i]<-0
  }
  b
}
