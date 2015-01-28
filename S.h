double S(double r, int typei, int typej) {
  double r0=1.536329; //nearest-neighbour atomic separation
  double n=2;
  double nc=6.5;
  double rc=2.18;
  double r1=2.45;
  double S;
  //constants 

  if (r<r1){
    S=pow(r0/r,n)*exp(n*(-pow(r/rc,nc)+pow(r0/rc,nc)));
  }
  else{
    S=ts(r);
  }
  return S;
}
