#ifndef HAMDER_H
#define HAMDER_H

#include <iostream>
#include <cmath>

double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr);

//Hamder() returns value of Hamilatonian matrix element differentiated wrt x,y or z.
double Hamder(int i, int j,int a, int b, std::vector<double>* d,double distr){

int k; //for looping
double h,Es,Ep,V[4];														//h,Es,Ep and V[4] is only used locally in Gethijab_der()
double Es_C=-2.99,Ep_C=3.71;											//C orbital energies Es and Ep
double V_CC[4];
V_CC[0]=-5;V_CC[1]=4.7;V_CC[2]=5.5;V_CC[3]=-1.55;				//CC interaction 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi

	Es=Es_C;Ep=Ep_C;						//set all parameters to the values for Xu's carbon
	for(k=0;k<4;k++){V[k]=V_CC[k];}

//start V&G routine
if(i==j){
	h=0;
}
else if(a*b==0){
	if(a==b){h=0;}										                                       //ss_sigma
	else if(a==0){h=V[1]*(1-pow((*d).at(b-1),2))/distr;}										//sp_sigma row
	else if(b==0){h=-V[1]*(1-pow((*d).at(a-1),2))/distr;}										//sp_sigma column
}
else if(a==b){h=2*(V[2]-V[3])*(1-pow((*d).at(a-1),2))/distr;}		                  //pp_sigma and pp_pi diagonal
else{h=(*d).at(b-1)*(V[2]-V[3])*(1-pow((*d).at(b-1),2))/distr;}							//pp_sigma and pp_pi off-diagonal
//V&G routine ends

return h;
} //Hamder() ends

	
#endif
