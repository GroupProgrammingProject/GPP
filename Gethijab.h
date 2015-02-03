#ifndef GETHIJAB_H
#define GETHIJAB_H

#include <iostream>
#include <cmath>

double Gethijab(int i, int j,int a, int b, double* d,int typei,int typej);

//Gethijab() returns value of non-scaled Hamiltionian matrix element h
double Gethijab(int i, int j,int a, int b, double* d,int typei,int typej){

int k; //for looping
double h,Es,Ep,V[4];														//h,Es,Ep and V[4] is only used locally in Gethijab()
double Es_C=-2.99,Ep_C=3.71;											//C orbital energies Es and Ep
double V_CC[4],V_CH[4],V_HH[4];
V_CC[0]=-5;V_CC[1]=4.7;V_CC[2]=5.5;V_CC[3]=-1.55;				//CC interaction 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi

if(typei==6 && typej==6){												//add more if statements, when we have interaction parameters for other atoms
	Es=Es_C;Ep=Ep_C;
	for(k=0;k<4;k++){V[k]=V_CC[k];}
}

//start V&G routine
if(i==j){
	if(a==b){
		if(a==0){h=Es;}
		else{h=Ep;}
	}
	else{h=0;}
}
else if(a*b==0){
	if(a==b){h=V[0];}														//ss_sigma
	else if(a==0){h=V[1]*d[b-1];}										//sp_sigma row
	else if(b==0){h=-V[1]*d[a-1];}										//sp_sigma column
}
else if(a==b){h=V[2]*pow(d[a-1],2)+V[3]-V[3]*pow(d[a-1],2);}		//pp_sigma and pp_pi diagonal
else{h=d[a-1]*d[b-1]*(V[2]-V[3]);}										//pp_sigma and pp_pi off-diagonal
//V&G routine ends

return h;
} //Gethijab() ends

#endif	
