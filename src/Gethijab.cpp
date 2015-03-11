#include "../include/Gethijab.h"

//Gethijab() returns value of non-scaled Hamiltionian matrix element h
double Gethijab(int i, int j,int a, int b, std::vector<double>* d,std::vector<double>* TBparam){

int k; //for looping
double h,Es,Ep,V[4];														//h,Es,Ep and V[4] is only used locally in Gethijab()
Es=TBparam->at(0);
Ep=TBparam->at(1);
//CC interaction in V: 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi
for(k=0;k<4;k++){V[k]=TBparam->at(k+2);}

//start V&G routine
if(i==j){
	if(a==b){
		if(a==0){h=Es;}
		else{h=Ep;}
	}
	else{h=0;}
}
else if(a*b==0){
	if(a==b){h=V[0];}											//ss_sigma
	else if(a==0){h=V[1]*(*d).at(b-1);}					//sp_sigma row
	else if(b==0){h=-V[1]*(*d).at(a-1);}				//sp_sigma column
}
else if(a==b){h=V[2]*pow((*d).at(a-1),2)+V[3]-V[3]*pow((*d).at(a-1),2);}		//pp_sigma and pp_pi diagonal
else{h=(*d).at(a-1)*(*d).at(b-1)*(V[2]-V[3]);}											//pp_sigma and pp_pi off-diagonal
//V&G routine ends

return h;
} //Gethijab() ends

	
