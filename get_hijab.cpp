#include <iostream>
#include <vector>
#include <cmath>
#include <string>

double Gethijab(int a, int b, double* d,double r,int typei,int typej);

//std::string pointer
//const std::string& name

//Gethijab() returns value of non-scaled Hamiltionian matrix element h
double Gethijab(int a, int b, double* d,double r,int typei,int typej){

int i;
double h,V[4];																//h and V[4] is only used locally in Gethijab()
double r_cutoff=2.6;														//where should I take r_cutoff from?
double V_CC[4],V_CH[4],V_HH[4];
V_CC[0]=-5;V_CC[1]=4.7;V_CC[2]=5.5;V_CC[3]=-1.55;				//CC interaction 0=ss_sigma, 1=sp_sigma, 2=pp_sigma, 3=pp_pi

if(typei==6 && typej==6)											//add more if statements, when we have interaction parameters for other atoms
{for(i=0;i<4;i++){V[i]=V_CC[i];}}

//start V&G routine
if(r>r_cutoff){h=0;}														//check if r within cut-off radius
else if(a*b==0){
	if(a==b){h=V[0];}														//ss_sigma
	else if(a==0){h=V[1]*d[b-1];}										//sp_sigma row
	else if(b==0){h=-V[1]*d[a-1];}										//sp_sigma column
}
else if(a==b){h=V[2]*pow(d[a-1],2)+V[3]-V[3]*pow(d[a-1],2);}		//pp_sigma and pp_pi diagonal
else{h=d[a-1]*d[b-1]*(V[2]-V[3]);}										//pp_sigma and pp_pi off-diagonal
//V&G routine ends

return h;
} //get_hijab() ends

	
int main() {	//main does what hamiltonian.cpp would ask of Gethijab()

//input from hamiltonian to Gethijab()
int i,j,a,b;
double h,d[3],r=1.73;
int typei=6,typej=6;
for(i=0;i<3;i++){d[i]=1;}	//dummy d vector
//end listing input

for(a=0;a<4;a++){
	for(b=0;b<4;b++){
		h=Gethijab(a,b,d,r,typei,typej);
		std::cout << h << "	";
	} //end b loop
std::cout << std::endl;
} //end a loop

return 0;
} //main ends


	
