#include <iostream>
#include <vector>
#include <cmath>
#include <string>

//double Gethijab(int a, int b, double* d,double r,std::string typei,std::string typej);
double Gethijab(int a, int b, double* d,double r);	//without types

//Gethijab() takes a,b,d and types of i and j and returns value of non-scaled Hamiltionian matrix element
//double Gethijab(int a, int b, double* d,double r,std::string typei,std::string typej){
double Gethijab(int a, int b, double* d,double r){	//without types

double h,r_cutoff,V_CC[4];	//parameters

//here would be an if statement so that typei=typej="C"
V_CC[0]=-5;V_CC[1]=4.7;V_CC[2]=5.5;V_CC[3]=-1.55;
r_cutoff=2.6;	//where do I read this from???

//start V&G routine
if(r<r_cutoff){h=0;}
else if(a*b==0;){h=V_CC[0];}
else

h=1;

return h;
} //get_hijab() ends

	
int main() {	//main does what hamiltonian.cpp would ask of Gethijab()

int i,j,a,b;
double h,d[3],r=1.73;
//char typei="C",typej="C";

for(i=0;i<3;i++){d[i]=1;}	//dummy d vector

for(a=0;a<4;a++){
	for(b=0;b<4;b++){
		std::cout << "a is " << a << " and b is " << b << std::endl;
//		h=Gethijab(a,b,d,r,typei,typej);
		h=Gethijab(a,b,d,r);
		std::cout << "h=" << h << std::endl;
	} //end b loop
} //end a loop
return 0;
} //main ends


	
