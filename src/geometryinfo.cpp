#include "../include/geometryinfo.h"

// Function to calculate one component of the distances between atoms
void DistanceComp(Eigen::MatrixXd* r, std::vector<double>* pos){
	int numatoms=pos->size();
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){(*r)(i, j)=pos->at(i)-pos->at(j);}
	}
}

// Function to calculate distance components between atoms
void GetDistanceComponents(Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	DistanceComp(rx, posx);
	DistanceComp(ry, posy);
	DistanceComp(rz, posz);
}


// Function to calculate the distances between atoms
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	int numatoms=posx->size();
	GetDistanceComponents(rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
		double distx=(*rx)(i, j);
		double disty=(*ry)(i, j);
		double distz=(*rz)(i, j);
		(*modr)(i, j)=sqrt(distx*distx+disty*disty+distz*distz);
		}
	}
}

// Function to calculate the distances between atoms and set all values above rc to zero
void GetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc){
	int numatoms=posx->size();
	GetAllDistances(modr, rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
			if ((*modr)(i,j)>rc){(*modr)(i, j)=rc+3;}
		}
	}
}

// Function to calculate one component of the distances between atoms
void PbcDistanceComp(Eigen::MatrixXd* r, std::vector<double>* pos, double latticeconst, double rv){
	int numatoms=pos->size();
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
			double result=rv+3.;
			double xi = pos->at(i), xj=pos->at(j);
			double dist=xi-xj;
			double dist1= xi-xj-latticeconst;
			double dist2= xi-xj+latticeconst;
			if (Abs(dist)<rv){result=dist;}
			else if (Abs(dist1)<rv){result=dist1;}
			else if (Abs(dist2)<rv){result=dist2;}
			else{;}
			(*r)(i, j)=result;
		}
	}
}

// Function to calculate distance components between atoms
void PbcGetDistanceComponents(Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double a, double b, double c, double rv){
	int numatoms=posx->size();
	PbcDistanceComp(rx, posx, a, rv);
	PbcDistanceComp(ry, posy, b, rv);
	PbcDistanceComp(rz, posz, c, rv);
}

// Function which may be always called and chooses between PBC or conventional GetAllDistances
void GetDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* lats, double rv, bool pbc){
  if (pbc == 1) {PbcGetAllDistances(modr,rx,ry,rz,posx,posy,posz,lats,rv);}
  else {GetAllDistances(modr,rx,ry,rz,posx,posy,posz,rv);}
}

// Function to calculate the distances between atoms
void PbcGetAllDistances (Eigen::MatrixXd* modr, Eigen::MatrixXd* rx, Eigen::MatrixXd* ry, Eigen::MatrixXd* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, std::vector<double>* lats, double rv){
	int numatoms=posx->size();
	double a = (*lats).at(0);
	double b = (*lats).at(1);
	double c = (*lats).at(2);
	PbcGetDistanceComponents(rx, ry, rz, posx, posy, posz, a, b, c, rv);
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
		double distx=(*rx)(i, j);
		double disty=(*ry)(i, j);
		double distz=(*rz)(i, j);
		double modulusr=sqrt(distx*distx+disty*disty+distz*distz);
		if (modulusr< rv){(*modr)(i, j)=modulusr;}
		else {(*modr)(i, j)=rv+3;}
		}
	}
}



// Function that calculates the nearest neighbours of every atom
// Nearest neighbours are the neighbours within a radius rv
// inear is the vector in which the indices will be printed
void NearestNeighbours(Eigen::MatrixXi* inear, std::vector<int>* nnear, Eigen::MatrixXd* modr, double rv){
	int numatoms=nnear->size();
	for (int i=0; i<numatoms; i++){
		int nnearcounter=0;
		int k=0;
		for (int j =0; j<numatoms; j++){
			double dist = (*modr)(i, j);
			if (dist<rv && dist!=0){
				nnearcounter++;
				(*inear)(i, k)=j;
				k++;
			}
			else {}	
		}
		nnear->at(i)=nnearcounter;
	}
}


// Simple absolute value function
double Abs(double c){
	double result=c;
	if (c<0){result=-c;}
	else{;}
return result;
}

// Function to set all elements of 2 equal length vectors equal too each other
// equivalent to A=B for std::vector<double> objects
void SetEqual(std::vector<double>* a, std::vector<double>* b){
	int n=a->size();
	for (int i=0; i<n; i++){a->at(i)=b->at(i);}
}

// Function that decides whether the Nearest neighbours lists need to be recalculated between MD steps
bool RecalculateNearestNeighbours(std::vector<double>* refposx, std::vector<double>* refposy, std::vector<double>* refposz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc, double rv){
	int numatoms=posx->size();
	bool recalc=0;
	std::vector<double> disp(numatoms);
	for (int i=0; i<numatoms; i++){
		double xdisp=posx->at(i)-refposx->at(i);
		double ydisp=posy->at(i)-refposy->at(i);
		double zdisp=posz->at(i)-refposz->at(i);
		disp.at(i)=sqrt(xdisp*xdisp+ydisp*ydisp+zdisp*zdisp);
	}
	double dmax=*(max_element(disp.begin(), disp.end()));
	//std::cout<<dmax<<std::endl;
	if (dmax>0.4*Abs(rc-rv)){
		SetEqual(refposx, posx);
		SetEqual(refposy, posy);
		SetEqual(refposz, posz);
		recalc=1;
	}
	else{;}
return recalc;
}

//Translates all atom positions to within unit cell
void pbcshift(std::vector<double>* x, std::vector<double>* y, std::vector<double>* z, std::vector<double>* lats)
{
	int N=(*x).size(); //number of atoms
	for(int i=0; i<N; i++){
		if((*x).at(i)>(*lats).at(0)/2){ //the cell is centred on (0,0,0) with side lengths (lats.at(0),lats.at(1),lats.at(2))
			while((*x).at(i)>(*lats).at(0)/2){(*x).at(i)=(*x).at(i)-(*lats).at(0);}
		}
		else if((*x).at(i)<-(*lats).at(0)/2){
			while((*x).at(i)<-(*lats).at(0)/2){(*x).at(i)=(*x).at(i)+(*lats).at(0);}
		}	
		if((*y).at(i)>(*lats).at(1)/2){
			while((*y).at(i)>(*lats).at(1)/2){(*y).at(i)=(*y).at(i)-(*lats).at(1);}
		}
		else if((*y).at(i)<-(*lats).at(1)/2){
			while((*y).at(i)<-(*lats).at(1)/2){(*y).at(i)=(*y).at(i)+(*lats).at(1);}
		}
		if((*z).at(i)>(*lats).at(2)/2){
			while((*z).at(i)>(*lats).at(2)/2){(*z).at(i)=(*z).at(i)-(*lats).at(2);}
		}
		else if((*z).at(i)<-(*lats).at(2)/2){
			while((*z).at(i)<-(*lats).at(2)/2){(*z).at(i)=(*z).at(i)+(*lats).at(2);}
		}
	}
}

