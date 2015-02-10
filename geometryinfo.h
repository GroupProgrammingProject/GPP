#ifndef GEOMETRYINFO_H
#define GEOMETRYINFO_H

#include <vector>
#include <cmath>

// Function to calculate the number of nearest neighbours
void NumNearestNeighbours(std::vector<double>* nnear, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	
}

// Function that calculates the nearest neighbour of every atom
// inear is the vector in which the indices will be printed
void NearestNeighbours(std::vector<double>* inear, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rv){
	// Delete all elements in inear
	inear->clear();
	int numatoms=posx->size();
	for (int i=0; i<numatoms; i++){
		
	}
}

// Function to generate sign for distances
// Sign convention to have distances signed according to sign a-b, ie ij element has sign pos[i]-pos[j]
double GetSign(double a, double b){
	double result=1;
	if (a-b<0){result=-1;}
	return(result);
}

// Function to calculate one component of the distances between atoms
void DistanceComp(std::vector<double>* r, std::vector<double>* pos){
	int numatoms=pos->size();
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){
			double sign = GetSign(pos->at(i), pos->at(j));
			r->at(numatoms*i+j)=sign*sqrt((pos->at(i)-pos->at(j))*(pos->at(i)-pos->at(j)));
		}
	}
}

// Function to calculate distance components between atoms
void GetDistanceComponents(std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	DistanceComp(rx, posx);
	DistanceComp(ry, posy);
	DistanceComp(rz, posz);
}

/*

void near_neigh(int N, double *x, double *y, double *z, double rc, int *nnear, int *inear, double sx, double sy, double sz)
{  //determine the nearest neighbours for each atom
	double dx,dy,dz,dist;
	for (int i=0; i<N; i++) { nnear[i]=0; }
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dx=x[i]-x[j];
			dx=dx-sx*round(dx/sx);//if dx>sx, the distance to N.N. is dx-sx
			dy=y[i]-y[j];
			dy=dy-sy*round(dy/sy);
			dz=z[i]-z[j];
			dz=dz-sz*round(dz/sz);
			dist=sqrt(dx*dx+dy*dy+dz*dz);
			if (dist<rc && i!=j) //add the atom j to the nearest neighbour list of i if this holds
			{
				nnear[i]=nnear[i]+1;
				inear[i*N+nnear[i]]=j; //a matrix with i rows, nnear[i] (no of nearest neighbours) columns
			}
		}
	}
}
*/

#endif
