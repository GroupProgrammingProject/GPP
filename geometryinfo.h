#ifndef GEOMETRYINFO_H
#define GEOMETRYINFO_H

#include <vector>
#include <cmath>

// Function declarations
void DistanceComp(std::vector<double>* r, std::vector<double>* pos);
void GetDistanceComponents(std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz);
void GetAllDistances (std::vector<double>* modr, std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz);
void GetAllDistances (std::vector<double>* modr, std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc);


// Function to calculate one component of the distances between atoms
void DistanceComp(std::vector<double>* r, std::vector<double>* pos){
	int numatoms=pos->size();
	for (int i=0; i<numatoms; i++){
		for (int j=0; j<numatoms; j++){r->at(numatoms*i+j)=pos->at(i)-pos->at(j);}
	}
}

// Function to calculate distance components between atoms
void GetDistanceComponents(std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	DistanceComp(rx, posx);
	DistanceComp(ry, posy);
	DistanceComp(rz, posz);
}

// Function to calculate the distances between atoms
void GetAllDistances (std::vector<double>* modr, std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz){
	int numatoms=posx->size();
	GetDistanceComponents(rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms*numatoms; i++){
		double distx=rx->at(i);
		double disty=ry->at(i);
		double distz=rz->at(i);
		modr->at(i)=sqrt(distx*distx+disty*disty+distz*distz);
	}
}

// Function to calculate the distances between atoms and set all values above rc to zero
void GetAllDistances (std::vector<double>* modr, std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, std::vector<double>* posx, std::vector<double>* posy, std::vector<double>* posz, double rc){
	int numatoms=posx->size();
	GetDistanceComponents(rx, ry, rz, posx, posy, posz);
	for (int i=0; i<numatoms*numatoms; i++){
		double distx=rx->at(i);
		double disty=ry->at(i);
		double distz=rz->at(i);
		modr->at(i)=sqrt(distx*distx+disty*disty+distz*distz);
		if (modr->at(i)>rc){modr->at(i)=0;}
		else{;}
	}
}

// Function to calculate the number of nearest neighbours
void NumNearestNeighbours(std::vector<int>* nnear){
	
}


// Function that calculates the nearest neighbours of every atom
// Nearest neighbours are the neighbours within a radius rv
// inear is the vector in which the indices will be printed
// 
void NearestNeighbours(std::vector<int>* inear, std::vector<int>* nnear, std::vector<double>* rx, std::vector<double>* ry, std::vector<double>* rz, double rv){
	// Set all elements in inear and nnear to 0
	std::fill(inear->begin(), inear->end(), 0);
	std::fill(nnear->begin(), inear->end(), 0);

	int numatoms=rx->size();
	for (int i=0; i<numatoms; i++){
		for (int j =0; j<numatoms && j!=i; )
		
	}
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
