#include <iostream>
#include <cmath>

void near_neigh(int N, double *x, double *y, double *z, double rc, double *nnear, double *inear, double sx, double sy, double sz)
{
	double dist;
	for(int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			dx=x[i]-x[j];
			if(dx/sx>0.5)
			{
				dx=dx-sx;
			}
			dy=y[i]-y[j];
			dz=z[i]-z[j];
			dist=(dx*dx+dy*dy+dz*dz);
			if (dist<rc && i!=j)
			{
				nnear[i]=nnear[i]+1;
				inear
			}
		}
	}
}
