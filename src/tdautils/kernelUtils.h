#include <cmath>


using namespace std;

// read element of matrix
double ReadMat(double*XX, int *pNN, int *pDD, int i, int d){
	double out=0.0;
	out=XX[(d-1)*(pNN[0])+i-1];
	return out;
}

// write element of matrix
void WriteMat(double*XX, int *pNN, int *pDD, int i, int d, double input){
	XX[(d-1)*(pNN[0])+i-1]=input;
}


// oneKernel
double oneKernel(double *pp, double* XX, int *pNN, int *pDD, double *hh){
	double sum=0.0;
	double tmp=0.0;	
	double out=0.0;	
	
	for (int i=1; i<=pNN[0]; i++){
		sum=0.0;
		for (int d=1; d<=pDD[0]; d++){
			tmp= pp[d-1] - ReadMat(XX, pNN, pDD, i, d);
			sum=sum+ tmp*tmp;		
		}
		sum= exp( - sum/ (2*hh[0]*hh[0]));
	out=out+sum;
	}
	out=out/pNN[0];
	return out;
}

