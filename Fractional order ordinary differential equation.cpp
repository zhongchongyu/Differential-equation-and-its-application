#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
using namespace std;

double memo(double* r,double* c,int k)
{

	double temp = 0;
	for (int j=0; j<=(k-2); j++)
	{
		temp = temp + c[j]*r[k-j-2];
	}
	return temp;
}


int main(void)
{
	int i,j;
	double TSim = 500;
	double h = 0.005;
	int n = (int)(TSim/h + 1);

	//分数阶偏导:
	double q1 = 0.8;
	double q2 = 0.8;
	//Lotka-Volterra系统常数:
	double r = 0.5;
	double delta = 0.4;
	double beta = 0.8;
	//二项式系数的计算:
	double cp1 = 1;
	double cp2 = 1;
	double* c1 = (double*)malloc(sizeof(double)*n);
	double* c2 = (double*)malloc(sizeof(double)*n);
	for (j=0; j<=100000; j++)
	{
		c1[j] = (1-(1+q1)/(j+1))*cp1;
		c2[j] = (1-(1+q2)/(j+1))*cp2;

		cp1 = c1[j];
		cp2 = c2[j];
	}

	// 初始条件的设定:
	double* x = (double*)malloc(sizeof(double)*n);
	double* y = (double*)malloc(sizeof(double)*n);
	x[0] = 0.027;
	y[0] = 0.024;

	for (i=1; i<n; i++)
	{
		x[i] = ( r*x[i-1]*(1-x[i-1]) - sqrt(x[i-1])*y[i-1] )*pow(h,q1) - memo(x, c1, i+1);
		y[i] = ( beta*y[i-1]*(sqrt(x[i-1]) - delta))*pow(h,q2) - memo(y, c2, i+1);
	}

	
	
	ofstream outfile;
	outfile.open("result.txt");

	for (int i = 0; i < n; i++)
	{
		
		outfile << x[i] << " " << y[i] << endl;

	}
	outfile.close();
	
	/*
	for (i=0; i<n; i++)
	{
		printf("%f %f \n", x[i], y[i]);
	}
	*/
	return 0;
}