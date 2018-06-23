//显式向前欧拉法求解抛物型偏微分方程：A*(u_xx)=(u_t)
//其中A=0.5
//x的取值范围：(0,xf) = (0,2)
//t的取值范围：(0,T) = (0,0.1)
#include <stdio.h>
#include <stdlib.h>				//动态分配数组
#include <math.h>				//求sin值
#include<fstream>
using namespace std;

#define A 0.5					//给定条件
#define xf 2
#define T 0.1
#define M 50
#define N 100
#define pi 3.14159265358979323846

int main(void)
{

	double dx = (double)xf / M, dt = (double)T / N;				//计算dx和dt
	int i; int k;
	double* it0 = (double*)malloc(sizeof(double)* (M + 1));	//动态分配数组it0（初始条件）
	for (i = 0; i <= M; i++)
	{
		*(it0 + i) = sin(pi * (i*dx));
	}

	double* bx0 = (double*)malloc(sizeof(double)* (N + 1));	//动态分配数组bx0和bxf（边界条件，本例中均为0）
	double* bxf = (double*)malloc(sizeof(double)* (N + 1));
	for (i = 0; i <= N; i++)
	{
		*(bx0 + i) = 0 * (i*dt);
		*(bxf + i) = 0 * (i*dt);
	}

	double** u;												//动态分配二维数组u的空间（用来存放函数的值）
	u = (double**)malloc(sizeof(double)*(M + 1));
	for (i = 0; i <= M; i++)
	{
		*(u + i) = (double*)malloc(sizeof(double)*(N + 1));
	}


	for (i = 0; i <= M; i++)			//将三个边界条件it0,bx0,bxf存入u矩阵中
	{
		u[i][0] = it0[i];
	}
	for (i = 0; i <= N; i++)
	{
		u[0][i] = bx0[i];
		u[M][i] = bxf[i];
	}
	double r = A*((dt / dx) / dx);			//计算r和r1
	double r1 = 1 - 2 * r;

	for (int k = 0; k <= N; k++)		//向前欧拉法计算u矩阵
	{
		for (i = 1; i<M; i++)
		{
			u[i][k + 1] = r*(u[i + 1][k] + u[i - 1][k]) + r1 * u[i][k];
		}
	}
	ofstream outfile;
	outfile.open("resule.txt");

	for (int i = 0; i < M; i++)
	{
		for (int k = 0; k < N; k++)
		{
			outfile << u[i][k] << " ";
		}
		outfile << endl;

	}
	outfile.close();

}