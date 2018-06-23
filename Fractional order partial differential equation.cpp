#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>



#define  pi 3.1415926

using namespace std;

//创建可变的二维数组指针
template <typename T>
T** new_Array2D(int row, int col);
//释放二维数组指针
template <typename T>
void delete_Array2D(T **arr, int row, int col);

//  inline function declarations
inline double f1(double, double); // reactoin function of prey
inline double f2(double, double); // reactoin function of predator

// memory function declarations
double memo(double *, double *, int);


// global parameter's defines
double r = 0.5, delta = 0.4, beta = 0.8;// other parameters 

int main()
{
	// initialization of parameters 
	const double Tsim = 0.5;// Length of simulation time
	const double Tstep = 0.2;//Time step
	const int col = static_cast<int>(Tsim / Tstep + 1); /*number of calculated mesh points of time axis,
														alternatively col number of time--space mesh matrix */
	const double L = 1;   // Length of spatial region
	const int row = 6; /* number of calculated mesh points of space axis,
					   alternatively row number of time--space mesh matrix */
	const double Xstep = (double)(L / (row - 1));    // Step length of space division
	double q1 = 0.8, q2 = 0.8;//orders of derivatives, respectively
	double d1 = 0.1, d2 = 0.2;//diffusion coefficients of prey and predator 
	// Length of one-dimension space and mesh 

	//分配 u、v变量的时空储存二维空间
	double **u = new_Array2D<double>(row, col);// attribute memory  of prey
	double **v = new_Array2D<double>(row, col); // attribute memory  of prey;// attribute memory  of predator
	// u和v的初值
	for (int i = 0; i < row; i++)
	{
		u[i][0] = 0.027;
		v[i][0] = 0.024;
	}

	//calculation of phase portraits /numerical solution/:
	double cp1 = 1, cp2 = 1;
	double *c1 = new double[col], *c2 = new double[col];
	for (int j = 0; j < col; j++)
	{
		c1[j] = (1 - (1 + q1) / (j + 1))*cp1;
		c2[j] = (1 - (1 + q2) / (j + 1))*cp2;
		cp1 = c1[j]; cp2 = c2[j];
	}

	// main loop
	for (int i = 1; i < col; i++)// loop of time axis
	{
		// loop of left boundary points 
		u[0][i] = (2 * d1*(u[1][i - 1] - u[0][i - 1]) / pow(Xstep, 2) + f1(u[0][i - 1], v[0][i - 1]))*pow(Tstep, q1) - memo(u[0], c1, i+1);
		v[0][i] = (2 * d2*(v[1][i - 1] - v[0][i - 1]) / pow(Xstep, 2) + f2(u[0][i - 1], v[0][i - 1]))*pow(Tstep, q2) - memo(v[0], c2, i+1);
		//loop of right boundary points
		u[row - 1][i] = (2 * d1*(u[row - 2][i - 1] - u[row - 1][i - 1]) / pow(Xstep, 2) + f1(u[row - 1][i - 1], v[row - 1][i - 1]))*pow(Tstep, q1) - memo(u[row - 1], c1, i+1);
		v[row - 1][i] = (2 * d2*(v[row - 2][i - 1] - v[row - 1][i - 1]) / pow(Xstep, 2) + f2(u[row - 1][i - 1], v[row - 1][i - 1]))*pow(Tstep, q2) - memo(v[row - 1], c2, i+1);
		//  loop of interior region
		for (int k = 1; k < (row - 1); k++) // loop of space axis
		{
			u[k][i] = (d1*(u[k + 1][i - 1] - 2 * u[k][i - 1] + u[k - 1][i - 1]) / pow(Xstep, 2) + f1(u[k][i - 1], v[k][i - 1]))*pow(Tstep, q1) - memo(u[k], c1, i+1);
			v[k][i] = (d2*(v[k + 1][i - 1] - 2 * v[k][i - 1] + v[k - 1][i - 1]) / pow(Xstep, 2) + f2(u[k][i - 1], v[k][i - 1]))*pow(Tstep, q1) - memo(v[k], c2, i+1);

		}
	}

	
	
	
	ofstream outfile;
	outfile.open("result_u.txt");
	for (int i = 0; i < row; i++)
	{
		for (int k = 0; k < col; k++)
		{
			outfile << u[i][k] << " ";
		}
		outfile << endl;

	}
	outfile.close();
	
	
	outfile.open("result_v.txt");
	for (int i = 0; i < row; i++)
	{
		for (int k = 0; k < col; k++)
		{
			outfile << v[i][k] << " ";
		}
		outfile << endl;

	}
	outfile.close();
	

	/*
	cout << c1[0] << endl;
	cout << u[0][0] << endl;
	*/

	delete_Array2D(u, row, col);
	delete_Array2D(v, row, col);
	delete[] c1; c1 = nullptr;
	delete[] c2; c2 = nullptr;
	return 0;
	system("pause");

}

//// memory function define
double memo(double *r, double *c, int k)
{
	double temp = 0;
	for (int j = 0; j < k - 1; j++)
		temp = temp + c[j] * r[k - j -2];
	return temp;
}

// reactoin field of prey
inline double f1(double u, double v)
{
	double result;
	result = r*u*(1 - v) - sqrt(u)*v;
	return result;
}
// reactoin field of predator
inline double f2(double u, double v)
{
	double result;
	result = beta*v*(sqrt(u) - delta);
	return result;
}






template <typename T>
T** new_Array2D(int row, int col)
{
	int size = sizeof(T);
	int point_size = sizeof(T*);
	//先申请内存，其中sizeof(T*) * row表示存放row个行指针
	T **arr = (T **)malloc(point_size * row + size * row * col);
	if (arr != NULL)
	{
		T *head = (T*)((int)arr + point_size * row);
		for (int i = 0; i < row; ++i)
		{
			arr[i] = (T*)((int)head + i * col * size);
			for (int j = 0; j < col; ++j)
				new (&arr[i][j]) T;
		}
	}
	return (T**)arr;
}
//释放二维数组
template <typename T>
void delete_Array2D(T **arr, int row, int col)
{
	for (int i = 0; i < row; ++i)
	for (int j = 0; j < col; ++j)
		arr[i][j].~T();
	if (arr != NULL)
		free((void**)arr);
}