/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Householder Reflection QR Decomposition
	Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

	Solve Ax = B(i..n)

	Input:
	N M
	a11 a12 a13 ... a1n b11 b12 b13 ... b1m
	a21 a22 a23 ... a2n b21 b22 b23 ... b2m
	 .   .   .       .   .   .   .       .
	an1 an2 an3 ... ann bn1 bn2 bn3 ... bnm 

	Example:
	3 3
	1 2 3 10 11 12
	4 5 6 13 14 15
	7 8 9 16 17 18

	Solves:
	1x+2y+3z = 10
	4x+5y+6z = 13
	7x+8y+9z = 16

	1x+2y+3z = 11
	4x+5y+6z = 14
	7x+8y+9z = 15

	1x+2y+3z = 12
	4x+5y+6z = 15
	7x+8y+9z = 18
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

int N, M;
void printMatrix(double_matrix, string);
void printVector(double_vector, string);
void printVector(double_matrix, int, string);
double_matrix multiplyMatrix(double_matrix, double_matrix);
double_matrix transposeMatriz(double_matrix);
double_matrix identity(int);
double_matrix matrixMinusMatrix(double_matrix, double_matrix);
double_matrix numberMultiplyMatrix(double_matrix, double);
double_matrix expandMatrix(double_matrix, int);
double norm(double_vector);
double innerProductSpace(double_vector, double_vector);
double signal(double);
const clock_t begin_time = clock();

int main(){
	double_matrix A, H, Q, R, B, X;
	double el;

	// NxN matriz with M results
	cin >> N >> M;

	for(int i = 0; i < N; i++)
	{
		double_vector row;
		A.push_back(row);
		B.push_back(row);
		X.push_back(row);
		
		for(int j = 0; j < N; j++)
		{
			cin >> el;
			A[i].push_back(el);
		}

		for(int j=0; j < M; j++)
		{
			cin >> el;
			B[i].push_back(el);
			X[i].push_back(0);
		}
	}

	printMatrix(A, "A");

	Q = identity(N);

	int i = 0;
	while(true)
	{	
		if(N-i < 2)
		{
			break;
		}

		// Ai
		double_vector Ai;
		for(int j = i; j < N; j++)
		{
			Ai.push_back(A[j][i]);
		}

		// Sign and Norm
		double sign = signal(A[i][i]);
		double normValue = norm(Ai);

		cout << "Ai = (" << endl;
		for(int l = 0; l < Ai.size() ; l++)
		{
			cout << Ai[l] << endl;
		}
		cout << ")" << endl;

		cout << "Signal = " << sign << endl;
		cout << "Norm = " << normValue << endl;

		// Vi = Ai - sign * norm * e1 (1, 0, . . . . ., 0)
		double_matrix Vi;
		for(int j = 0; j < Ai.size(); j++)
		{
			double_vector row;
			if(j == 0){
				row.push_back(Ai[j] - (sign * normValue));
			}
			else
			{
				row.push_back(Ai[j]);
			}
			Vi.push_back(row);
		}

		printMatrix(Vi, "Vi");

		// Calculating Hi. Need I, Vi in Matrix and Vector Formats
		double multiple = 2/innerProductSpace(transposeMatriz(Vi)[0], transposeMatriz(Vi)[0]);
		cout << "Multiple = " << multiple << endl;

		double_matrix ViViT = multiplyMatrix(Vi, transposeMatriz(Vi));
		printMatrix(transposeMatriz(Vi), "Vi_t");
		printMatrix(ViViT, "Vi x Vi_t");

		printMatrix(identity(N-i), "I(N-i)");
		printMatrix(numberMultiplyMatrix(ViViT, multiple), "Multiple x (Vi x Vi_t)");

		H = matrixMinusMatrix( identity(N-i), numberMultiplyMatrix(ViViT, multiple) );

		printMatrix(H, "Hi Temp");

		H = expandMatrix(H, i);

		printMatrix(H, "Hi Complete");

		A = multiplyMatrix(H, A);

		printMatrix(A, "New A");

		Q = multiplyMatrix(Q, H);

		i++;
	}

	//Creating R by A
	for(int i = 0; i < N; i++)
	{
		double_vector row;
		R.push_back(row);

		for(int j = 0; j < N; j++)
		{
			R[i].push_back(A[i][j]);
		}
	}

	printMatrix(Q, "Matrix Q");
	printMatrix(R, "Matriz R");

	double_matrix Q_t = transposeMatriz(Q);
	//Calculating X values
	for(int i = 0; i < M; i++)
	{
		cout << "CASE #" << (i+1) << ":" << endl << endl;
		printVector(B, i, "B" + to_string(i+1) + ": ");

		double_matrix B_Matrix;

		for (int j = 0; j < N; j++)
		{
			double_vector row;
			row.push_back(B[j][i]);
			B_Matrix.push_back(row);
		}

		double_matrix Y = multiplyMatrix(Q_t, B_Matrix);
		printVector(Y, 0, "Matrix Y (y = Q_t*B)");

		//y = Rx - Finding out X
		for(int l=N-1; l >= 0; l--)
		{
			double value = Y[l][0];

			for(int c=l+1; c < N; c++)
			{
				value -= R[l][c]*X[i][c];
			}

			X[i][l] = value/R[l][l];
		}

		printVector(X[i],"Result (Rx = y)");
	}
	
	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}

void printMatrix(double_matrix U, string title)
{
	cout << title << endl;

	for(int i = 0; i < U.size(); i++)
	{
		for(int j = 0; j < U[i].size(); j++)
		{
			cout << setw(15) << fixed << U[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
}

void printVector(double_vector U, string title)
{
	cout << title << endl;

	for(int j = 0; j < U.size(); j++)
	{
		cout << (j+1) << ": " << setw(15) << fixed << U[j] << endl;
	}

	cout << endl;
}

void printVector(double_matrix U, int L, string title)
{
	cout << title << endl;

	for(int j = 0; j < U.size(); j++)
	{
		cout << (j+1) << ": " << setw(15) << fixed << U[j][L] << endl;
	}

	cout << endl;
}


double_matrix multiplyMatrix(double_matrix A, double_matrix B)
{
	double_matrix C;

	for(int i = 0; i < A.size() ; i++)
	{
		double_vector row;
		C.push_back(row);

		for(int j = 0; j < B[0].size() ; j++)
		{
			double sum = 0;

			for(int k = 0; k < A[0].size() ; k++)
			{
				sum += A[i][k] * B[k][j];
			}

			C[i].push_back(sum);
		}
	}

	return C;
}

double_matrix transposeMatriz(double_matrix A)
{
	double_matrix B;

	int l = A[0].size();
	int c = A.size();

	for(int i = 0; i < l; i++)
	{
		double_vector row;
		B.push_back(row);

		for(int j = 0; j < c; j++)
		{
			B[i].push_back(A[j][i]);
		}
	}

	return B;
}

double innerProductSpace(double_vector v, double_vector u){
	double tot = 0;
	for(int i = 0; i < N; i++)
	{
		tot += v[i]*u[i];
	}
	return tot;
}

double norm(double_vector u){
	double sum = 0;

	for(int i = 0; i < u.size(); i++)
	{
		sum += pow(u[i], 2);
	}

	return sqrt(sum);
}

double signal(double n){
	if(n == 0)
	{
		return 1;
	}

	return n/abs(n);
}

double_matrix identity(int n){
	double_matrix I;

	for(int i = 0; i < n ; i++)
	{
		double_vector row;
		I.push_back(row);

		for(int j = 0; j < n ; j++)
		{
			if(i == j){
				I[i].push_back(1);
			}
			else
			{
				I[i].push_back(0);
			}
		}
	}

	return I;
}

double_matrix matrixMinusMatrix(double_matrix A, double_matrix B)
{
	double_matrix C;

	for(int i = 0; i < A.size(); i++)
	{
		double_vector row;
		C.push_back(row);

		for(int j = 0; j < A[i].size(); j++)
		{
			C[i].push_back(A[i][j] - B[i][j]);
		}
	}

	return C;
}

double_matrix numberMultiplyMatrix(double_matrix A, double M)
{
	double_matrix B;

	for(int i = 0; i < A.size(); i++)
	{
		double_vector row;
		B.push_back(row);

		for(int j = 0; j < A[i].size(); j++)
		{
			B[i].push_back(A[i][j] * M);
		}
	}

	return B;
}

double_matrix expandMatrix(double_matrix A, int I){
	double_matrix B;

	for(int i = 0; i < N; i++)
	{
		double_vector row;
		B.push_back(row);

		for(int j = 0; j < N; j++)
		{
			if((i < I) || (j < I))
			{
				if(j == i)
				{
					B[i].push_back(1);
				}
				else
				{
					B[i].push_back(0);
				}
			}
			else{
				B[i].push_back(A[i-I][j-I]);
			}
		}
	}

	return B;
}