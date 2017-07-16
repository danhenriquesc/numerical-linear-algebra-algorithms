/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Gram-Schmidt QR Decomposition
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
double innerProductSpace(double_vector, double_vector);
double_vector projection(double_vector, double_vector);
double norm(double_vector);
double_matrix multiplyMatrix(double_matrix, double_matrix);
double_matrix transposeMatriz(double_matrix);
void printMatrix(double_matrix, string);
void printVector(double_vector, string);
void printVector(double_matrix, int, string);
const clock_t begin_time = clock();


int main(){
	double_matrix V_t, V, U, Q, R, Ql, Rl, D, D_t, B, X;
	double el;

	//reading input
	cin >> N >> M;

	for(int i = 0; i < N; i++)
	{
		double_vector row;
		V_t.push_back(row);
		V.push_back(row);
		U.push_back(row);
		Q.push_back(row);
		R.push_back(row);
		Ql.push_back(row);
		Rl.push_back(row);
		D.push_back(row);
		D_t.push_back(row);
		B.push_back(row);
		X.push_back(row);

		for(int j = 0; j < N; j++)
		{
			cin >> el;
			V_t[i].push_back(el);
		}

		for(int j = 0; j < N; j++)
		{
			Q[i].push_back(0);
			Ql[i].push_back(0);
			R[i].push_back(0);
			Rl[i].push_back(0);
			D[i].push_back(0);
			D_t[i].push_back(0);
			V[i].push_back(0);
		}

		for(int j=0; j < M; j++)
		{
			cin >> el;
			B[i].push_back(el);
			X[i].push_back(0);
		}
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			V[j][i] = V_t[i][j];
		}
	}
	
	printMatrix(V_t, "Original T Matrix");
	printMatrix(V, "Original Matrix");

	//Gram-Schmidt Process
	U[0] = V[0];

	for(int i = 1; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			el = V[i][j];

			for(int t=0;t<i;t++)
			{
				el -= projection(U[t],V[i])[j];
			}

			U[i].push_back(el);
		}
	}

	printMatrix(U, "Orthogonal Matrix");

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			Ql[j][i] = U[i][j];
		}
	}

	printMatrix(Ql, "Matrix Q'");

	for(int i = 0; i < N; i++)
	{
		D[i][i] = 1/norm(U[i]);
	}

	printMatrix(D, "Matrix D");


	for(int n = 0; n < N; n++)
	{
		for(int m = 0; m < N; m++)
		{
			if(n > m)
			{
				Rl[n][m] = 0;
			}
			else
			{
				Rl[n][m] = innerProductSpace(V[m],U[n]);
			}
		}
	}

	printMatrix(Rl, "Matrix R'");

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			D_t[j][i] = D[i][j];
		}
	}

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			double sum = 0;

			for(int k = 0; k < N; k++)
			{
				sum += Ql[i][k] * D_t[k][j];
			}

			Q[i][j] = sum;
		}
	}

	printMatrix(Q, "Matrix Q");

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			double sum = 0;

			for(int k = 0; k < N; k++)
			{
				sum += D[i][k] * Rl[k][j];
			}

			R[i][j] = sum;
		}
	}
	
	printMatrix(R, "Matrix R");

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
				value -= R[l][c]*X[c][i];
			}

			X[l][i] = value/R[l][l];
		}

		printVector(X, i, "Result (Rx = y)");
	}
	
	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}


double innerProductSpace(double_vector v, double_vector u){
	double tot = 0;
	for(int i = 0; i < N; i++)
	{
		tot += v[i]*u[i];
	}
	return tot;
}

double_vector projection(double_vector u, double_vector v){
	double_vector p;
	double el;

	double tmp_div = innerProductSpace(v,u)/innerProductSpace(u,u);

	for(int i = 0; i<N; i++)
	{
		el = tmp_div*u[i];
		p.push_back(el);
	}

	return p;
}

double norm(double_vector u){
	double sum = 0;

	for(int i = 0; i < N; i++)
	{
		sum += pow(u[i], 2);
	}

	return sqrt(sum);
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

	for(int i = 0; i < A.size(); i++)
	{
		double_vector row;
		B.push_back(row);

		for(int j = 0; j < A[i].size(); j++)
		{
			B[i].push_back(A[j][i]);
		}
	}

	return B;
}