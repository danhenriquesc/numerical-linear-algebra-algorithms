/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Given's Rotation QR Decomposition
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
const clock_t begin_time = clock();

int main(){
	double_matrix A, G, R, Q, B, X;
	double el;
	
	// NxN matriz with M results
	cin >> N >> M;

	//reading input
	for(int i=0; i < N; i++)
	{
		double_vector row;
		A.push_back(row);
		G.push_back(row);
		R.push_back(row);
		Q.push_back(row);
		B.push_back(row);
		X.push_back(row);

		for(int j=0; j < N; j++)
		{
			cin >> el;
			A[i].push_back(el);
			G[i].push_back(0);
			Q[i].push_back(0);
		}

		for(int j=0; j < M; j++)
		{
			cin >> el;
			B[i].push_back(el);
			X[i].push_back(0);
		}
	}

	printMatrix(A, "A");

	// Givens Rotation
	for(int j = 0; j < N; j++)
	{
		for(int i = N-1; i > j; i--)
		{
			double r;
			r = sqrt( pow(A[i][j], 2) + pow(A[i-1][j], 2) );

			for(int n = 0; n < N; n++)
			{
				for(int m = 0; m < N; m++)
				{
					if( (n == i && m == i) || (n == (i-1) && m == (i-1)) )
					{
						G[n][m] = A[i-1][j]/r; //cos
					}
					else if( (n == i && m == (i-1)) || (n == (i-1) && m == i) )
					{
						G[n][m] = A[i][j]/r; //sin

						if(n > m)
						{
							G[n][m] *= -1; // -sin
						}
					}
					else if( n == m)
					{
						G[n][m] = 1;
					}
					else
					{
						G[n][m] = 0;
					}
				}
			}

			//first G
			if( (j == 0) && (i == N-1 ))
			{
				Q = transposeMatriz(G);
			}
			else
			{
				Q = multiplyMatrix(Q, transposeMatriz(G));
			}

			A = multiplyMatrix(G, A);
		}
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