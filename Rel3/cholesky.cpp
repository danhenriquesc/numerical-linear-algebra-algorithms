/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Cholesky
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
const clock_t begin_time = clock();

int main(){
	double_matrix A, B, X, Y, G, G_t;
	double el;

	//N is the NxN size of square matrix
	cin >> N >> M;

	//reading input
	for(int i=0; i < N; i++)
	{
		double_vector row;
		A.push_back(row);
		B.push_back(row);
		X.push_back(row);
		Y.push_back(row);
		G.push_back(row);
		G_t.push_back(row);


		for(int j=0; j < N; j++)
		{
			cin >> el;
			A[i].push_back(el);

			G[i].push_back(0);
			G_t[i].push_back(0);
		}

		for(int j=0; j < M; j++)
		{
			cin >> el;
			B[i].push_back(el);
			X[i].push_back(0);
			Y[i].push_back(0);
		}
	}

	//cholesky

	// j = column | i = line
	double sum;
	for(int j = 0; j < N; j++)
	{
		for(int i = j; i < N; i++)
		{
			//diagonal
			if(j == i)
			{
				sum = 0;
				
				for(int k = 0; k < i; k++)
				{
					sum += pow(G[i][k], 2);
				}

				G[i][i] = sqrt(A[i][i] - sum);
				G_t[i][i] = G[i][i];
			}
			else
			{
				sum = 0;
				
				for(int k = 0; k < j; k++)
				{
					sum += G[i][k] * G[j][k];
				}

				G[i][j] = (A[i][j] - sum)/G[j][j];
				G_t[j][i] = G[i][j];
			}
		}
	}

	printMatrix(A, "Matrix A");
	printMatrix(G, "Matrix G");
	printMatrix(G_t, "Matrix Gt");

	//Calculating X values
	for(int i = 0; i < M; i++)
	{
		cout << "CASE #" << (i+1) << ":" << endl << endl;
		printVector(B[i], "B" + to_string(i+1) + ": ");

		//Gy = B - Finding out y
		for(int l=0; l < N; l++)
		{
			double value = B[l][i];

			for(int c=0; c < l; c++)
			{
				value -= G[l][c]*Y[i][c];
			}

			Y[i][l] = value/G[l][l];
		}

		printVector(Y[i],"Matrix Y (Gy = B)");

		//y = G_t*x - Finding out X
		for(int l=N-1; l >= 0; l--)
		{
			double value = Y[i][l];

			for(int c=l+1; c < N; c++)
			{
				value -= G_t[l][c]*X[i][c];
			}

			X[i][l] = value/G_t[l][l];
		}

		printVector(X[i],"Result (Ux = y)");
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