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
double_matrix multiplyMatrix(double_matrix, double_matrix);
double_matrix transposeMatriz(double_matrix);
const clock_t begin_time = clock();

int main(){
	double_matrix A, G, R, Q;
	double el;
	
	// N lines and M columns matrix
	cin >> N >> M;

	//reading input
	for(int i=0; i < N; i++)
	{
		double_vector row;
		A.push_back(row);
		G.push_back(row);
		R.push_back(row);
		Q.push_back(row);

		for(int j=0; j < M; j++)
		{
			cin >> el;
			A[i].push_back(el);
		}

		for(int j=0; j < N; j++)
		{
			G[i].push_back(0);
			Q[i].push_back(0);
		}
	}

	printMatrix(A, "A");

	// Givens Rotation
	for(int j = 0; j < M; j++)
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
	for(int i = 0; i < M; i++)
	{
		double_vector row;
		R.push_back(row);

		for(int j = 0; j < M; j++)
		{
			R[i].push_back(A[i][j]);
		}
	}

	printMatrix(Q, "Matrix Q");
	printMatrix(R, "Matriz R");


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