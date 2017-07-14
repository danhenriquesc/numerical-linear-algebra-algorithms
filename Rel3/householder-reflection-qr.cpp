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
double_matrix identity(int);
double_matrix matrixMinusMatrix(double_matrix, double_matrix);
double_matrix numberMultiplyMatrix(double_matrix, double);
double_matrix expandMatrix(double_matrix, int);
double norm(double_vector);
double innerProductSpace(double_vector, double_vector);
double signal(double);
const clock_t begin_time = clock();

int main(){
	double_matrix A, H, Q, R;
	double el;

	// N lines and M columns matrix
	cin >> N >> M;

	for(int i = 0; i < N; i++)
	{
		double_vector row;
		A.push_back(row);
		
		for(int j = 0; j < M; j++)
		{
			cin >> el;
			A[i].push_back(el);
		}
	}

	printMatrix(A, "A");

	Q = identity(N);

	int i = 0;
	while(true)
	{	
		if( (N-i < 2) || (M-i < 2) )
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