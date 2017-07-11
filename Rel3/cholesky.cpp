#include <iostream>
#include <cstring>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

#define MAX_SIZE 305
#define DECIMAL_THRESHOLD 0.00000000000000001

int N;
void printMatrix(double A[MAX_SIZE][MAX_SIZE], string);
const clock_t begin_time = clock();

int main(){
	double A[MAX_SIZE][MAX_SIZE], G[MAX_SIZE][MAX_SIZE], G_t[MAX_SIZE][MAX_SIZE];
	
	memset(G, 0, sizeof(G));
	memset(G_t, 0, sizeof(G_t));

	//N is the NxN size of square matrix
	cin >> N;

	//reading input
	for(int i=0; i < N; i++){
		for(int j=0; j < N; j++){
			cin >> A[i][j];
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

	printMatrix(G, "Matrix G");
	printMatrix(G_t, "Matrix Gt");

	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}

void printMatrix(double A[MAX_SIZE][MAX_SIZE], string title){
	cout << title << endl;

	for(int i=0; i < N ; i++)
	{
		for(int j=0; j < N; j++)
		{
			if(abs(A[i][j]) < DECIMAL_THRESHOLD)
				cout << setw(12) << "0";
			else
				cout << setw(12) << A[i][j];
		}
		cout << endl;
	}

	cout << endl;
}