#include <iostream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

int N;
void printMatrix(double_matrix, string);
const clock_t begin_time = clock();

int main(){
	double_matrix A, G, G_t;
	double el;

	//N is the NxN size of square matrix
	cin >> N;

	//reading input
	for(int i=0; i < N; i++){
		double_vector row;
		A.push_back(row);
		G.push_back(row);
		G_t.push_back(row);


		for(int j=0; j < N; j++){
			cin >> el;
			A[i].push_back(el);

			G[i].push_back(0);
			G_t[i].push_back(0);
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

void printMatrix(double_matrix U, string title){
	cout << title << endl;

	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < N; j++)
		{
			cout << setw(15) << fixed << U[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
}