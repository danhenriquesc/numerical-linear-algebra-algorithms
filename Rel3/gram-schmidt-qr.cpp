#include <iostream>
#include <vector>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

int N;
double innerProductSpace(double_vector, double_vector);
double_vector projection(double_vector, double_vector);
double norm(double_vector);
void printMatrix(double_matrix, string);
const clock_t begin_time = clock();


int main(){
	double_matrix V_t, V, U, Q, R, Ql, Rl, D, D_t;
	double el;

	//reading input
	cin >> N;

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