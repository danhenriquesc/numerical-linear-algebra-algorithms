/*
	Algorithm to Orthogonalization of Vectors
	Gram-Schmidt Process
	Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

	Input

	N M (N = Number of vectors | M = Number os elements of each vectors)
	Vec1El1 Vec1El2 Vec1El3 Vec1El4 ... Vec1ElM
	Vec2El1 Vec2El2 Vec2El3 Vec2El4 ... Vec2ElM
	Vec3El1 Vec3El2 Vec3El3 Vec3El4 ... Vec3ElM
	  ...	  ...	  ...	  ...	...	  ...
	VecNEl1 VecNEl2 VecNEl3 VecNEl4 ... Vec3ElM
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <ctime>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

int N, M;
double innerProductSpace(double_vector, double_vector);
double_vector projection(double_vector, double_vector);
void printMatrix(double_matrix, string);
const clock_t begin_time = clock();


int main(){
	double_matrix V, U;
	double el;

	//reading input
	cin >> N >> M;

	for(int i = 0; i < N; i++){
		double_vector row;
		V.push_back(row);
		U.push_back(row);


		for(int j = 0; j < M; j++){
			cin >> el;
			V[i].push_back(el);
		}
	}

	printMatrix(V, "Original Matrix");

	//Gram-Schmidt Process
	U[0] = V[0];

	for(int i = 1; i < N; i++){
		for(int j = 0; j < M; j++){
			el = V[i][j];

			for(int t=0;t<i;t++){
				el -= projection(U[t],V[i])[j];
			}

			U[i].push_back(el);
		}
	}

	printMatrix(U, "Orthogonal Matrix");
	
	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}


double innerProductSpace(double_vector v, double_vector u){
	double tot = 0;
	for(int i = 0; i < M; i++){
		tot += v[i]*u[i];
	}
	return tot;
}

double_vector projection(double_vector u, double_vector v){
	double_vector p;
	double el;

	double tmp_div = innerProductSpace(v,u)/innerProductSpace(u,u);

	for(int i = 0; i<M; i++){
		el = tmp_div*u[i];
		p.push_back(el);
	}

	return p;
}

void printMatrix(double_matrix U, string title){
	cout << title << endl;

	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			cout << setw(15) << fixed << U[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
}