/*
	Algorithm to Calculate Matrix Determinant
	Laplace Expansion

	Input

	N (Square Matrix Order)
	Elem11 Elem12 Elem13 ... Elem1N
	Elem21 Elem22 Elem23 ... Elem2N
	Elem31 Elem32 Elem33 ... Elem3N
	  ...	...	    ...  ...   ...
	ElemN1 ElemN2 ElemN3 ... ElemNN
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

const clock_t begin_time = clock();

double laplace(double_matrix V, int N){
	double result = 0, el;

	if(N == 2){
		return (V[0][0]*V[1][1]) - (V[0][1]*V[1][0]);
	}else{
		for(int i=0; i<N;i++){
			double_matrix N_V;
			
			//ignores line 1 (index 0)
			for(int k = 1; k < N; k++){
				double_vector row;
				
				N_V.push_back(row);

				for(int j = 0; j < N; j++){
					//ignores current column
					if(i==j) continue;

					el = V[k][j];

					N_V[k-1].push_back(el);
				}
			}

			result += pow(-1,i) * V[0][i] * laplace(N_V, N-1);
		}

		return result;
	}
}

int main(){
	double_matrix V;
	double el, det;
	int N;

	//reading input
	cin >> N;

	for(int i = 0; i < N; i++){
		double_vector row;
		V.push_back(row);

		for(int j = 0; j < N; j++){
			cin >> el;
			V[i].push_back(el);
		}
	}

	//Laplace Expansion
	det = laplace(V, N);

	//printing result
	cout << "Determinant: " << fixed << det << endl;

	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;
}