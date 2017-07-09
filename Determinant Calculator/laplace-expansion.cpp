//ENTRADA
// N (Ordem de uma matriz quadrada)
// el1,1 el1,2 el1,3 .. el1,N
// el2,1 el2,2 el2,3 .. el2,N
// el3,1 el3,2 el3,3 .. el3,N
//  ...	  ... 	... 	 ...             
// elN,1 elN,2 elN,3 .. elN,N

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
			
			//pula linha 1 (index 0)
			for(int k = 1; k < N; k++){
				double_vector row;
				
				N_V.push_back(row);

				for(int j = 0; j < N; j++){
					//se for coluna atual, remove
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

	//Leitura
	cin >> N;

	for(int i = 0; i < N; i++){
		double_vector row;
		V.push_back(row);

		for(int j = 0; j < N; j++){
			cin >> el;
			V[i].push_back(el);
		}
	}

	//Laplace
	det = laplace(V, N);

	//Impressão
	cout << fixed << det << endl;

	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;
}