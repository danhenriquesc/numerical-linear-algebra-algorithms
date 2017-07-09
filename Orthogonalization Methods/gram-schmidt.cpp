//ENTRADA
// N (numero de vet.) M (numero de elementos de cada vetor)
// Vetor1Elem1 Vetor1Elem2 Vetor1Elem3 Vetor1Elem4 ... Vetor1ElemS
// Vetor2Elem1 Vetor2Elem2 Vetor2Elem3 Vetor2Elem4 ... Vetor2ElemS
// Vetor3Elem1 Vetor3Elem2 Vetor3Elem3 Vetor3Elem4 ... Vetor3ElemS
// Vetor4Elem1 Vetor4Elem2 Vetor4Elem3 Vetor4Elem4 ... Vetor4ElemS
//     ...		   ... 		   ... 		   ...             ...
// VetorNElem1 VetorNElem2 VetorNElem3 VetorNElem4 ... VetorNElemS

#include <iostream>
#include <vector>
#include <ctime>

using namespace std;

typedef vector<double> double_vector;
typedef vector<double_vector> double_matrix;

int N, M;
const clock_t begin_time = clock();

// Calculando produto interno
double innerProductSpace(double_vector v, double_vector u){
	double tot = 0;
	for(int i = 0; i < M; i++){
		tot += v[i]*u[i];
	}
	return tot;
}

// Calculando projeção
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

int main(){
	double_matrix V, U;
	double el;

	//Leitura
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

	//Realizando o Processo de Gram-Schmidt
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

	//Imprimindo
	for(int i = 0; i < N; i++){
		for(int j = 0; j < M; j++){
			cout << U[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}
