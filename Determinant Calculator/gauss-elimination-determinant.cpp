/*
	Algorithm to Calculate Matrix Determinant
	Gauss Elimination
	Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

	Input

	N (Square Matrix Order)
	Elem11 Elem12 Elem13 ... Elem1N
	Elem21 Elem22 Elem23 ... Elem2N
	Elem31 Elem32 Elem33 ... Elem3N
	  ...	...	    ...  ...   ...
	ElemN1 ElemN2 ElemN3 ... ElemNN

	Example:
	3 3
	1 2 3
	4 5 6
	7 8 9
*/

#include <iostream>
#include <cstring>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

#define MAX_SIZE 305
#define DECIMAL_THRESHOLD 0.0000000001

int n;
void printMatrix(double A[MAX_SIZE][MAX_SIZE], string);
const clock_t begin_time = clock();

int main(){
	double M[MAX_SIZE][MAX_SIZE];
	double temp[MAX_SIZE];
	double mult;
	double result;
	int signal;

	//n is the nxn size of square matrix
	cin >> n;

	//reading input
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			cin >> M[i][j];
		}
	}

	//row reduction
	signal = 1;
	for(int c=1;c<=n;c++){
		int l = c;
		while(M[l][c] == 0){
			l++;
		}
		if(l!=c){
			memcpy(temp,M[c],sizeof(temp));
			memcpy(M[c],M[l],sizeof(M[c]));
			memcpy(M[l],temp,sizeof(M[l]));
			signal *= -1;
		}

		for(l=c+1;l<=n;l++){
			mult = M[l][c]/M[c][c];

			for(int k=1;k<=n;k++){
				M[l][k] = M[l][k]-mult*M[c][k];
			}
		}
	}

	printMatrix(M, "Row Echelon Form");

	//calculing result
	result = signal;
	for(int i=1;i<=n;i++){
		result *= M[i][i];
	}
	cout << "Result: " << fixed << result << endl;

	cout << "Total Time: " << fixed << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " s" << endl;

	return 0;
}

void printMatrix(double A[MAX_SIZE][MAX_SIZE], string title){
	cout << title << endl;

	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(abs(A[i][j]) < DECIMAL_THRESHOLD)
				cout << setw(12) << "0";
			else
				cout << setw(12) << A[i][j];
		}
		cout << "  |  " << setw(12) << A[i][n+1] << endl;
	}

	cout << endl;
}