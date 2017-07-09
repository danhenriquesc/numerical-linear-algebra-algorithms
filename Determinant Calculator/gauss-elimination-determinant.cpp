/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Determinant of Matrix using Gauss Elimination
	Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

	Solve Det(A)

	Input:
	N M
	a11 a12 a13 ... a1n
	a21 a22 a23 ... a2n
	 .   .   .       .  
	an1 an2 an3 ... ann

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
void printMatrix(double A[MAX_SIZE][MAX_SIZE],string);
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