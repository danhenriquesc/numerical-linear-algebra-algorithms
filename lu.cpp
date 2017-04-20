/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	LU Decomposition
	Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

	Solve Ax = B(i..n)

	Input:
	N M
	a11 a12 a13 ... a1n b11 b12 b13 ... b1m
	a21 a22 a23 ... a2n b21 b22 b23 ... b2m
	 .   .   .       .   .   .   .       .
	an1 an2 an3 ... ann bn1 bn2 bn3 ... bnm 

	Example:
	3 3
	1 2 3 10 11 12
	4 5 6 13 14 15
	7 8 9 16 17 18

	Solves:
	1x+2y+3z = 10
	4x+5y+6z = 13
	7x+8y+9z = 16

	1x+2y+3z = 11
	4x+5y+6z = 14
	7x+8y+9z = 15

	1x+2y+3z = 12
	4x+5y+6z = 15
	7x+8y+9z = 18
*/

#include <iostream>
#include <cstring>
#include <iomanip>
#include <ctime>
#include <cmath>

using namespace std;

#define MAX_SIZE 105
#define DECIMAL_THRESHOLD 0.0000000001

int n, m;
void printMatrix(double A[MAX_SIZE][MAX_SIZE],string);
void printMatrix(double B[MAX_SIZE],string,string);
const clock_t begin_time = clock();

int main(){
	double A[MAX_SIZE][MAX_SIZE];
	double L[MAX_SIZE][MAX_SIZE];
	double U[MAX_SIZE][MAX_SIZE];

	double B[MAX_SIZE][MAX_SIZE];
	double temp[MAX_SIZE];
	double X[MAX_SIZE][MAX_SIZE];
	double Y[MAX_SIZE][MAX_SIZE];
	double mult;

	//n is the nxn size of square matrix
	//m is the numbers of B's in Ax = B
	cin >> n >> m;

	//reading input
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			cin >> A[i][j];
		}
		for(int j=1;j<=m;j++){
			cin >> B[j][i];	
		}
	}

	printMatrix(A, "Matrix A (Ax = B)");
	
	memset(L,0,sizeof(L));
	for(int i=1;i<=n;i++) L[i][i] = 1;

	//LU decomposition
	for(int c=1;c<=n;c++){
		int l = c;
		while(A[l][c] == 0){
			l++;
		}
		if(l!=c){
			memcpy(temp,A[c],sizeof(temp));
			memcpy(A[c],A[l],sizeof(A[c]));
			memcpy(A[l],temp,sizeof(A[l]));
		}

		for(l=c+1;l<=n;l++){
			mult = A[l][c]/A[c][c];
			L[l][c] = mult;

			for(int j=1;j<=n+1;j++){
				A[l][j] = A[l][j]-mult*A[c][j];
			}
		}
	}

	memcpy(U,A,sizeof(U));
	printMatrix(L, "Matrix Lower (LUx = B)");
	printMatrix(U, "Matrix Upper (LUx = B)");


	//proccessing all B's in LUx = B
	for(int j=1;j<=m;j++){
		cout << "CASE #" << j << ":" << endl << endl;

		printMatrix(B[j],"Matrix B","b");
		
		//Ly = B - Finding out y
		for(int l=1;l<=n;l++){
			double value = B[j][l];

			for(int c=1;c<l;c++){
				value -= L[l][c]*Y[j][c];
			}

			Y[j][l] = value/L[l][l];
		}

		printMatrix(Y[j],"Matrix Y (Ly = B)","y");

		//y = Ux - Finding out X
		for(int l=n;l>=1;l--){
			double value = Y[j][l];

			for(int c=l+1;c<=n;c++){
				value -= U[l][c]*X[j][c];
			}

			X[j][l] = value/U[l][l];
		}

		printMatrix(X[j],"Result (Ux = y)","x");
	}

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
		cout << endl;
	}

	cout << endl;
}

void printMatrix(double B[MAX_SIZE], string title, string var){
	cout << title << endl;

	for(int i=1;i<=n;i++){
		cout << var << i << ": " << setw(10) << B[i] << endl;
	}

	cout << endl;
}