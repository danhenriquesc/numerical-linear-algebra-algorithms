/*
	Algorithm to solve Linear Equations Algorithms (Just for square matrix)
	Gauss Elimination
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
	double M[MAX_SIZE][MAX_SIZE];
	double AM[MAX_SIZE][MAX_SIZE];
	double B[MAX_SIZE][MAX_SIZE];
	double temp[MAX_SIZE];
	double mult;

	//n is the nxn size of square matrix
	//m is the numbers of B's in Ax = B
	cin >> n >> m;

	//reading input
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			cin >> M[i][j];
		}
		for(int j=1;j<=m;j++){
			cin >> B[j][i];	
		}
	}

	//proccessing all B's in Ax = B
	for(int j=1;j<=m;j++){
		cout << "CASE #" << j << ":" << endl << endl;

		//creating the Augmented matrix
		memcpy(AM,M,sizeof(M));
		for (int i = 1; i <= n; i++)
		{
			AM[i][n+1] = B[j][i];
		}

		printMatrix(AM, "Augmented matrix");

		//row reduction
		for(int c=1;c<=n;c++){
			int l = c;
			while(AM[l][c] == 0){
				l++;
			}
			if(l!=c){
				memcpy(temp,AM[c],sizeof(temp));
				memcpy(AM[c],AM[l],sizeof(AM[c]));
				memcpy(AM[l],temp,sizeof(AM[l]));
			}

			for(l=c+1;l<=n;l++){
				mult = AM[l][c]/AM[c][c];

				for(int j=1;j<=n+1;j++){
					AM[l][j] = AM[l][j]-mult*AM[c][j];
				}
			}
		}

		printMatrix(AM, "Row Echelon Form");
		
		//finding out variables
		for(int l=n;l>=1;l--){
			double value = AM[l][n+1];

			for(int c=l+1;c<=n;c++){
				value -= AM[l][c]*B[j][c];
			}

			B[j][l] = value/AM[l][l];
		}

		printMatrix(B[j],"Result","x");
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
		cout << "  |  " << setw(12) << A[i][n+1] << endl;
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