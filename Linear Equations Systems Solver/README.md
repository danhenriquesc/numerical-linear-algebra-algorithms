# Linear Equations Algorithms

Algorithm to solve Linear Equations Algorithms (Gauss Elimination, Gauss-Jordan Elimination, LU Decomposition, LDU Decomposition, Crout Algorithm) [Just for square matrix]

Implemented by: Daniel Henrique (daniel.henrique.sc@gmail.com | daniel.henrique@ime.uerj.br) - 2017

## Solve Ax = B(i..n)

#### Input:
```
N M
a11 a12 a13 ... a1n b11 b12 b13 ... b1m
a21 a22 a23 ... a2n b21 b22 b23 ... b2m
.   .   .       .   .   .   .       .
an1 an2 an3 ... ann bn1 bn2 bn3 ... bnm 
```

#### Example:
```
3 2
3 7 12 9 5
1 5 7 12 1
3 3 3 7 1
```

#### Solves:

|  3x<sub>1</sub> + 7x<sub>2</sub> + 12x<sub>3</sub> =  9<br>
|   x<sub>1</sub> + 5x<sub>2</sub> +  7x<sub>3</sub> = 12<br>
|  3x<sub>1</sub> + 3x<sub>2</sub> +  3x<sub>3</sub> =  7<br>
<br><br>
|  3x<sub>1</sub> + 7x<sub>2</sub> + 12x<sub>3</sub> =  5<br>
|   x<sub>1</sub> + 5x<sub>2</sub> +  7x<sub>3</sub> =  1<br>
|  3x<sub>1</sub> + 3x<sub>2</sub> +  3x<sub>3</sub> =  1<br>


### Gauss Elimination OUTPUT of Example INPUT:
```
CASE #1:

Augmented matrix
           3           7          12  |             9
           1           5           7  |            12
           3           3           3  |             7

Row Echelon Form
           3           7          12  |             9
           0     2.66667           3  |             9
           0           0        -4.5  |          11.5

Result
x1:   -1.36111
x2:       6.25
x3:   -2.55556

CASE #2:

Augmented matrix
           3           7          12  |             5
           1           5           7  |             1
           3           3           3  |             1

Row Echelon Form
           3           7          12  |             5
           0     2.66667           3  |     -0.666667
           0           0        -4.5  |            -5

Result
x1:   0.722222
x2:       -1.5
x3:    1.11111

Total Time: 0.000633 s
```

### Gauss-Jordan Elimination OUTPUT of Example INPUT:
```
CASE #1:

Augmented matrix
           3           7          12  |             9
           1           5           7  |            12
           3           3           3  |             7

Row Echelon Form
           3           7          12  |             9
           0     2.66667           3  |             9
           0           0        -4.5  |          11.5

Matrix with Diagonal 1
           1     2.33333           4  |             3
           0           1       1.125  |         3.375
           0           0           1  |      -2.55556

Reduced Row Echelon Form
           1           0           0  |      -1.36111
           0           1           0  |          6.25
           0           0           1  |      -2.55556

Result
x1:   -1.36111
x2:       6.25
x3:   -2.55556

CASE #2:

Augmented matrix
           3           7          12  |             5
           1           5           7  |             1
           3           3           3  |             1

Row Echelon Form
           3           7          12  |             5
           0     2.66667           3  |     -0.666667
           0           0        -4.5  |            -5

Matrix with Diagonal 1
           1     2.33333           4  |       1.66667
           0           1       1.125  |         -0.25
           0           0           1  |       1.11111

Reduced Row Echelon Form
           1           0           0  |      0.722222
           0           1           0  |          -1.5
           0           0           1  |       1.11111

Result
x1:   0.722222
x2:       -1.5
x3:    1.11111

Total Time: 0.000453 s
```

### LU Decomposition OUTPUT of Example INPUT:
```
Matrix A (Ax = B)
           3           7          12
           1           5           7
           3           3           3

Matrix Lower (LUx = B)
           1           0           0
    0.333333           1           0
           1        -1.5           1

Matrix Upper (LUx = B)
           3           7          12
           0     2.66667           3
           0           0        -4.5

CASE #1:

Matrix B
b1:          9
b2:         12
b3:          7

Matrix Y
y1:          9
y2:          9
y3:       11.5

Result
x1:   -1.36111
x2:       6.25
x3:   -2.55556

CASE #2:

Matrix B
b1:          5
b2:          1
b3:          1

Matrix Y
y1:          5
y2:  -0.666667
y3:         -5

Result
x1:   0.722222
x2:       -1.5
x3:    1.11111

Total Time: 0.000290 s
```

### LDU Decomposition OUTPUT of Example INPUT:
```
Matrix A (Ax = B)
           3           7          12
           1           5           7
           3           3           3

Matrix Lower (LU'x = B)
           1           0           0
    0.333333           1           0
           1        -1.5           1

Matrix Upper' (LU'x = B)
           3           7          12
           0     2.66667           3
           0           0        -4.5

Matrix Lower (LDUx = B)
           1           0           0
    0.333333           1           0
           1        -1.5           1

Matrix Diagonal' (LDUx = B)
           3           0           0
           0     2.66667           0
           0           0        -4.5

Matrix Upper' (LDUx = B)
           1     2.33333           4
           0           1       1.125
           0           0           1

CASE #1:

Matrix B
b1:          9
b2:         12
b3:          7

Matrix Z (Lz = B)
z1:          9
z2:          9
z3:       11.5

Matrix Y (Dy = z)
z1:          3
z2:      3.375
z3:   -2.55556

Result (Ux = y)
x1:   -1.36111
x2:       6.25
x3:   -2.55556

CASE #2:

Matrix B
b1:          5
b2:          1
b3:          1

Matrix Z (Lz = B)
z1:          5
z2:  -0.666667
z3:         -5

Matrix Y (Dy = z)
z1:    1.66667
z2:      -0.25
z3:    1.11111

Result (Ux = y)
x1:   0.722222
x2:       -1.5
x3:    1.11111

TEMPO TOTAL: 0.000441 seg
```

### Crout Algorithm OUTPUT of Example INPUT:
```
Matrix A (Ax = B)
           3           7          12
           1           5           7
           3           3           3

Matrix Lower (LUx = B)
           3           0           0
           1     2.66667           0
           3          -4        -4.5

Matrix Upper (LUx = B)
           1     2.33333           4
           0           1       1.125
           0           0           1

CASE #1:

Matrix B
b1:          9
b2:         12
b3:          7

Matrix Y (Ly = B)
y1:          3
y2:      3.375
y3:   -2.55556

Result (Ux = y)
x1:   -1.36111
x2:       6.25
x3:   -2.55556

CASE #2:

Matrix B
b1:          5
b2:          1
b3:          1

Matrix Y (Ly = B)
y1:    1.66667
y2:      -0.25
y3:    1.11111

Result (Ux = y)
x1:   0.722222
x2:       -1.5
x3:    1.11111

Total Time: 0.000232 s
```
