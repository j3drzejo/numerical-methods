
#include "vectalg.h"
using namespace std;

// swaping rows
void swapRows(Matrix& A, Vector& b, int row1, int row2, int n) {
  double temp;
  for(int i = 0; i < n; i++) {
    temp = A(row1, i);
    A(row1, i) = A(row2, i);
    A(row2, i) = temp;
  }

  temp = b[row1];
  b[row1] = b[row2];
  b[row2] = temp;
}

// looking for max of scaled 
int Scaling(Matrix& A, Vector& maxNorm, int kol, int startRow, int n) {
  double m = maxNorm[startRow];
  double currMax = A(startRow, kol) / maxNorm[startRow];
  int currMaxIndex = startRow;
  for(int i = startRow; i < n; i++) {
    if(maxNorm[i] != 0) {
      double check = (A(i, kol) / maxNorm[i]);
      if(check > currMax) {
        currMaxIndex = i;
        currMax = check;
      }
    }
  }

  return currMaxIndex;
}

Vector MatrixTimesVector(const Matrix& A, const Vector& b, int n) {
  Vector sol(n);
  for(int i = 0; i < n; i++) {
    // zeroing forgot
    sol[i] = 0;
    for(int j = 0; j < n; j++) {
      sol[i] += A(i, j) * b[j];
    }
  }
  return sol;
}

Vector VectorDiff(const Vector& A, const Vector& B, int n) {
  Vector sol(n);

  for(int i = 0; i < n; i++) {
    sol[i] = A[i] - B[i];
  }

  return sol;
}

Vector VectorSum(const Vector& A, const Vector& B, int n) {
  Vector sol(n);

  for(int i = 0; i < n; i++) {
    sol[i] = A[i] + B[i];
  }

  return sol;
}

Vector solveEquations( 
	const Matrix & Aa, 
	const Vector & bb, 
	double  epss) {
	
	Matrix A(Aa);
  Vector b(bb);  

	int n = A.size(); // size of n 
  Vector x(n);

	// zeroing solution vector
	for(int i = 0; i < n; i++) {
		x[i] = 0;
	}
    
  Vector maxNorm(n);

  // caluclating forms for scalings
  for(int i = 0; i < n; i++) {
    double norm = 0;
    for(int j = 0; j < n; j++) {
      norm = max(norm, abs(A(i,j)));
    }
    maxNorm[i] = norm;
  }
  for(int k = 0; k < n - 1; k++) {
    // loking for scaling where is max of ... using logic from lecture
    int currMax = Scaling(A, maxNorm, k, k, n);
    // swapping max row with curr row
    swapRows(A, b, k, currMax, n);
    for(int i = k + 1; i < n; i++) {
      // fidning factor for zeroing
      double factor = A(i, k) / A(k, k);
      for(int j = k; j < n; j++) {
        A(i, j) -= factor * A(k, j);

      }
      b[i] -= factor * b[k];
    }
  }

  // backfactoring solution
  for(int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for(int j = i + 1; j < n; j++) {
      sum += A(i, j) * x[j];
    }
    x[i] = (b[i] - sum) / A(i,i);
  }
  b = bb;
  // getting U matrix;
  A = Aa;

  Vector Rk = VectorDiff(b, MatrixTimesVector(A, x, n), n);
  while(abs(Rk.max_norm()) > epss) {
    Vector Ek = solveEquations(A, Rk, epss);
    x = VectorSum(x, Ek, n);
    Rk = VectorDiff(b, MatrixTimesVector(A, x, n), n);
  }
	return x;
}
/*
int main(){
    cout.precision(17);
    int n = 0;
    double eps = 0;

    // wczytywanie danych
    cin >> n;
    Matrix a(n);
    Vector b(n);
    cin >> a >> b >> eps;

    Vector x = solveEquations(a, b, eps);

    auto residual = residual_vector(a, b, x);
    cout << "rozwiazanie = " << x << endl;
    cout << "rezydualny = " << residual << endl;
    cout << "blad = " << residual.max_norm()
         << " limit = " << eps << endl ;
    cout << "Test " << (residual.max_norm() < eps ? "":"nie ") << "zaliczony" << endl; 
    return 0;
}	
*/
/*
4
2 3 0 0 
6 3 9 0 
0 2 5 2
0 0 4 3
21 69 34 22
1e-15
*/