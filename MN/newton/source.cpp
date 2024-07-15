
#include <iostream>  
#include <math.h>  
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

typedef void (*FuncPointer)(const double* x, double* y, double* Df);

void printVector(const double* x, unsigned N){
  for(unsigned i=0;i<N;++i)
    printf("%17.17f ",x[i]);
  printf("\n");
}

int findCurve(FuncPointer f, double* x, unsigned k, double h) {
  double y[2]; // values of f at xs
  double Df[6]; // derrivativs 

  for(int i = 1; i <= k; i++) {
    // step of parameter
    x[2] += h;

    for(int j = 0; j < 1000 ; j++) {
      // to many steps
      if(j + 1 == 1000) {
        return i;
      }

      // function with given parameters returning values and derivatives
      f(x, y, Df);

      // checking accuracy
      if(max(abs(y[0]), abs(y[1])) <= 1.0e-14) {
        printVector(x, 3);
        break;
      }
      
      /*
      Jac = ( a b )     Jac ^-1 * 1 / detJac = ( d -b )
            ( c d )                            ( -c a )

      det = a * d - c * d 
      multiplaying by vector of ys 
      accidentaly changed here using ctrl f and change all xd
      */
      x[0] -= (y[0] * Df[4] - y[1] * Df[1]) / (Df[0] * Df[4] - Df[1] * Df[3]);
      x[1] -= (y[1] * Df[0] - y[0] * Df[3]) / (Df[0] * Df[4] - Df[1] * Df[3]);
    }
  }

  return 0;
}

int findSurface(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2) {
  double y[1];
  double Df[3];

  double x2 = x[2];

  for(int i = 1; i <= k1; i++) {
    x[1] += h1;
    x[2] = x2;
    for(int j = 1; j <= k2; j++) {
      x[2] += h2;

      for(int k = 0; k < 1000; k++) {
        if(k + 1 == 1000) {
          return i * k1 + j;
        }

        f(x, y, Df);

        if(abs(y[0]) <= 1.0e-14) {
          printVector(x, 3);
          break;
        }
        // only first derivative needed same logic
        x[0] -= y[0] / Df[0];
      }
    }
    if(i != k1)
     cout << endl;
  }
  return 0;
}

int findFixedPoints(FuncPointer f, double* x, unsigned k1, unsigned k2, double h1, double h2) {
  double y[2];
  double Df[8];

  double x3 = x[3];

  for(int i = 1; i <= k1; i++) {
    x[2] += h1;
    x[3] = x3;
    for(int j = 1; j <= k2; j++) {
      x[3] += h2;

      for(int k = 0; k < 1000; k++) {
        if(k + 1 == 1000) {
          return i * k1 + j;
        }

        f(x, y, Df);
        if(max(abs(y[0] - x[0]), abs(y[1] - x[1])) <= 1.0e-14) {
          printVector(x, 4);
          break;
        }

        // f(x,y,a,b) = (f1(x,y,a,b), f2(x,y,a,b))
        // g(x,y,a,b) = (f1(x,y,a,b) - x, f2(x,y,a,b) - y)
        // when f(x,y,a,b) = (x,y) then g(x,y,a,b) = (0,0)
        // g is just slight modification of f so need to subtract 1 from df[0] df[5]
        Df[0] -= 1;
        Df[5] -= 1;
        // dont need to multiply by df[1] always 1
        // forgot that dont need to always be 1, might be checked on different function 
        x[0] -= ((y[0] - x[0]) * Df[5] - ((y[1] - x[1]) * Df[1])) / (Df[0] * Df[5] - Df[4] * Df[1]);
        x[1] -= ((y[1] - x[1]) * Df[0] - ((y[0] - x[0]) * Df[4])) / (Df[0] * Df[5] - Df[4] * Df[1]);  
      }
    }
    if(i != k1)
     cout << endl;
  }
  return 0;
}

/*
  void henon(const double* x, double* y, double* Df){
    // funkcja dana jest wzorem henon(x,y,a,b) = (1+y-a*x^2,b*x)
    const double x2 = x[0]*x[0];
    
    y[0] = 1 + x[1] - x[2]*x2;
    y[1] = x[3]*x[0];
    
    //obliczam pierwszy wiersz macierzy
    Df[0] = -2*x[2]*x[0];
    Df[1] = 1.;
    Df[2] = -x2;
    Df[3] = 0.;
    
    //obliczam drugi wiersz macierzy
    Df[4] = x[3];
    Df[5] = 0.;
    Df[6] = 0.;
    Df[7] = x[0];
  }

  int main(){
    double x[4] = {-1.2807764064044151, -0.6403882032022076, 1.0000000000000000, 0.50000000000000000};
    findFixedPoints(henon,x,4,4,1./16,1./16);
    return 0;
  }

/**
spodziewane wyjscie:
-1.19763031176984502 -0.67366705037053787 1.06250000000000000 0.56250000000000000 
-1.16253262436707128 -0.72658289022941946 1.06250000000000000 0.62500000000000000 
-1.12828395985954577 -0.77569522240343769 1.06250000000000000 0.68750000000000000 
-1.09489692504918534 -0.82117269378688895 1.06250000000000000 0.75000000000000000 

-1.15709574728685860 -0.65086635784885805 1.12500000000000000 0.56250000000000000 
-1.12409377442300484 -0.70255860901437794 1.12500000000000000 0.62500000000000000 
-1.09187315546713570 -0.75066279438365580 1.12500000000000000 0.68750000000000000 
-1.06044486059083676 -0.79533364544312757 1.12500000000000000 0.75000000000000000 

-1.12017996019888111 -0.63010122761187060 1.18750000000000000 0.56250000000000000 
-1.08904242173442811 -0.68065151358401754 1.18750000000000000 0.62500000000000000 
-1.05862710283202821 -0.72780613319701948 1.18750000000000000 0.68750000000000000 
-1.02894361972548642 -0.77170771479411493 1.18750000000000000 0.75000000000000000 

-1.08638630667790914 -0.61109229750632399 1.25000000000000000 0.56250000000000000 
-1.05691785736085264 -0.66057366085053293 1.25000000000000000 0.62500000000000000 
-1.02811959340942205 -0.70683222046897776 1.25000000000000000 0.68750000000000000 
-1.00000000000000000 -0.75000000000000000 1.25000000000000000 0.75000000000000000 
*/

/*
#include "source.cpp"

  void henon(const double* x, double* y, double* Df){
    // funkcja dana jest wzorem henon(x,y,a,b) = (1+y-a*x^2,b*x)
    const double x2 = x[0]*x[0];
    
    y[0] = 1 + x[1] - x[2]*x2;
    y[1] = x[3]*x[0];
    
    //obliczam pierwszy wiersz macierzy
    Df[0] = -2*x[2]*x[0];
    Df[1] = 1.;
    Df[2] = -x2;
    Df[3] = 0.;
    
    //obliczam drugi wiersz macierzy
    Df[4] = x[3];
    Df[5] = 0.;
    Df[6] = 0.;
    Df[7] = x[0];
  }

  int main(){
    double x[4] = {-1.2807764064044151, -0.6403882032022076, 1.0000000000000000, 0.50000000000000000};
    findFixedPoints(henon,x,4,4,1./16,1./16);
    return 0;
  }

/**
spodziewane wyjscie:
-1.19763031176984502 -0.67366705037053787 1.06250000000000000 0.56250000000000000 
-1.16253262436707128 -0.72658289022941946 1.06250000000000000 0.62500000000000000 
-1.12828395985954577 -0.77569522240343769 1.06250000000000000 0.68750000000000000 
-1.09489692504918534 -0.82117269378688895 1.06250000000000000 0.75000000000000000 

-1.15709574728685860 -0.65086635784885805 1.12500000000000000 0.56250000000000000 
-1.12409377442300484 -0.70255860901437794 1.12500000000000000 0.62500000000000000 
-1.09187315546713570 -0.75066279438365580 1.12500000000000000 0.68750000000000000 
-1.06044486059083676 -0.79533364544312757 1.12500000000000000 0.75000000000000000 

-1.12017996019888111 -0.63010122761187060 1.18750000000000000 0.56250000000000000 
-1.08904242173442811 -0.68065151358401754 1.18750000000000000 0.62500000000000000 
-1.05862710283202821 -0.72780613319701948 1.18750000000000000 0.68750000000000000 
-1.02894361972548642 -0.77170771479411493 1.18750000000000000 0.75000000000000000 

-1.08638630667790914 -0.61109229750632399 1.25000000000000000 0.56250000000000000 
-1.05691785736085264 -0.66057366085053293 1.25000000000000000 0.62500000000000000 
-1.02811959340942205 -0.70683222046897776 1.25000000000000000 0.68750000000000000 
-1.00000000000000000 -0.75000000000000000 1.25000000000000000 0.75000000000000000 
*/
