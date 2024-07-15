
#include "source.h"
#include <iostream>

using namespace std;

spline::spline(int n) : n(n) {
  x = new double[n];
  y = new double[n];
  h = new double[n];
  a = new double[n];
  b = new double[n];
  c = new double[n];
  d = new double[n];
}

void spline::set_points(double x[], double y[]) {
  double u[n];
  double z[n];
  double v[n - 1];
  u[0] = 1;
  z[0] = 0;
  u[n - 1] = 1;
  z[n - 1] = 0;
  c[n - 1] = 0;

  for(int i = 0; i < n; i++) {
    this->x[i] = x[i];
    this->y[i] = y[i];
  }

  for(int i = 0; i < n - 1; i++) {
    h[i] = x[i + 1] - x[i];
}

  for(int i = 1; i < n - 1; i++) {
    if(i == 1) {
      u[i] = 2 * (x[i + 1] - x[i - 1]);
    } else {
      u[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * (h[i- 1] / u[i - 1]);
    }
    z[i] = (((3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1])) - h[i - 1] * z[i - 1]) / u[i];
  }

  for (int j = n - 2; j >= 0; j--) {
    if(j == 0) {
      c[j] = z[j];
    } else {
      c[j] = z[j] - (h[j] / u[j]) * c[j + 1];
    }
    b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
    d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    a[j] = y[j];
  }
}

double spline::operator()(double z) const {
  int i;
  for (i = 0; i < n - 1 && z > x[i + 1]; i++) {}
  double dx = z - x[i];
  return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}


spline::~spline() {
  delete[] x;
  delete[] y;
  delete[] h;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] d;
}
