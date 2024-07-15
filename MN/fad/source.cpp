
#include "funkcja.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

class Jet {
  public:
  long double ders[6];

  Jet() {
    for(int i = 0; i < 6; i++) {
      ders[i] = 0;
    }
  }

  Jet(int val) {
    ders[0] = val;
    for(int i = 1; i < 6; i++) {
      ders[i] = 0;
    }
  }

  Jet(double val) {
    ders[0] = val;
    for(int i = 1; i < 6; i++) {
      ders[i] = 0;
    }
  }

  Jet(double x, bool isX) {
    ders[0] = x;
    for(int i = 1; i < 6; i++) {
      ders[i] = 0;
    }

    if(isX) {
      ders[1] = 1;
    } else if(!isX) {
      ders[2] = 1;
    }
  }

  Jet(const Jet& other) {
    for (int i = 0; i < 6; ++i) {
      this->ders[i] = other.ders[i];
    }
  }

  Jet& operator=(const Jet& other) {
    if (this != &other) {
      for (int i = 0; i < 6; ++i) {
        this->ders[i] = other.ders[i];
      }
    }
    return *this;
  }

  friend Jet operator +(const Jet& f, const Jet& g) {
    Jet result;
    for(int i = 0; i < 6; i++) {
      result.ders[i] = f.ders[i] + g.ders[i];
    }
    return result;
  }
  
  friend Jet operator -(const Jet& f) {
    Jet result;
    for(int i = 0; i < 6; i++) {
      result.ders[i] = - f.ders[i];
    }
    return result;
  }

  friend Jet operator -(const Jet& f, const Jet& g) {
    Jet result = f + -(g);
    return result;
  }

  friend Jet operator *(const Jet& f, const Jet& g) {
    Jet result;

    result.ders[0] = f.ders[0] * g.ders[0];
    result.ders[1] = f.ders[1] * g.ders[0] + f.ders[0] * g.ders[1];
    result.ders[2] = f.ders[2] * g.ders[0] + f.ders[0] * g.ders[2];
    result.ders[3] = f.ders[3] * g.ders[0] + f.ders[1] * g.ders[1] + f.ders[1] * g.ders[1] + f.ders[0] * g.ders[3];
    result.ders[4] = f.ders[4] * g.ders[0] + f.ders[1] * g.ders[2] + f.ders[2] * g.ders[1] + f.ders[0] * g.ders[4];
    result.ders[5] = f.ders[0] * g.ders[5] + f.ders[2] * g.ders[2] + f.ders[2] * g.ders[2] + f.ders[5] * g.ders[0];

    return result;
  }


  friend Jet operator /(const Jet& f, const Jet& g) {
    Jet result;

    result.ders[0] = f.ders[0] / g.ders[0];
    result.ders[1] = (f.ders[1] - result.ders[0] * g.ders[1]) / g.ders[0];
    result.ders[2] = (f.ders[2] - result.ders[0] * g.ders[2]) / g.ders[0];
    result.ders[3] = (f.ders[3] - 2 * result.ders[1] * g.ders[1] - result.ders[0] * g.ders[3]) / g.ders[0];
    result.ders[4] = (f.ders[4] - result.ders[2] * g.ders[1] - result.ders[0] * g.ders[4] - result.ders[1] * g.ders[2]) / g.ders[0];
    result.ders[5] = (f.ders[5] - 2 * result.ders[2] * g.ders[2] - result.ders[0] * g.ders[5]) / g.ders[0];

    return result;
  }

  friend Jet sin(const Jet& f) {
    Jet result;

    result.ders[0] = sin(f.ders[0]);
    result.ders[1] = cos(f.ders[0]) * f.ders[1];
    result.ders[2] = cos(f.ders[0]) * f.ders[2];
    result.ders[3] = - sin(f.ders[0]) * f.ders[1] * f.ders[1] + cos(f.ders[0]) * f.ders[3];
    result.ders[4] = - sin(f.ders[0]) * f.ders[1] * f.ders[2] + cos(f.ders[0]) * f.ders[4];
    result.ders[5] = - sin(f.ders[0]) * f.ders[2] * f.ders[2] + cos(f.ders[0]) * f.ders[5];

    return result;
  }

  friend Jet cos(const Jet& f) {
    Jet result;

    result.ders[0] = cos(f.ders[0]);
    result.ders[1] = - (sin(f.ders[0]) * f.ders[1]);
    result.ders[2] = - (sin(f.ders[0]) * f.ders[2]);
    result.ders[3] = - (cos(f.ders[0]) * f.ders[1] * f.ders[1] + sin(f.ders[0]) * f.ders[3]);
    result.ders[4] = - (cos(f.ders[0]) * f.ders[1] * f.ders[2] + sin(f.ders[0]) * f.ders[4]);
    result.ders[5] = - (cos(f.ders[0]) * f.ders[2] * f.ders[2] + sin(f.ders[0]) * f.ders[5]);

    return result;
  }

  friend Jet exp(const Jet& f) {
    Jet result;

    result.ders[0] = exp(f.ders[0]);
    result.ders[1] = exp(f.ders[0]) * f.ders[1];
    result.ders[2] = exp(f.ders[0]) * f.ders[2];
    result.ders[3] = exp(f.ders[0]) * f.ders[1] * f.ders[1] + exp(f.ders[0]) * f.ders[3];
    result.ders[4] = exp(f.ders[0]) * f.ders[1] * f.ders[2] + exp(f.ders[0]) * f.ders[4];
    result.ders[5] = exp(f.ders[0]) * f.ders[2] * f.ders[2] + exp(f.ders[0]) * f.ders[5];

    return result;
  }

};

int main() {
  int M = 0;
  double x = 0;
  double y = 0;
  Jet result;
  
  cin >> M;
  for(int i = 0; i < M; i++) {
    cin >> x >> y;
    result = funkcja(Jet(x, true), Jet(y, false));
    for(int j = 0; j < 6; j++) {
      cout << fixed << setprecision(16) << result.ders[j];
      if(j != 5) {
        cout << " ";
      }
    }
    if(i != M - 1) {
      cout << endl;
    }
  }
}