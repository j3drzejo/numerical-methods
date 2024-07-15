
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
// S = a + aq + aq^2
// P = a * aq * aq^2, so q = (P^(1/3))/(a), after some manipulation 
// 0 = a^2 + (P^(1/3) - S)a + (P^(1/3))^2

float getAns(float firstAns, float secondAns, float cbrtprod, float sum) {
  float q1 = cbrtprod / firstAns;
  float q2 = cbrtprod / secondAns;

  if(firstAns > powf(q1, 2) * firstAns) {
    return firstAns;
  } else {
    return secondAns;
  }
}

int main() {
  int numOfSessions = 0;
  cin >> numOfSessions;

  for(int i = 0; i < numOfSessions; i++) {
    float sum = 0.0;
    float product = 0.0;

    cin >> product >> sum;

    // formula from lecture, but now using delta to find roots 
    float cbrtprod = cbrtf(product);
    float p = (-(cbrtprod) + (sum));
    float q = powf(cbrtprod, 2);

    
    float sqrDelta = sqrtf(fabs(sum - 3.0 * cbrtprod)) * sqrtf(fabs(sum + cbrtprod));

    if(p != 0 && q != 0 && (((sum - 3.0 * cbrtprod >= 0) && (sum + cbrtprod) >= 0)) || ((sum - 3.0 * cbrtprod <= 0) && (sum + cbrtprod) <= 0) ) {
      float x1 = 0.0;
      float x2 = 0.0;

      if(p >= 0) {
        x1 = (p + sqrDelta) / 2;
        x2 = q / x1;
      } else {
        x2 = (p - sqrDelta) / 2;
        x1 = q / x2;
      }
      float Ans = getAns(x1, x2, cbrtprod, sum);

      float geom = cbrtprod / Ans;

      if(isnan(Ans) || isinf(Ans) || ((Ans == 0 || Ans * geom == 0 || Ans * geom * geom == 0) && (Ans != 0 || Ans * geom != 0 || Ans * geom * geom != 0))) {
        Ans = 0.0;
        geom = 0.0;
      }

      cout << fixed << scientific << setprecision(10) << Ans << " " << Ans * geom << " " << Ans * geom * geom;
    } else {
      cout << fixed << scientific << setprecision(10) << 0.0 << " " << 0.0 << " " << 0.0;
    }
    if(i != numOfSessions - 1) {
      cout << endl;
    }
  }
}