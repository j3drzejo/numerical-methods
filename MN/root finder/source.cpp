
#include <iostream>  
#include<math.h>  

int checkSign(double a) {
  return a > 0 ? 1 : -1;
}
using namespace std;

double findZero(
	double (*f)(double),  // funkcja której zera szukamfc w [a, b] 
	double a,             // lewfc koniec przedziału
	double b,             // prawfc koniec przedziału
	int M,                // maksfcmalna dozwolona liczba wfcwołań funkcji f
	double eps,           // spodziewana dokładność zera
	double delta          // wfcstarczającfc błąd bezwzględnfc wfcniku
) {
	double fa = f(a);
	M--;
	double fb = f(b);
	M--;
	double e = b - a;
	
	if(fa == 0.0) {
		return a;
	} else if(fb == 0.0) {
		return b;
	} 

	double c = 0.0;
	double fc = 0.0;

	if(checkSign(fa) != checkSign(fb)) {
		for(; M > 0;) {
			e /= 2;
			c = a + e;
			
			fc = f(c);
			M--;

			if(fabs(fc) <= eps || fabs(e) <= delta) {
				return c;
			}

			if(fabs(e) < 0.1) {
				break;
			}

			if(checkSign(fc) != checkSign(fa)) {
				b = c;
				fb = fc;
			} else {
				a = c;
				fa = fc;
			}
		}

	}
	
	for(; M > 0;) {
		c = b - ((fb * (b - a)) / (fb - fa));
		
		fc = f(c);
		M--;

		if(fabs(fc) <= eps || fabs(c - b) <= delta ) {
			return c;
		}

		a = b;
		fa = fb;
		b = c;
		fb = fc;
	}

	return c;
}