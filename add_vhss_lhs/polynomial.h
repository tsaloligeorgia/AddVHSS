#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iostream>
#include <gmpxx.h>
#include <vector>

class Polynomial {
public:
	Polynomial(int degree, mpz_class p);
	virtual ~Polynomial();

	mpz_class getCoeff(int pos);

	void setCoeff(mpz_class value, int pos);

	int getDegree();

	void print_polynomial();
private:
	std::vector<mpz_class> coef;
	int deg;

};

#endif /* POLYNOMIAL_H_ */
