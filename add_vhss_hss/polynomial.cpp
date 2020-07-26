#include "polynomial.h"

Polynomial::Polynomial(int degree, mpz_class p) {
	this->coef.resize(degree + 1);
	this->deg = degree;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	unsigned long tmp = time(NULL) + degree;
	gmp_randseed_ui(state, tmp);

	for (int i = 0; i < degree + 1; i++) {
		mpz_class z(200);
		mpz_class c;
		//mpz_inits(z, c);
		mpz_urandomb(z.get_mpz_t(), state, 2 * GMP_LIMB_BITS);
		/*std::cout << "Z: " << z << std::endl;
		 std::cout << "P: " << p << std::endl;*/
		c = z % p;
		coef[i] = c;
	}
	gmp_randclear(state);

}

Polynomial::~Polynomial() {
	// TODO Auto-generated destructor stub
}

mpz_class Polynomial::getCoeff(int pos) {
	return this->coef[pos];
	//return this->coef[pos];
}

void Polynomial::setCoeff(mpz_class value, int pos) {
	if (pos <= deg) {
		this->coef[pos] = value;
	} else {
		this->coef.resize(pos + 1);
		this->coef[pos] = value;
	}
}

int Polynomial::getDegree() {
	return this->deg;
}

void Polynomial::print_polynomial() {
	for (int i = 0; i < this->deg; i++) {
		if (coef[i])
			std::cout << coef[i] << "x^" << i << " + ";

	}
	std::cout << std::endl;
}

