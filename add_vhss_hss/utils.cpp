#include "utils.h"

Utils::Utils() {

}

Utils::Utils(mpz_class k) {
	this->k = k;

}

Utils::~Utils() {
}

void Utils::computeProducts(Polynomial polynomial, int numberOfServers,
		mpz_class pre_computed[]) {
	mpz_class evaluation_theta[numberOfServers];

	for (int i = 1; i < numberOfServers + 1; i++) {
		mpz_class tmp(i);
		evaluation_theta[i - 1] = polynomialEvaluation(polynomial, tmp);

	}
	for (int j = 1; j < numberOfServers + 1; j++) {
		mpz_class lambda_ij(1);
		mpz_class tmp_j(j);
		for (int l = 1; l < numberOfServers + 1; l++) {
			if (l != j) {
				mpz_class tmp_l(l);
				mpz_class tmp_l_2(l);

				tmp_l_2 = tmp_l_2 - tmp_j;

				tmp_l_2 = mod_inv(tmp_l_2);

				tmp_l = (tmp_l * tmp_l_2) % this->k;
				lambda_ij = (lambda_ij * tmp_l) % this->k;

			}
		}
		lambda_ij = (lambda_ij * evaluation_theta[j - 1]) % this->k;
		pre_computed[j - 1] = lambda_ij;

	}

}

mpz_class Utils::mod_inv(mpz_class element) {
	mpz_class out;
	mpz_invert(out.get_mpz_t(), element.get_mpz_t(), this->k.get_mpz_t());
	return out;
}

mpz_class Utils::polynomialEvaluation(Polynomial polynomial, mpz_class x) {

	mpz_class tmp;
	int degree = polynomial.getDegree();
	tmp = polynomial.getCoeff(degree);

	// Evaluate value of polynomial using Horner's method
	for (int i = degree - 1; i >= 0; i--) {
		tmp = (tmp * x) % this->k;

		mpz_class coeff;
		coeff = polynomial.getCoeff(i);
		tmp = tmp + coeff;

	}
	return tmp;


}
