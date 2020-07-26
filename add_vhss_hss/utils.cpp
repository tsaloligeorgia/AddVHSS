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

	/*std::cout << "POlynomial: ";
	 mpfq_p_127_1_poly_print( polynomial);
	 std::cout << std::endl;*/
	for (int i = 1; i < numberOfServers + 1; i++) {
		//mpz_init(evaluation_theta[i - 1]);

		mpz_class tmp(i);
		/*mpz_init(tmp);
		 mpz_set_ui(tmp, i);
		 */
		evaluation_theta[i - 1] = polynomialEvaluation(polynomial, tmp);
		/*std::cout << "evaluation of : " << i << " ";
		std::cout << evaluation_theta[i - 1];
		std::cout << std::endl;*/

	}

	for (int j = 1; j < numberOfServers + 1; j++) {
		mpz_class lambda_ij(1);
		/*mpz_init(lambda_ij);
		 mpz_set_ui(lambda_ij, 1);*/
		mpz_class tmp_j(j);
		/*mpz_init(tmp_j);
		 mpz_set_ui(tmp_j, j);*/
		for (int l = 1; l < numberOfServers + 1; l++) {
			if (l != j) {
				mpz_class tmp_l(l);
				/*mpz_init(tmp_l);
				 mpz_set_ui(tmp_l, l);
				 */
				mpz_class tmp_l_2(l);
				/*mpz_init(tmp_l_2);
				 mpz_set_ui(tmp_l_2, l);*/
				tmp_l_2 = tmp_l_2 - tmp_j;
				//mpz_sub(tmp_l_2, tmp_l_2, tmp_j);
				//mpfq_p_127_1_sub(tmp_l_2, tmp_l_2, tmp_j);

				tmp_l_2 = mod_inv(tmp_l_2);
				//mpfq_p_127_1_inv(tmp_l_2, tmp_l_2);

				tmp_l = (tmp_l * tmp_l_2) % this->k;
				/*	mpz_mul(tmp_l, tmp_l, tmp_l_2);
				 mpz_mod(tmp_l, tmp_l, this->k);
				 */
				//mpfq_p_127_1_mul(tmp_l, tmp_l, tmp_l_2);
				lambda_ij = (lambda_ij * tmp_l) % this->k;
				/*mpz_mul(lambda_ij, lambda_ij, tmp_l);
				 mpz_mod(lambda_ij, lambda_ij, this->k);*/
				//mpfq_p_127_1_mul(lambda_ij, lambda_ij, tmp_l);
				/*std::cout << "lambda_ij of : " << j << " ";
				 mpfq_p_127_1_print( lambda_ij);
				 std::cout << std::endl;*/

			}
		}
		lambda_ij = (lambda_ij * evaluation_theta[j - 1]) % this->k;
		/*mpz_mul(lambda_ij, lambda_ij, evaluation_theta[j - 1]);
		 mpz_mod(lambda_ij, lambda_ij, this->k);*/
		//mpfq_p_127_1_mul(lambda_ij, lambda_ij, evaluation_theta[j - 1]);
		pre_computed[j - 1] = lambda_ij;
		//mpz_set(pre_computed[j - 1], lambda_ij);
		/*	std::cout << "pre_computed[j - 1] : " << j << " ";
		 mpfq_p_127_1_print( pre_computed[j - 1]);
		 std::cout << std::endl;*/

	}

}

mpz_class Utils::mod_inv(mpz_class element) {
	mpz_class out;
	mpz_invert(out.get_mpz_t(), element.get_mpz_t(), this->k.get_mpz_t());
	return out;
}

mpz_class Utils::polynomialEvaluation(Polynomial polynomial, mpz_class x) {

	mpz_class tmp;
	//mpz_init(tmp);
	int degree = polynomial.getDegree(); //mpfq_p_127_1_poly_deg(this->polynomial);
	//mpfq_p_127_1_poly_getcoef(tmp, polynomial, degree); //
	tmp = polynomial.getCoeff(degree);

	// Evaluate value of polynomial using Horner's method
	for (int i = degree - 1; i >= 0; i--) {
		tmp = (tmp * x) % this->k;
		/*mpz_mul(tmp, tmp, x);
		 mpz_mod(tmp, tmp, this->k);
		 //mpfq_p_127_1_mul(tmp, tmp, x);*/

		mpz_class coeff;
		//mpz_init(coeff);
		coeff = polynomial.getCoeff(i);
		//	mpfq_p_127_1_poly_getcoef(coeff, polynomial, i);
		tmp = tmp + coeff;
		//mpz_add(tmp, tmp, coeff);
		//mpfq_p_127_1_add(tmp, tmp, coeff);

	}
	//*result = tmp;
	return tmp;
	//mpz_set(*result, tmp);

//result = &tmp;

	return 1;

}
