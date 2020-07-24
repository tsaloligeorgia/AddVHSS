#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "polynomial.h"

class Utils {

public:
	Utils();
	Utils(mpz_class k);
	virtual ~Utils();

	void computeProducts(Polynomial polynomial, int nr_of_servers,
			mpz_class pre_computed[]);
	mpz_class polynomialEvaluation(Polynomial polynomial, mpz_class x);

	static mpz_class random_prime(int size) {
		gmp_randstate_t state;
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL);
		gmp_randseed_ui(state, tmp);

		mpz_t result, prime;
		mpz_init(result);
		mpz_init(prime);
		mpz_urandomb(result, state, size);

		mpz_nextprime(prime, result);
		gmp_randclear(state);

		mpz_class to_return(prime);

		mpz_clear(result);
		mpz_clear(prime);

		return to_return;
	}
	;

	static mpz_class generate_random(int size, mpz_class field) {
		gmp_randstate_t state;
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL)+size;
		gmp_randseed_ui(state, tmp);
		mpz_t result;
		mpz_init(result);
		mpz_urandomb(result, state, size);

		mpz_class to_return(result);

		gmp_randclear(state);
		mpz_clear(result);

		return to_return % field;
	}
	;

private:

	mpz_class mod_inv(mpz_class element);
	mpz_class k;
};

#endif /* UTILS_H_ */
