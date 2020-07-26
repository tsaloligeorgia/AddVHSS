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

	static mpz_class generate_safe_prime(int size, mpz_class pp) {

		gmp_randstate_t state;
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL) + mpz_get_ui(pp.get_mpz_t());
		gmp_randseed_ui(state, tmp);

		mpz_t result, prime;

		mpz_init(result);
		mpz_init(prime);
		mpz_urandomb(result, state, size);

		mpz_nextprime(prime, result);

		while ((mpz_cmp(prime, pp.get_mpz_t()) == 0)) {
			mpz_urandomb(result, state, size);
			mpz_nextprime(prime, result);

		}

		gmp_randclear(state);

		mpz_class to_return(prime);

		mpz_clear(result);
		mpz_clear(prime);

		return to_return;
	}
	;

	static void generate_primeN(int size, mpz_class N, mpz_class *q,
			mpz_class *p, mpz_class *n_hat) {
		gmp_randstate_t state;
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL);
		gmp_randseed_ui(state, tmp);

		mpz_t result, tmp_prime, p_h, q_h, phi_n, gcd_r, p_h_1, q_h_1, n_hat_t;
		mpz_class one(1);

		mpz_init(result);
		mpz_init(tmp_prime);
		mpz_init(p_h);
		mpz_init(q_h);
		mpz_init(phi_n);
		mpz_init(gcd_r);
		mpz_init(p_h_1);
		mpz_init(q_h_1);
		mpz_init(n_hat_t);
		mpz_urandomb(result, state, size);

		mpz_nextprime(tmp_prime, result);

		mpz_mul_ui(p_h, tmp_prime, 2);
		mpz_add_ui(p_h, p_h, 1);
		do {
			do {

				mpz_urandomb(result, state, size);
				mpz_nextprime(tmp_prime, result);
				mpz_mul_ui(p_h, tmp_prime, 2);
				mpz_add_ui(p_h, p_h, 1);

			} while (mpz_probab_prime_p(p_h, 100) != 2);

			do {
				mpz_urandomb(result, state, size);
				mpz_nextprime(tmp_prime, result);
				mpz_mul_ui(q_h, tmp_prime, 2);
				mpz_add_ui(q_h, q_h, 1);

			} while (mpz_probab_prime_p(q_h, 100) != 2 || mpz_cmp(q_h, p_h) == 0);
			mpz_sub_ui(p_h_1, p_h, 1);
			mpz_sub_ui(q_h_1, q_h, 1);
			mpz_mul(phi_n, p_h_1, q_h_1);

			mpz_gcd(gcd_r, N.get_mpz_t(), phi_n);
		} while (mpz_cmp(gcd_r, one.get_mpz_t()) != 0);

		mpz_mul(n_hat_t, q_h, p_h);
		*q = mpz_class(q_h);
		*p = mpz_class(p_h);
		*n_hat = mpz_class(n_hat_t);

		gmp_randclear(state);

		mpz_clear(result);
		mpz_clear(tmp_prime);
		mpz_clear(p_h);
		mpz_clear(q_h);
		mpz_clear(phi_n);
		mpz_clear(gcd_r);
		mpz_clear(p_h_1);
		mpz_clear(q_h_1);
		mpz_clear(n_hat_t);

	}
	;

	static mpz_class generate_random(int size, mpz_class field) {
		gmp_randstate_t state;
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL) + size;
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
