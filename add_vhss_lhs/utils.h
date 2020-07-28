#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "polynomial.h"
#include "params.h"

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

		mpz_probab_safe_prime_p_next(prime, result, 50);

		while ((mpz_cmp(prime, pp.get_mpz_t()) == 0)) {
			mpz_urandomb(result, state, size);
			mpz_probab_safe_prime_p_next(prime, result, 50);

		}

		gmp_randclear(state);

		mpz_class to_return(prime);

		mpz_clear(result);
		mpz_clear(prime);

		return to_return;
	}
	;

	static int mpz_probab_safe_prime_p(mpz_t n, int reps) {
		mpz_t m;
		int res;

		if (!mpz_probab_prime_p(n, reps + 1)) {
			return 0;
		}

		mpz_init(m);
		mpz_sub_ui(m, n, 1L);
		mpz_div_ui(m, m, 2L);

		res = mpz_probab_prime_p(m, reps + 1);

		mpz_clear(m);

		return res;
	}

	static void mpz_probab_safe_prime_p_next(mpz_t rop, mpz_t n, int reps) {
		int increased = 0;

		mpz_set(rop, n);

		if (mpz_cmp_ui(rop, 5) < 0) {
			mpz_set_ui(rop, 5);
		} else if (mpz_cmp_ui(rop, 7) < 0) {
			mpz_set_ui(rop, 7);
		} else {

			/* Make sure that rop is odd. */
			if (!mpz_tstbit(rop, 0)) {
				mpz_add_ui(rop, rop, 1L);
				increased = 1;
			}

			/* Make sure that m is odd, where rop=2m+1. */
			if (!mpz_tstbit(rop, 1)) {
				mpz_add_ui(rop, rop, 2L);
				increased = 1;
			}

			/* If both rop and m were already odd, then we add 4. */
			if (!increased) {
				mpz_add_ui(rop, rop, 4L);
			}

			while (mpz_probab_safe_prime_p(rop, reps) == 0) {
				mpz_add_ui(rop, rop, 4L);
			}
		}
	}

	static void generate_primeN(int size, mpz_class N, mpz_class *q,
			mpz_class *p, mpz_class *n_hat) {

		gmp_randstate_t state;
		gmp_randinit_mt(state);
		unsigned long tmp = rand();
		gmp_randseed_ui(state, tmp);
		mpz_class safe_prime_p;
		mpz_class safe_prime_q;
		mpz_class one(1);

		mpz_t result, gcd_r, phi_n, p_t, q_t;
		mpz_init(result);
		mpz_init(gcd_r);
		mpz_init(phi_n);
		mpz_init(p_t);
		mpz_init(q_t);

		do {
			mpz_urandomb(result, state, size);
			mpz_probab_safe_prime_p_next(safe_prime_p.get_mpz_t(), result, 50);

			do {
				mpz_urandomb(result, state, size);
				mpz_probab_safe_prime_p_next(safe_prime_q.get_mpz_t(), result,
						50);
			} while (mpz_cmp(safe_prime_p.get_mpz_t(), safe_prime_q.get_mpz_t())
					== 0);

			mpz_sub_ui(p_t, safe_prime_p.get_mpz_t(), 1L);
			mpz_sub_ui(q_t, safe_prime_q.get_mpz_t(), 1L);
			mpz_mul(phi_n, p_t, q_t);

			mpz_gcd(gcd_r, N.get_mpz_t(), phi_n);

		} while (mpz_cmp(gcd_r, one.get_mpz_t()) != 0);

		*p = mpz_class(safe_prime_p);
		*q = mpz_class(safe_prime_q);
		*n_hat = mpz_class(safe_prime_p*safe_prime_q);

		mpz_clear(result);
		mpz_clear(gcd_r);
		mpz_clear(phi_n);
		mpz_clear(p_t);
		mpz_clear(q_t);
		gmp_randclear(state);
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

	static mpz_class random_element() {

		gmp_randstate_t state;
		mpz_class field(FINITE_FIELD);
		gmp_randinit_default(state);
		unsigned long tmp = time(NULL) + 32;
		gmp_randseed_ui(state, tmp);
		mpz_t result;
		mpz_init(result);
		mpz_urandomb(result, state, 32);

		mpz_class to_return(result);

		gmp_randclear(state);
		mpz_clear(result);

		return to_return % field;

	}

private:

	mpz_class mod_inv(mpz_class element);
	mpz_class k;
};

#endif /* UTILS_H_ */
