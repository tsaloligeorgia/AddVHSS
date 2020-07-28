#ifndef VHSSLHS_H_
#define VHSSLHS_H_

#include <gmpxx.h>
#include <iostream>
#include <map>
#include <time.h>

#include <vector>

#include "polynomial.h"
#include "utils.h"
#include "params.h"
#include "VerificationKey.h"
#include "SecretKey.h"
#include "Proof.h"

class VHSS_LHS {
public:
	VHSS_LHS();
	virtual ~VHSS_LHS();

	void key_gen(mpz_class p, mpz_class q, SecretKey *sk, VerificationKey *vk);

	void gen_shares_threshold(mpz_class secret_input, mpz_class random_e,
			std::map<int, mpz_class> *shares);

	mpz_class partial_eval(int j, std::map<int, mpz_class> shares_from_clients);

	mpz_class final_eval(std::map<int, mpz_class> partial_evaluations);

	void partial_proof(SecretKey *sk, VerificationKey *vk, mpz_class fid,
			mpz_class x_i, int i, Proof *sigma);

	void final_proof(VerificationKey *vk, mpz_class fid,
			std::vector<Proof> sigmas, Proof *final_proof);

	int verify(VerificationKey vk, Proof final_proof, mpz_class y);

	static mpz_class HashF(mpz_class element, mpz_class prime) {
		mpz_class tmp(element);
		mpz_class g(G);
		mpz_class result;
		//std::cout << "q  " << prime << std::endl;

		mpz_powm(result.get_mpz_t(), g.get_mpz_t(), element.get_mpz_t(),
				prime.get_mpz_t());
		//std::cout << "temp: "<< result << std::endl;
		mpz_class two(2);
		while (mpz_probab_prime_p(result.get_mpz_t(), 10) == 0
				|| (mpz_cmp(result.get_mpz_t(), two.get_mpz_t()) == 0)) {
			tmp += mpz_class(1);
			mpz_powm(result.get_mpz_t(), g.get_mpz_t(), tmp.get_mpz_t(),
					prime.get_mpz_t());
			//std::cout << "temp: "<< result << std::endl;

		}

		return result;
	}

	static mpz_class random_z_star(mpz_class n_hat, int i) {
		mpz_class random;

		unsigned long seed = rand() % 10 + 1+i;
		gmp_randstate_t rstate;
		gmp_randinit_mt(rstate);
		gmp_randseed_ui(rstate, seed);
		gmp_randclass rr(gmp_randinit_default);
		rr.seed(seed);
		random = rr.get_z_range(n_hat);
		mpz_class one(1);
		mpz_class zero(0);
		mpz_class result_gcd;

		mpz_gcd(result_gcd.get_mpz_t(), n_hat.get_mpz_t(), random.get_mpz_t());
		while (mpz_cmp(result_gcd.get_mpz_t(), one.get_mpz_t()) != 0
				&& (mpz_cmp(random.get_mpz_t(), zero.get_mpz_t()) == 0)) {
			random = rr.get_z_range(n_hat);
			mpz_gcd(result_gcd.get_mpz_t(), n_hat.get_mpz_t(),
					random.get_mpz_t());

		}
		gmp_randclear(rstate);
		return random;
	}

};

#endif /* VHSSLHS_H_ */
