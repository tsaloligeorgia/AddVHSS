
#ifndef VHSSHSS_H_
#define VHSSHSS_H_

#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <map>

#include "polynomial.h"
#include "utils.h"
#include "params.h"

class VHSS_HSS {
public:
	VHSS_HSS();
	virtual ~VHSS_HSS();

	void gen_shares(mpz_class secret_input, mpz_class random_e,
			std::map<int, mpz_class> *shares, mpz_class *tau);

	mpz_class partial_proof(std::map<int, mpz_class> shares);

	mpz_class final_proof(std::vector<mpz_class> sigmas);

	mpz_class partial_eval(int j, std::map<int, mpz_class> shares_from_clients);

	mpz_class final_eval(std::vector<mpz_class> partial_evaluations);

	int verify(std::vector<mpz_class> taus, mpz_class sigma, mpz_class y);
};

#endif /* VHSSHSS_H_ */
