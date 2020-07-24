#ifndef VHSSTSS_H_
#define VHSSTSS_H_

#include <gmpxx.h>
#include <iostream>
#include <map>
#include <time.h>

#include <vector>

#include "Matrix.h"
#include "polynomial.h"
#include "utils.h"
#include "params.h"

class VHSS_TSS {
public:
	VHSS_TSS();
	virtual ~VHSS_TSS();

	void key_gen(mpz_class p, mpz_class q, mpz_class *pk, mpz_class *sk);

	void gen_shares_threshold(mpz_class secret_input, mpz_class sk,
			mpz_class random_e, mpz_class pub_key,
			std::map<int, mpz_class> *shares,
			std::map<int, mpz_class> *shared_key, Matrix *A_i, mpz_class *hash_h);
	mpz_class partial_eval(int j, std::map<int, mpz_class> shares_from_clients);

	mpz_class final_eval(std::map<int, mpz_class> partial_evaluations);

	std::map<int, std::map<int, mpz_class>> partial_proof(
			std::map<int, std::map<int, mpz_class>> omegas,
			std::vector<mpz_class> hash_hs, std::vector<Matrix> matrices,
			mpz_class N, mpz_class phi_N);

	mpz_class final_proof(
			std::map<int, std::map<int, mpz_class>> partial_proofs,
			std::vector<mpz_class> public_keys, std::vector<mpz_class> hash_hs,
			std::vector<Matrix> matrices, mpz_class N, mpz_class phi_N);

	int verify(std::vector<mpz_class> hash_hs, mpz_class final_proof, mpz_class y);
private:
	std::map<int, mpz_class> partial_proof_i(std::map<int, mpz_class> omega,
			mpz_class hash_H, Matrix A, mpz_class N, mpz_class phi_n);

	mpz_class final_proof_i(mpz_class pub_key, mpz_class hash,
			Matrix A, std::map<int, mpz_class> sigma, mpz_class N,
			mpz_class phi_n);

};

#endif /* VHSSTSS_H_ */
