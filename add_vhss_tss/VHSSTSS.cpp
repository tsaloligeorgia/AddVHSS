
#include "VHSSTSS.h"

VHSS_TSS::VHSS_TSS() {

}

VHSS_TSS::~VHSS_TSS() {

}

void VHSS_TSS::key_gen(mpz_class p, mpz_class q, mpz_class *pk, mpz_class *sk) {

	unsigned long seed = rand() % 10 + 1;

	gmp_randstate_t rstate;
	// initialize state for a Mersenne Twister algorithm. This algorithm is fast and has good randomness properties.
	gmp_randinit_mt(rstate);

	// create the generator seed for the random engine to reference
	gmp_randseed_ui(rstate, seed);

	mpz_class p_prime((p - 1) / 2);
	mpz_class q_prime((q - 1) / 2);
	mpz_class N_prime(2 * (q_prime * p_prime));
	gmp_randclass rr(gmp_randinit_default);
	rr.seed(seed);
	mpz_class result_gcd;
	mpz_class one(1);

	mpz_class e = rr.get_z_range(N_prime);
	mpz_gcd(result_gcd.get_mpz_t(), e.get_mpz_t(), N_prime.get_mpz_t());
	while (mpz_cmp(result_gcd.get_mpz_t(), one.get_mpz_t()) != 0) {
		e = rr.get_z_range(N_prime);
		mpz_gcd(result_gcd.get_mpz_t(), e.get_mpz_t(), N_prime.get_mpz_t());
	}
	mpz_class d;
	mpz_invert(d.get_mpz_t(), e.get_mpz_t(), N_prime.get_mpz_t());

	*pk = e;
	*sk = d;
	gmp_randclear(rstate);
}

void VHSS_TSS::gen_shares_threshold(mpz_class secret_input, mpz_class sk,
		mpz_class random_e, mpz_class pub_key, std::map<int, mpz_class> *shares,
		std::map<int, mpz_class> *shared_key, Matrix *A_i, mpz_class *hash_h) {

	Polynomial polynomial(THRESHOLD_CLIENTS, mpz_class(FINITE_FIELD));
	polynomial.setCoeff(secret_input, 0);

	//polynomial.print_polynomial();

	//mpfq_p_127_1_poly_print(this->k, polynomial);

	Utils util(FINITE_FIELD);

	mpz_class s[NR_SERVERS];
	for (int i = 0; i < NR_SERVERS; i++) {
		s[i] = mpz_class(0);
	}

	util.computeProducts(polynomial, NR_SERVERS, s);

	std::map<int, mpz_class> s_tmp;
	for (int i = 1; i < NR_SERVERS + 1; i++) {
		s_tmp.insert(std::pair<int, mpz_class>(i, s[i - 1]));
	}

	Matrix A_i_tmp = Matrix(NR_SERVERS, THRESHOLD, true);
	Matrix A_i_t = A_i_tmp.submatrix(0, THRESHOLD, 0, THRESHOLD);

	mpz_class det = 2 * A_i_t.determinant();

	mpz_class result_gcd;
	mpz_class one(1);

	mpz_gcd(result_gcd.get_mpz_t(), det.get_mpz_t(), pub_key.get_mpz_t());
	while (mpz_cmp(result_gcd.get_mpz_t(), one.get_mpz_t()) != 0) {
		A_i_tmp = Matrix(NR_SERVERS, THRESHOLD, true);
		A_i_t = A_i_tmp.submatrix(0, THRESHOLD, 0, THRESHOLD);

		det = 2 * A_i_t.determinant();

		mpz_gcd(result_gcd.get_mpz_t(), det.get_mpz_t(), pub_key.get_mpz_t());

	}
	gmp_randclass rr(gmp_randinit_default);

	std::vector<mpz_class> vec_d;
	vec_d.push_back(sk);
	for (int i = 0; i < THRESHOLD; i++) {
		vec_d.push_back(rr.get_z_range(pub_key + sk));

	}

	std::vector<mpz_class> omega = A_i_tmp.multiplyByVector(vec_d);
	std::map<int, mpz_class> omegas;
	int i = 1;
	for (mpz_class o : omega) {
		omegas.insert(std::pair<int, mpz_class>(i, o));
		i++;
	}

	mpz_class exp = (secret_input + random_e) % mpz_class(FINITE_FIELD - 1);
	mpz_class h_i;
	mpz_powm(h_i.get_mpz_t(), mpz_class(G).get_mpz_t(), exp.get_mpz_t(),
			mpz_class(FINITE_FIELD).get_mpz_t());

	*hash_h = h_i;
	*A_i = A_i_tmp;
	*shared_key = omegas;
	*shares = s_tmp;

}

mpz_class VHSS_TSS::partial_eval(int j,
		std::map<int, mpz_class> shares_from_clients) {
	mpz_class partial_eval(0);
	for (int i = 0; i < NR_CLIENTS; i++) {
		partial_eval += shares_from_clients[i + 1];
	}
	return partial_eval;

}

mpz_class VHSS_TSS::final_eval(std::map<int, mpz_class> partial_evaluations) {
	mpz_class final_eval(0);
	for (int i = 0; i < NR_SERVERS; i++) {
		final_eval += partial_evaluations[i + 1];

	}
	final_eval = final_eval % mpz_class(FINITE_FIELD);
	return final_eval;

}
std::map<int, mpz_class> VHSS_TSS::partial_proof_i(
		std::map<int, mpz_class> omega, mpz_class hash_H, Matrix A, mpz_class N,
		mpz_class phi_n) {
	std::map<int, mpz_class> result;

	Matrix matrix_a_s = A.submatrix(0, THRESHOLD, 0, THRESHOLD);
	Matrix matrix_c_s_adjugate = matrix_a_s.adjoint();
	for (int j = 1; j < THRESHOLD + 1; j++) {
		mpz_class tmp_val = matrix_c_s_adjugate.getElement(0, j - 1);
		mpz_class tmp = omega[j];
		mpz_class exponent = (2 * tmp_val * tmp) % phi_n;
		mpz_class sigma;
		mpz_powm(sigma.get_mpz_t(), hash_H.get_mpz_t(), exponent.get_mpz_t(),
				N.get_mpz_t());
		result.insert(std::pair<int, mpz_class>(j, sigma));
	}

	return result;

}

std::map<int, std::map<int, mpz_class>> VHSS_TSS::partial_proof(
		std::map<int, std::map<int, mpz_class>> omegas,
		std::vector<mpz_class> hash_hs, std::vector<Matrix> matrices,
		mpz_class N, mpz_class phi_N) {
	//std::cout << "Starting partial_proof" << std::endl;
	std::map<int, std::map<int, mpz_class>> result;

	for (int i = 0; i < NR_CLIENTS; i++) {
		std::map<int, mpz_class> r = this->partial_proof_i(omegas[i + 1],
				hash_hs[i], matrices[i], N, phi_N);
		result.insert(std::pair<int, std::map<int, mpz_class>>(i + 1, r));

	}
	//std::cout << "End partial_proof" << std::endl;

	return result;

}

mpz_class VHSS_TSS::final_proof_i(mpz_class pub_key, mpz_class hash, Matrix A,
		std::map<int, mpz_class> sigma, mpz_class N, mpz_class phi_n) {

	mpz_class sigma_bar(1);
	for (int j = 1; j < THRESHOLD + 1; j++) {
		sigma_bar *= sigma[j] % N;

	}
	Matrix matrix_a_s = A.submatrix(0, THRESHOLD, 0, THRESHOLD);
	mpz_class detal_a_is = matrix_a_s.determinant();
	mpz_class tmp(2 * detal_a_is);
	mpz_class alpha, beta, rest;
	//a*s + b*t = g
	//mpz_gcdext (mpz_t g, mpz_t s, mpz_t t, const mpz_t a, const mpz_t b);
	mpz_gcdext(rest.get_mpz_t(), alpha.get_mpz_t(), beta.get_mpz_t(),
			tmp.get_mpz_t(), pub_key.get_mpz_t());
	alpha = alpha % phi_n;
	beta = beta % phi_n;

	mpz_class sigma_tmp_1;
	mpz_powm(sigma_tmp_1.get_mpz_t(), sigma_bar.get_mpz_t(), alpha.get_mpz_t(),
			N.get_mpz_t());

	mpz_class sigma_tmp_2;
	mpz_powm(sigma_tmp_2.get_mpz_t(), hash.get_mpz_t(), beta.get_mpz_t(),
			N.get_mpz_t());

	mpz_class mpz_result = (sigma_tmp_1 * sigma_tmp_2) % N;

	return mpz_result;

}

mpz_class VHSS_TSS::final_proof(
		std::map<int, std::map<int, mpz_class>> partial_proofs,
		std::vector<mpz_class> public_keys, std::vector<mpz_class> hash_hs,
		std::vector<Matrix> matrices, mpz_class N, mpz_class phi_N) {

	std::vector<mpz_class> tmp;
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		mpz_class r_tmp = this->final_proof_i(public_keys[i - 1],
				hash_hs[i - 1], matrices[i - 1], partial_proofs[i], N, phi_N);
		tmp.push_back(r_tmp);
	}

	mpz_class final_result(1);
	int i = 0;
	for (mpz_class sigma : tmp) {
		mpz_class tmp_pow;
		mpz_powm(tmp_pow.get_mpz_t(), sigma.get_mpz_t(),
				public_keys[i].get_mpz_t(), N.get_mpz_t());
		final_result = final_result * tmp_pow;
		i++;
	}
	return final_result;

}

int VHSS_TSS::verify(std::vector<mpz_class> hash_hs, mpz_class final_proof,
		mpz_class y) {
	mpz_class prod(1);
	for (mpz_class hash : hash_hs) {
		prod = prod * hash;
	}
	prod = prod % mpz_class(FINITE_FIELD);
	mpz_class sigma = final_proof % mpz_class(FINITE_FIELD);
	mpz_class h_y;
	mpz_powm(h_y.get_mpz_t(), mpz_class(G).get_mpz_t(), y.get_mpz_t(),
			mpz_class(FINITE_FIELD).get_mpz_t());
	/*
	 std::cout << " prod: " << prod << std::endl;

	 std::cout << "h_y : " << h_y << std::endl;

	 std::cout << "sigma : " << sigma << std::endl;*/
	if (prod == sigma) {
		if (h_y == sigma) {
			return 1;
		}
	}

	return 0;

}

