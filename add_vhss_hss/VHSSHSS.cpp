#include "VHSSHSS.h"

VHSS_HSS::VHSS_HSS() {
	// TODO Auto-generated constructor stub

}

VHSS_HSS::~VHSS_HSS() {
	// TODO Auto-generated destructor stub
}

void VHSS_HSS::gen_shares(mpz_class secret_input, mpz_class random_e,
		std::map<int, mpz_class> *shares, mpz_class *tau) {

	Polynomial polynomial(THRESHOLD_CLIENTS, mpz_class(FINITE_FIELD));
	polynomial.setCoeff(secret_input, 0);

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

	mpz_class exp = secret_input + random_e;

	mpz_class tau_t;
	mpz_class g(G);
	mpz_class prime(FINITE_FIELD);
	mpz_powm(tau_t.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(),
			prime.get_mpz_t());

	*shares = s_tmp;
	*tau = tau_t;

}

mpz_class VHSS_HSS::partial_eval(int j,
		std::map<int, mpz_class> shares_from_clients) {
	mpz_class partial_eval(0);
	for (int i = 0; i < NR_CLIENTS; i++) {
		partial_eval += shares_from_clients[i + 1];
	}
	return partial_eval;

}

mpz_class VHSS_HSS::final_eval(std::vector<mpz_class> partial_evaluations) {
	mpz_class final_eval(0);
	for (int i = 0; i < NR_SERVERS; i++) {
		final_eval += partial_evaluations[i];

	}
	return final_eval;

}

mpz_class VHSS_HSS::partial_proof(std::map<int, mpz_class> shares) {
	mpz_class y_j(0);

	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		y_j = y_j + shares[i];
	}
	mpz_class sigma_temp;
	mpz_class g(G);
	mpz_class prime(FINITE_FIELD);

	//mpz_pow_ui(sigma_temp.get_mpz_t(), g.get_mpz_t(), y_j.get_ui());
	mpz_powm(sigma_temp.get_mpz_t(), g.get_mpz_t(), y_j.get_mpz_t(),
	 prime.get_mpz_t());

	return sigma_temp;

}

mpz_class VHSS_HSS::final_proof(std::vector<mpz_class> sigmas) {
	mpz_class final_proof(1);
	for (mpz_class sigma : sigmas) {
		final_proof = final_proof * sigma;
	}
	return final_proof;
}

int VHSS_HSS::verify(std::vector<mpz_class> taus, mpz_class sigma,
		mpz_class y) {
	mpz_class prod(1);
	for (mpz_class tau : taus) {
		prod = prod * tau;

	}

	mpz_class modQ(FINITE_FIELD);

	prod = prod % modQ;

	mpz_class y_mod = y % modQ;
	mpz_class sigma_mod = sigma % modQ;

	mpz_class hash_y;
	mpz_class g(G);
	mpz_powm(hash_y.get_mpz_t(), g.get_mpz_t(), y.get_mpz_t(),
			modQ.get_mpz_t());

	mpz_class hash_y_mod;

	mpz_powm(hash_y_mod.get_mpz_t(), g.get_mpz_t(), y_mod.get_mpz_t(),
			modQ.get_mpz_t());

	if (mpz_cmp(sigma_mod.get_mpz_t(), hash_y.get_mpz_t()) == 0) {
		if (mpz_cmp(prod.get_mpz_t(), hash_y_mod.get_mpz_t()) == 0) {
			return 1;
		} else {
			std::cout << "prod ne hash_y" << std::endl;
		}
	} else {
		std::cout << "sigma ne hash_y" << std::endl;
	}

	return 0;

}

