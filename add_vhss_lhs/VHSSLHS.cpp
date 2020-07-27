#include "VHSSLHS.h"

VHSS_LHS::VHSS_LHS() {

}

VHSS_LHS::~VHSS_LHS() {

}

void VHSS_LHS::key_gen(mpz_class p, mpz_class q, SecretKey *secret_key,
		VerificationKey *vk) {
	mpz_class N = p * q;

	mpz_class p_hat_prime, q_hat_prime;

	mpz_class p_hat;
	mpz_class q_hat;
	mpz_class n_hat;
	Utils::generate_primeN(SECURITY/2, N, &p_hat, &q_hat, &n_hat);

	mpz_class g = VHSS_LHS::random_z_star(n_hat, 1);

	mpz_class g1 = VHSS_LHS::random_z_star(n_hat, 2);

	std::vector<mpz_class> hs;
	for (int i = 0; i < NR_CLIENTS; i++) {
		mpz_class random = VHSS_LHS::random_z_star(n_hat, i + 100);
		hs.push_back(random);
	}

	/*VerificationKey vk_t;
	SecretKey sk_t;*/

	vk->setN(N);

	vk->setNHat(n_hat);

	vk->setG(g);

	vk->setG1(g1);

	vk->setHs(hs);

	//*vk = vk_t;


	secret_key->setPHat(p_hat);
	secret_key->setQHat(q_hat);
	//*secret_key = sk_t;

}

void VHSS_LHS::gen_shares_threshold(mpz_class secret_input, mpz_class random_e,
		std::map<int, mpz_class> *shares) {

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

	*shares = s_tmp;

}

mpz_class VHSS_LHS::partial_eval(int j,
		std::map<int, mpz_class> shares_from_clients) {
	mpz_class partial_eval(0);
	for (int i = 0; i < NR_CLIENTS; i++) {
		partial_eval += shares_from_clients[i + 1];
	}
	return partial_eval;

}

mpz_class VHSS_LHS::final_eval(std::map<int, mpz_class> partial_evaluations) {
	mpz_class final_eval(0);
	for (int i = 0; i < NR_SERVERS; i++) {
		final_eval += partial_evaluations[i + 1];

	}
	final_eval = final_eval % mpz_class(FINITE_FIELD);
	return final_eval;

}

void VHSS_LHS::partial_proof(SecretKey sk, VerificationKey vk, mpz_class fid,
		mpz_class x_i_R, int i, Proof *sigma) {

	mpz_class g = vk.getG();
	mpz_class n_hat = vk.getNHat();
	mpz_class g1 = vk.getG1();
	mpz_class hi = vk.getHS()[i - 1];

	mpz_class phi = (sk.getP() - mpz_class(1)) * (sk.getQ() - mpz_class(1));

	mpz_class e = VHSS_LHS::HashF(fid, mpz_class(FINITE_FIELD));
	mpz_class e_n = e * vk.getN();

	unsigned long seed = rand() % 10 + 1;
	gmp_randstate_t rstate;
	gmp_randinit_mt(rstate);
	gmp_randseed_ui(rstate, seed);
	gmp_randclass rr(gmp_randinit_default);
	rr.seed(seed);
	mpz_class s_i = rr.get_z_bits(2048) % e_n;

	mpz_class s_i_pow;

	mpz_powm(s_i_pow.get_mpz_t(), g.get_mpz_t(), s_i.get_mpz_t(),
			n_hat.get_mpz_t());

	mpz_class g1_pow;

	mpz_powm(g1_pow.get_mpz_t(), g1.get_mpz_t(), x_i_R.get_mpz_t(),
			n_hat.get_mpz_t());

	mpz_class right_hand_side = s_i_pow * g1_pow;
	right_hand_side = (right_hand_side * hi) % n_hat;

	mpz_class inv_e_N;

	mpz_invert(inv_e_N.get_mpz_t(), e_n.get_mpz_t(), phi.get_mpz_t());

	mpz_class x;

	mpz_powm(x.get_mpz_t(), right_hand_side.get_mpz_t(), inv_e_N.get_mpz_t(),
			n_hat.get_mpz_t());

	/*	std::cout << " n_hat " << n_hat << " x = " << x << " e_n = " << e_n
	 << std::endl;
	 std::cout << " right_hand_side " << right_hand_side << " phi = " << phi
	 << " g = " << g << std::endl;
	 std::cout << " hi " << hi << " g1 = " << g1 << " si = " << s_i << std::endl;

	 std::cout << " x_i_r " << x_i_R << " g1_pow = " << g1_pow << " s_i_pow = "
	 << s_i_pow << std::endl;*/

	Proof proof;
	proof.setE(e);
	proof.setSi(s_i);
	proof.setX_i(x);

	*sigma = proof;

	gmp_randclear(rstate);

}

void VHSS_LHS::final_proof(VerificationKey vk, mpz_class fid,
		std::vector<Proof> sigmas, Proof *final_proof) {
	mpz_class final_proof_result;

	mpz_class e = VHSS_LHS::HashF(fid, mpz_class(FINITE_FIELD));
	mpz_class e_n = e * vk.getN();

	mpz_class prod_partial_proof(1);
	for (unsigned int i = 0; i < sigmas.size(); i++) {
		prod_partial_proof = (prod_partial_proof * sigmas[i].getX());
	}
	prod_partial_proof = prod_partial_proof % vk.getNHat();

	mpz_class sum_s_i(0);
	for (unsigned int i = 0; i < sigmas.size(); i++) {
		sum_s_i = sum_s_i + sigmas[i].getSi();
	}
	mpz_class s = sum_s_i % e_n;

	mpz_class tmp = sum_s_i - s;

	mpz_class inv_t;
	mpz_invert(inv_t.get_mpz_t(), e_n.get_mpz_t(), vk.getNHat().get_mpz_t());
	mpz_class s_prime = tmp * inv_t % vk.getNHat();

	mpz_class low_part;
	mpz_powm(low_part.get_mpz_t(), vk.getG().get_mpz_t(), s_prime.get_mpz_t(),
			vk.getNHat().get_mpz_t());

	mpz_class inv_low;
	mpz_invert(inv_low.get_mpz_t(), low_part.get_mpz_t(),
			vk.getNHat().get_mpz_t());

	mpz_class x_tilda = (prod_partial_proof * inv_low) % vk.getNHat();

	/*
	 std::cout << " e_n " << e_n << " x_tilda = " << x_tilda << " low_part = "
	 << low_part << std::endl;
	 std::cout << " inv_low " << inv_low << " sum_s_i = " << sum_s_i << " tmp = "
	 << tmp << std::endl;
	 std::cout << " prod_partial_proof " << prod_partial_proof << " s = " << s
	 << std::endl;
	 */

	Proof final;
	final.setE(e);

	final.setSi(s);

	final.setX_i(x_tilda);

	*final_proof = final;

}

int VHSS_LHS::verify(VerificationKey vk, Proof final_proof, mpz_class y) {
	mpz_class prod(1);
	for (int i = 0; i < NR_CLIENTS; i++) {
		prod = prod * vk.getHS()[i];
	}

	mpz_class g_power_s;
	mpz_powm(g_power_s.get_mpz_t(), vk.getG().get_mpz_t(),
			final_proof.getSi().get_mpz_t(), vk.getNHat().get_mpz_t());

	mpz_class g1_power_y;
	mpz_powm(g1_power_y.get_mpz_t(), vk.getG1().get_mpz_t(), y.get_mpz_t(),
			vk.getNHat().get_mpz_t());

	mpz_class right_part = (g_power_s * g1_power_y * prod) % vk.getNHat();

	mpz_class e_n = final_proof.getE() * vk.getN();
	mpz_class left_part;
	mpz_powm(left_part.get_mpz_t(), final_proof.getX().get_mpz_t(),
			e_n.get_mpz_t(), vk.getNHat().get_mpz_t());

	std::cout << "Left part: " << left_part << std::endl;
	std::cout << "right_part: " << right_part << std::endl;

	if (mpz_cmp(left_part.get_mpz_t(), right_part.get_mpz_t()) == 0)
		return 1;

	return 0;

}

