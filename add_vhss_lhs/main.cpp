#include<bits/stdc++.h>
#include <gmpxx.h>

#include "params.h"
#include "Client.h"
#include "Server.h"
#include "utils.h"
#include "SecretKey.h"
#include "VerificationKey.h"
#include "VHSSLHS.h"
#include "Proof.h"

using namespace std;

void print_parameters() {
	std::cout << "-------------------------" << std::endl;
	std::cout << "--- \tParameters\t ---" << std::endl;
	std::cout << "--- \tNR_CLIENTS\t ---\t" << NR_CLIENTS << std::endl;
	std::cout << "--- \tNR_SERVERS\t ---\t" << NR_SERVERS << std::endl;
	std::cout << "--- \tFINITE_FIELD\t ---\t" << FINITE_FIELD << std::endl;
	std::cout << "--- \tSECURITY\t ---\t" << SECURITY << std::endl;
	std::cout << "-------------------------" << std::endl;
}

void print_secret(SecretKey *sk) {
	std::cout << "----Secret Key ----" << std::endl;
	std::cout << "----P_hat : " << sk->getP() << std::endl;
	std::cout << "----Q_hat : " << sk->getQ() << std::endl;
	std::cout << "---- --- ----" << std::endl;

}

void print_verification(VerificationKey *vk) {

	std::cout << "----Verification Key ----" << std::endl;
	std::cout << "----N : " << vk->getN() << std::endl;
	std::cout << "----N_hat : " << vk->getNHat() << std::endl;
	std::cout << "----G : " << vk->getG() << std::endl;
	std::cout << "----G1 : " << vk->getG1() << std::endl;
	std::cout << "----Hs: " << std::endl;
	for (size_t i = 0; i < vk->getHS().size(); i++) {
		std::cout << vk->getHS()[i] << std::endl;
	}
	std::cout << "---- --- ----" << std::endl;

}

// Driver program
int main() {
	print_parameters();

	VHSS_LHS lhs;
	mpz_class p = Utils::generate_safe_prime(SECURITY, mpz_class(2)); //Utils::generate_safe_prime(32, mpz_class(2));
	cout << "p: " << p << endl;
	mpz_class q = Utils::generate_safe_prime(SECURITY, p);
	while (mpz_cmp(p.get_mpz_t(), q.get_mpz_t()) == 0) {
		q = Utils::generate_safe_prime(SECURITY, p);
	}
	cout << "q: " << q << endl;

	mpz_class phi_N((p - 1) * (q - 1));
	mpz_class N(p * q);

	std::vector<Client> clients;
	std::vector<Server> servers;

	SecretKey sk;
	VerificationKey vk;
	lhs.key_gen(p, q, &sk, &vk);

	print_secret(&sk);
	print_verification(&vk);

	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		if (i != (NR_CLIENTS)) {
			Client c(i, mpz_class(i + 1), mpz_class(1), lhs);
			clients.push_back(c);
		} else {
			Client c(i, mpz_class(i + 1), mpz_class(7665587881), lhs);
			clients.push_back(c);
		}
	}

	for (int j = 1; j < NR_SERVERS + 1; j++) {
		Server s(j);
		servers.push_back(s);
	}

	cout << "Generating Shares " << endl;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		std::map<int, mpz_class> shared_keys;
		std::map<int, mpz_class> shares;

		c.generate_shares(&shares);
		for (int j = 0; j < NR_SERVERS; j++) {
			Server s = servers[j];
			s.setShare(c.getI(), shares[s.getJ()]);
			servers[j] = s;
		}

	}
	cout << "Finished Shares " << endl;

	cout << "Generate Partial Evals " << endl;

	std::map<int, mpz_class> partial_evals;
	for (int j = 0; j < NR_SERVERS; j++) {
		Server s = servers[j];
		mpz_class partial = lhs.partial_eval(s.getJ(), s.getShares());

		partial_evals.insert(std::pair<int, mpz_class>(s.getJ(), partial));
	}

	cout << "Finish Partial Evals " << endl;

	cout << "Generate Final Evals " << endl;
	mpz_class y = lhs.final_eval(partial_evals);
	cout << "Finish Final Evals " << endl;
	cout << "y : " << y << std::endl;

	cout << "Generate Partial Proofs  " << endl;
	std::vector<Proof> sigmas;
	mpz_class ris = 0;
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		Proof sigma;
		if (i != (NR_CLIENTS)) {
			mpz_class ri(1);
			ris = ris + ri;
			lhs.partial_proof(sk, vk, mpz_class(3), mpz_class(i + 1), i,
					&sigma);
			sigmas.push_back(sigma);
		} else {

			lhs.partial_proof(sk, vk, mpz_class(3), mpz_class(i + 1), i,
					&sigma);
			sigmas.push_back(sigma);

		}
	}

	cout << "Finish Partial Proofs  " << endl;
	Proof final_p;
	lhs.final_proof(vk, mpz_class(2), sigmas, &final_p);

	cout << "Generate Final Proof  " << endl;

	cout << "Finish Final Proof  " << endl;

	int result = lhs.verify(vk, final_p, y);
	if (result == 1) {
		std::cout << "PERFECT!!! " << std::endl;
		std::cout << "y = " << y << std::endl;
	} else {
		std::cout << "WRONG!!! " << std::endl;
	}
	return 0;
}

