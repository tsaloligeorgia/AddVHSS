#include<bits/stdc++.h>
#include <gmpxx.h>
#include <vector>
#include <map>
#include <chrono>

#include "params.h"
#include "Client.h"
#include "Server.h"
#include "utils.h"
#include "SecretKey.h"
#include "VerificationKey.h"
#include "VHSSLHS.h"
#include "Proof.h"

using namespace std;
using namespace std::chrono;

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

	auto start_key_setup = high_resolution_clock::now();
	lhs.key_gen(p, q, &sk, &vk);
	auto final_key_setup = high_resolution_clock::now();
	auto duration_key_setup = duration_cast<microseconds>(
			final_key_setup - start_key_setup);
	cout << "Time taken by generate_shares: " << duration_key_setup.count()
			<< " microseconds" << endl;

	print_secret(&sk);
	print_verification(&vk);

	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		if (i != (NR_CLIENTS)) {
			Client c(i, mpz_class(i + 1), mpz_class(0), lhs);
			clients.push_back(c);
		} else {
			Client c(i, mpz_class(i + 1), mpz_class(0), lhs);
			clients.push_back(c);
		}
	}

	for (int j = 1; j < NR_SERVERS + 1; j++) {
		Server s(j);
		servers.push_back(s);
	}

	cout << "Generating Shares " << endl;
	std::vector<microseconds> gen_shares_timing;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		std::map<int, mpz_class> shared_keys;
		std::map<int, mpz_class> shares;


		auto start_gen_shares = high_resolution_clock::now();
		c.generate_shares(&shares);
		auto final_gen_shares = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(
				final_gen_shares - start_gen_shares);
		gen_shares_timing.push_back(duration);
		for (int j = 0; j < NR_SERVERS; j++) {
			Server s = servers[j];
			s.setShare(c.getI(), shares[s.getJ()]);
			servers[j] = s;
		}

	}

	cout << "Finished Shares " << endl;
	sort(gen_shares_timing.begin(), gen_shares_timing.end());
	cout << "Time taken by generate_shares: " << gen_shares_timing[0].count()
			<< " microseconds" << endl;

	cout << "Generate Partial Evals " << endl;

	std::map<int, mpz_class> partial_evals;
	for (int j = 0; j < NR_SERVERS; j++) {
		Server s = servers[j];
		auto start_partial_eval = high_resolution_clock::now();
		mpz_class partial = lhs.partial_eval(s.getJ(), s.getShares());
		auto final_partial_eval = high_resolution_clock::now();
		auto duration_eval = duration_cast<microseconds>(
				final_partial_eval - start_partial_eval);

		cout << "Time taken by partial_eval: " << duration_eval.count()
				<< " microseconds" << endl;

		partial_evals.insert(std::pair<int, mpz_class>(s.getJ(), partial));
	}

	cout << "Finish Partial Evals " << endl;

	cout << "Generate Final Evals " << endl;
	auto start_final_eval = high_resolution_clock::now();
	mpz_class y = lhs.final_eval(partial_evals);
	auto final_final_eval = high_resolution_clock::now();
	auto duration_fina_eval = duration_cast<nanoseconds>(
			final_final_eval - start_final_eval);
	cout << "Finish Final Evals " << endl;
	cout << "Time taken by final_eval: " << duration_fina_eval.count()
			<< " nanoseconds" << endl;

	cout << "y : " << y << std::endl;

	cout << "Generate Partial Proofs  " << endl;
	std::vector<Proof> sigmas;
	mpz_class ris = 0;
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		Proof sigma;
		if (i != (NR_CLIENTS)) {
			mpz_class ri(1);
			ris = ris + ri;
			auto start_partial_proof = high_resolution_clock::now();
			lhs.partial_proof(sk, vk, mpz_class(3), mpz_class(i + 1), i,
					&sigma);
			auto final_partial_proof = high_resolution_clock::now();
			auto duration_proof = duration_cast<microseconds>(
					final_partial_proof - start_partial_proof);

			cout << "Time taken by partial_proof: " << duration_proof.count()
					<< " microseconds" << endl;
			sigmas.push_back(sigma);
		} else {
			auto start_partial_proof = high_resolution_clock::now();
			lhs.partial_proof(sk, vk, mpz_class(3), mpz_class(i + 1), i,
					&sigma);
			auto final_partial_proof = high_resolution_clock::now();
			auto duration_proof = duration_cast<microseconds>(
					final_partial_proof - start_partial_proof);

			cout << "Time taken by partial_proof: " << duration_proof.count()
					<< " microseconds" << endl;
			sigmas.push_back(sigma);

		}
	}

	cout << "Finish Partial Proofs  " << endl;
	Proof final_p;
	auto start_final_proof = high_resolution_clock::now();
	lhs.final_proof(vk, mpz_class(2), sigmas, &final_p);
	auto final_final_proof = high_resolution_clock::now();
	auto duration_final_proof = duration_cast<microseconds>(
			final_final_proof - start_final_proof);

	cout << "Time taken by final_proof: " << duration_final_proof.count()
			<< " microseconds" << endl;

	cout << "Generate Final Proof  " << endl;

	cout << "Finish Final Proof  " << endl;

	auto start_verify = high_resolution_clock::now();
	int result = lhs.verify(vk, final_p, y);
	auto final_verify = high_resolution_clock::now();
	auto duration_verify = duration_cast<microseconds>(
			final_verify - start_verify);

	cout << "Time taken by verify: " << duration_verify.count()
			<< " microseconds" << endl;
	if (result == 1) {
		std::cout << "PERFECT!!! " << std::endl;
		std::cout << "y = " << y << std::endl;
	} else {
		std::cout << "WRONG!!! " << std::endl;
	}
	return 0;
}

