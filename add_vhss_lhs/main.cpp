#include<bits/stdc++.h>
#include <gmpxx.h>
#include <vector>
#include <map>
#include <chrono>
#include <string>
#include <fstream>

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
		std::cout << vk->getHS().at(i) << std::endl;
	}
	std::cout << "---- --- ----" << std::endl;

}

void read_file(std::vector<mpz_class> *to_save) {
	std::vector<mpz_class> to_save_t;

	ifstream myfile("result.txt");
	if (myfile.is_open()) {
		string line;
		getline(myfile, line);
		for (int i = 0; i < NR_CLIENTS + 10; i++) {
			getline(myfile, line);
			std::string::size_type sz;     // alias of size_t
			float converted = std::stof(line, &sz);
			mpz_class changed_value(converted * 100);
			to_save_t.push_back(changed_value);

		}
		myfile.close();
	}
	*to_save = to_save_t;

}
int main() {
	print_parameters();

	std::vector<mpz_class> input_data;

	read_file(&input_data);

	VHSS_LHS lhs;
	mpz_class p = Utils::generate_safe_prime(SECURITY, mpz_class(2)); //Utils::generate_safe_prime(32, mpz_class(2));
	//cout << "p: " << p << endl;
	mpz_class q = Utils::generate_safe_prime(SECURITY, p);
	while (mpz_cmp(p.get_mpz_t(), q.get_mpz_t()) == 0) {
		q = Utils::generate_safe_prime(SECURITY, p);
	}
	//cout << "q: " << q << endl;

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
	cout << "Time taken by key_gen: " << duration_key_setup.count()
			<< " microseconds" << endl;

	//print_secret(&sk);
	//print_verification(&vk);
	//print_verification(&vk);
	mpz_class R_is(0);

	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		if (i != (NR_CLIENTS)) {
			mpz_class r_i = Utils::random_element();
			R_is = R_is + r_i;
			Client c(i, input_data[i - 1], r_i, lhs);
			clients.push_back(c);
		} else {
			mpz_class r_i = (R_is / (q - 1)) * (q - 1) - R_is;
			Client c(i, input_data[i - 1], r_i, lhs);
			clients.push_back(c);
		}
	}

	for (int j = 1; j < NR_SERVERS + 1; j++) {
		Server s(j);
		servers.push_back(s);
	}

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

	sort(gen_shares_timing.begin(), gen_shares_timing.end());
	cout << "Time taken by generate_shares: "
			<< gen_shares_timing[NR_CLIENTS / 2].count() << " microseconds"
			<< endl;

	std::vector<microseconds> partial_eval_timing;
	std::map<int, mpz_class> partial_evals;

	for (int j = 0; j < NR_SERVERS; j++) {
		Server s = servers[j];
		auto start_partial_eval = high_resolution_clock::now();
		mpz_class partial = lhs.partial_eval(s.getJ(), s.getShares());
		auto final_partial_eval = high_resolution_clock::now();
		auto duration_eval = duration_cast<microseconds>(
				final_partial_eval - start_partial_eval);

		partial_eval_timing.push_back(duration_eval);

		partial_evals.insert(std::pair<int, mpz_class>(s.getJ(), partial));
	}

	sort(partial_eval_timing.begin(), partial_eval_timing.end());

	cout << "Time taken by partial_eval: "
			<< partial_eval_timing[NR_SERVERS / 2].count() << " microseconds"
			<< endl;

	auto start_final_eval = high_resolution_clock::now();
	mpz_class y = lhs.final_eval(partial_evals);
	auto final_final_eval = high_resolution_clock::now();
	auto duration_fina_eval = duration_cast<nanoseconds>(
			final_final_eval - start_final_eval);
	cout << "Time taken by final_eval: " << duration_fina_eval.count()
			<< " nanoseconds" << endl;

	std::vector<microseconds> timing_partial_proofs;
	std::vector<Proof> sigmas;
	mpz_class finite(FINITE_FIELD);
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		Proof sigma;
		Client c = clients[i - 1];
		if (i != (NR_CLIENTS)) {
			mpz_class r_i = Utils::random_element();
			R_is = R_is + r_i;
			auto start_partial_proof = high_resolution_clock::now();
			lhs.partial_proof(&sk, &vk, mpz_class(2), c.getSecretInput(), i,
					&sigma);
			auto final_partial_proof = high_resolution_clock::now();
			auto duration_proof = duration_cast<microseconds>(
					final_partial_proof - start_partial_proof);

			timing_partial_proofs.push_back(duration_proof);

	/*		 cout << "Time taken by partial_proof: " << duration_proof.count()
			 << " microseconds" << endl;*/
			sigmas.push_back(sigma);
		} else {
			mpz_class r_i = (R_is / (finite - 1)) * (finite - 1) - R_is;
			auto start_partial_proof = high_resolution_clock::now();
			lhs.partial_proof(&sk, &vk, mpz_class(2), c.getSecretInput(), i,
					&sigma);
			auto final_partial_proof = high_resolution_clock::now();
			auto duration_proof = duration_cast<microseconds>(
					final_partial_proof - start_partial_proof);
/*
						cout << "Time taken by partial_proof: " << duration_proof.count()
			 << " microseconds" << endl;*/
			sigmas.push_back(sigma);
			timing_partial_proofs.push_back(duration_proof);

		}
	}

	//print_verification(&vk);

	sort(timing_partial_proofs.begin(), timing_partial_proofs.end());
	cout << "Time taken by partial_proof: "
			<< timing_partial_proofs[NR_CLIENTS / 2].count() << " microseconds"
			<< endl;

	Proof final_p;
	auto start_final_proof = high_resolution_clock::now();
	lhs.final_proof(&vk, mpz_class(2), sigmas, &final_p);
	auto final_final_proof = high_resolution_clock::now();
	auto duration_final_proof = duration_cast<microseconds>(
			final_final_proof - start_final_proof);

	cout << "Time taken by final_proof: " << duration_final_proof.count()
			<< " microseconds" << endl;

	auto start_verify = high_resolution_clock::now();
	int result = lhs.verify(vk, final_p, y);
	auto final_verify = high_resolution_clock::now();
	auto duration_verify = duration_cast<microseconds>(
			final_verify - start_verify);

	cout << "Time taken by verify: " << duration_verify.count()
			<< " microseconds" << endl;
	if (result == 1) {
		std::cout << "y = " << y << std::endl;
	} else {
		std::cout << "WRONG!!! " << std::endl;
	}
	return 0;
}

