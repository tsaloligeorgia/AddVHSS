#include<bits/stdc++.h>
#include <gmpxx.h>
#include <vector>
#include <map>
#include <chrono>

#include "Matrix.h"
#include "VHSSTSS.h"
#include "params.h"
#include "Client.h"
#include "Server.h"
#include "utils.h"

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

	VHSS_TSS tss;
	mpz_class p = Utils::random_prime(SECURITY);
	cout << "p: " << p << endl;
	mpz_class q = Utils::random_prime(SECURITY);
	while (mpz_cmp(p.get_mpz_t(), q.get_mpz_t()) == 0) {
		q = Utils::random_prime(SECURITY);
	}
	cout << "q: " << q << endl;

	//mpz_class p("115396624046509855500193098450942091340358892908888356557165444351499881695511");
	//mpz_class q("98738332208726129244226887848246886965304787986285813154910268348013862125541");
	mpz_class phi_N((p - 1) * (q - 1));
	mpz_class N(p * q);

	std::vector<Client> clients;
	std::vector<Server> servers;

	std::vector<mpz_class> pub_keys;
	mpz_class phi_field(FINITE_FIELD - 1);
	mpz_class R_is(0);
	std::vector<microseconds> key_gen_timing;
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		mpz_class pk;
		mpz_class sk;
		auto start_key_setup = high_resolution_clock::now();
		tss.key_gen(p, q, &pk, &sk);
		auto final_key_setup = high_resolution_clock::now();
		auto duration_key_setup = duration_cast<microseconds>(
				final_key_setup - start_key_setup);
		key_gen_timing.push_back(duration_key_setup);

		if (i != NR_CLIENTS) {
			mpz_class r_i = Utils::random_element();
			R_is = R_is + r_i;

			Client c(i, input_data[i-1], sk, pk, r_i, tss);
			pub_keys.push_back(pk);
			clients.push_back(c);
		} else {

			mpz_class R_i = (((R_is / phi_field) * phi_field) - R_is) % phi_field;
			//cout << "R_i " << R_i << endl;
			Client c1(i, input_data[i-1], sk, pk, R_i, tss);
			pub_keys.push_back(pk);
			clients.push_back(c1);

		}
	}
	sort(key_gen_timing.begin(), key_gen_timing.end());
	cout << "Time taken by key_gen: " << key_gen_timing[0].count()
			<< " microseconds" << endl;

	for (int i = 1; i < NR_SERVERS + 1; i++) {
		Server s(i);
		servers.push_back(s);

	}

	std::map<int, std::map<int, mpz_class>> omegas;
	std::vector<Matrix> matrix_as;
	std::vector<mpz_class> hash_Hs;

	cout << "Generating Shares " << endl;
	std::vector<microseconds> gen_shares_timing;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		mpz_class hash;
		std::map<int, mpz_class> shared_keys;
		std::map<int, mpz_class> shares;
		Matrix A(NR_SERVERS, THRESHOLD);
		auto start_gen_shares = high_resolution_clock::now();
		c.generate_shares(N, &shares, &shared_keys, &A, &hash);
		auto final_gen_shares = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(
				final_gen_shares - start_gen_shares);
		gen_shares_timing.push_back(duration);

		omegas.insert(
				std::pair<int, std::map<int, mpz_class>>(c.getI(),
						shared_keys));
		matrix_as.push_back(A);
		hash_Hs.push_back(hash);
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
		mpz_class partial = tss.partial_eval(s.getJ(), s.getShares());
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
	mpz_class y = tss.final_eval(partial_evals);
	auto final_final_eval = high_resolution_clock::now();
	auto duration_fina_eval = duration_cast<nanoseconds>(
			final_final_eval - start_final_eval);
	cout << "Finish Final Evals " << endl;
	cout << "Time taken by final_eval: " << duration_fina_eval.count()
			<< " nanoseconds" << endl;

	cout << "Generate Partial Proofs  " << endl;
	auto start_partial_proof = high_resolution_clock::now();
	std::map<int, std::map<int, mpz_class>> r = tss.partial_proof(omegas,
			hash_Hs, matrix_as, N, phi_N);
	auto final_partial_proof = high_resolution_clock::now();
	auto duration_proof = duration_cast<microseconds>(
			final_partial_proof - start_partial_proof);

	cout << "Time taken by partial_proof: " << duration_proof.count()
			<< " microseconds" << endl;
	cout << "Finish Partial Proofs  " << endl;

	cout << "Generate Final Proof  " << endl;
	auto start_final_proof = high_resolution_clock::now();
	mpz_class sigma = tss.final_proof(r, pub_keys, hash_Hs, matrix_as, N,
			phi_N);
	auto final_final_proof = high_resolution_clock::now();
	auto duration_final_proof = duration_cast<microseconds>(
			final_final_proof - start_final_proof);

	cout << "Time taken by final_proof: " << duration_final_proof.count()
			<< " microseconds" << endl;
	cout << "Finish Final Proof  " << endl;

	auto start_verify = high_resolution_clock::now();
	int result = tss.verify(hash_Hs, sigma, y);
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

