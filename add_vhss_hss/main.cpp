#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <map>
#include <fstream>
#include <string>

#include <chrono>
using namespace std::chrono;

#include "Client.h"
#include "Server.h"
#include "VHSSHSS.h"
#include "params.h"

using namespace std;

void print_parameters() {
	std::cout << "-------------------------" << std::endl;
	std::cout << "--- \tParameters\t ---" << std::endl;
	std::cout << "--- \tNR_CLIENTS\t ---\t" << NR_CLIENTS << std::endl;
	std::cout << "--- \tNR_SERVERS\t ---\t" << NR_SERVERS << std::endl;
	std::cout << "--- \tFINITE_FIELD\t ---\t" << FINITE_FIELD << std::endl;
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

	std::vector<Client> clients;
	std::vector<Server> servers;
	mpz_class q(FINITE_FIELD);
	mpz_class one(1);

	VHSS_HSS vhss;

	mpz_class R_is(0);
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		if (i != NR_CLIENTS) {
			mpz_class r_i = Utils::random_element();
			R_is = R_is + r_i;
			Client c(i, input_data[i - 1], r_i, vhss);
			clients.push_back(c);
		} else {
			mpz_class r_i = (R_is / (q - 1)) * (q - 1) - R_is;
			Client c(i, input_data[i - 1], r_i, vhss);
			clients.push_back(c);
		}

	}

	for (int j = 1; j < NR_SERVERS + 1; j++) {
		Server s(j);
		servers.push_back(s);
	}

	std::vector<microseconds> gen_shares_timing;
	std::vector<mpz_class> taus;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		std::map<int, mpz_class> shares;

		mpz_class tau;
		auto start_gen_shares = high_resolution_clock::now();
		c.generate_shares(&shares, &tau);
		auto final_gen_shares = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(
				final_gen_shares - start_gen_shares);
		gen_shares_timing.push_back(duration);

		taus.push_back(tau);
		for (int j = 0; j < NR_SERVERS; j++) {
			Server s = servers[j];
			s.setShare(c.getI(), shares[s.getJ()]);
			servers[j] = s;
		}

	}
	sort(gen_shares_timing.begin(), gen_shares_timing.end());
	cout << "Time taken by generate_shares: " << gen_shares_timing[0].count()
			<< " microseconds" << endl;

	/*for (microseconds duration : gen_shares_timing) {
	 cout << "Time taken by generate_shares: " << duration.count()
	 << " microseconds" << endl;
	 }*/

	std::vector<mpz_class> partial_evals;
	std::vector<mpz_class> partial_proofs;

	std::vector<microseconds> partial_evals_timing;
	std::vector<microseconds> partial_proofs_timing;
	for (Server s : servers) {
		std::map<int, mpz_class> shares = s.getShares();
		auto start_partial_proof = high_resolution_clock::now();
		mpz_class sigma = vhss.partial_proof(shares);

		auto final_partial_proof = high_resolution_clock::now();
		auto duration_proof = duration_cast<microseconds>(
				final_partial_proof - start_partial_proof);

		partial_proofs_timing.push_back(duration_proof);



		auto start_partial_eval = high_resolution_clock::now();
		mpz_class y_j = vhss.partial_eval(s.getJ(), shares);

		auto final_partial_eval = high_resolution_clock::now();
		auto duration_eval = duration_cast<microseconds>(
				final_partial_eval - start_partial_eval);



		partial_evals_timing.push_back(duration_eval);

		partial_evals.push_back(y_j);
		partial_proofs.push_back(sigma);

	}

	sort(partial_evals_timing.begin(), partial_evals_timing.end());
	sort(partial_proofs_timing.begin(), partial_proofs_timing.end());

	cout << "Time taken by partial_proof: " << partial_proofs_timing[NR_SERVERS/2].count()
					<< " microseconds" << endl;

	cout << "Time taken by partial_eval: " << partial_proofs_timing[NR_SERVERS/2].count()
					<< " microseconds" << endl;

	auto start_final_eval = high_resolution_clock::now();
	mpz_class y = vhss.final_eval(partial_evals);
	auto final_final_eval = high_resolution_clock::now();
	auto duration_fina_eval = duration_cast<nanoseconds>(
			final_final_eval - start_final_eval);

	cout << "Time taken by final_eval: " << duration_fina_eval.count()
			<< " nanoseconds" << endl;

	auto start_final_proof = high_resolution_clock::now();
	mpz_class sigma = vhss.final_proof(partial_proofs);
	auto final_final_proof = high_resolution_clock::now();
	auto duration_final_proof = duration_cast<nanoseconds>(
			final_final_proof - start_final_proof);

	cout << "Time taken by final_proof: " << duration_final_proof.count()
			<< " nanoseconds" << endl;

	auto start_verify = high_resolution_clock::now();
	int result = vhss.verify(taus, sigma, y);
	auto final_verify = high_resolution_clock::now();
	auto duration_verify = duration_cast<microseconds>(
			final_verify - start_verify);

	cout << "Time taken by verify: " << duration_verify.count()
			<< " microseconds" << endl;

	if (result == 1) {
		std::cout << "y = " << (y % mpz_class(FINITE_FIELD)) << std::endl;
	} else {
		std::cout << "fail!" << std::endl;
	}

	return 0;
}
