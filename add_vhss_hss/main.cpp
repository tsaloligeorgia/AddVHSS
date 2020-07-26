#include <iostream>
#include <gmpxx.h>
#include <vector>
#include <map>

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

int main() {
	print_parameters();

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
			Client c(i, mpz_class(2), r_i, vhss);
			clients.push_back(c);
		} else {
			mpz_class r_i = (R_is / (q - 1)) * (q - 1) - R_is;
			Client c(i, mpz_class(2), r_i, vhss);
			clients.push_back(c);
		}

	}

	for (int j = 1; j < NR_SERVERS + 1; j++) {
		Server s(j);
		servers.push_back(s);
	}

	cout << "Generating Shares " << endl;
	std::vector<mpz_class> taus;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		std::map<int, mpz_class> shares;

		mpz_class tau;
		c.generate_shares(&shares, &tau);
		taus.push_back(tau);
		for (int j = 0; j < NR_SERVERS; j++) {
			Server s = servers[j];
			s.setShare(c.getI(), shares[s.getJ()]);
			servers[j] = s;
		}

	}
	cout << "Finished Shares " << endl;

	std::cout << "starting paritals" << std::endl;
	std::vector<mpz_class> partial_evals;
	std::vector<mpz_class> partial_proofs;
	for (Server s : servers) {
		std::map<int, mpz_class> shares = s.getShares();
		mpz_class sigma = vhss.partial_proof(shares);
		mpz_class y_j = vhss.partial_eval(s.getJ(), shares);

		partial_evals.push_back(y_j);
		partial_proofs.push_back(sigma);

	}
	std::cout << "fnished paritals" << std::endl;

	std::cout << "starting finals" << std::endl;
	mpz_class y = vhss.final_eval(partial_evals);
	mpz_class sigma = vhss.final_proof(partial_proofs);
	std::cout << "end finals" << std::endl;

	std::cout << "starting verify" << std::endl;
	int result = vhss.verify(taus, sigma, y);
	if (result == 1) {
		std::cout << "PERFECT!" << std::endl;
		std::cout << "y = " << (y % mpz_class(FINITE_FIELD)) << std::endl;
	} else {
		std::cout << "fail!" << std::endl;
	}
	std::cout << "end verify" << std::endl;

	return 0;
}
