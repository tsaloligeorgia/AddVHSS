/*



 * test.cpp
 *
 *  Created on: Jul 23, 2020
 *      Author: yoda
*/

#include<bits/stdc++.h>
#include <gmpxx.h>

#include "Matrix.h"
#include "VHSSTSS.h"
#include "params.h"
#include "Client.h"
#include "Server.h"
#include "utils.h"

using namespace std;

// Driver program
int main() {

	VHSS_TSS tss;
	mpz_class p = Utils::random_prime(128);
	cout << "p: " << p << endl;
	mpz_class q = Utils::random_prime(128);
	while (mpz_cmp(p.get_mpz_t(), q.get_mpz_t()) == 0) {
		q = Utils::random_prime(128);
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
	for (int i = 1; i < NR_CLIENTS + 1; i++) {
		mpz_class pk;
		mpz_class sk;
		tss.key_gen(p, q, &pk, &sk);

		if (i != NR_CLIENTS) {
			mpz_class R_i = Utils::generate_random(30, phi_field);
			//cout << "R_i " << R_i << endl;
			R_is = R_is + R_i;

			Client c(i, mpz_class(2), sk, pk, R_i, tss);
			pub_keys.push_back(pk);
			clients.push_back(c);
		} else {
			mpz_class R_i = (((R_is / phi_field) * phi_field) - R_is) % phi_field;
			//cout << "R_i " << R_i << endl;
			Client c1(i, mpz_class(2), sk, pk, R_i, tss);
			pub_keys.push_back(pk);
			clients.push_back(c1);

		}
	}

	for (int i = 1; i < NR_SERVERS + 1; i++) {
		Server s(i);
		servers.push_back(s);

	}


	std::map<int, std::map<int, mpz_class>> omegas;
	std::vector<Matrix> matrix_as;
	std::vector<mpz_class> hash_Hs;

	cout << "Generating Shares " << endl;
	for (int i = 0; i < NR_CLIENTS; i++) {
		Client c = clients[i];
		mpz_class hash;
		std::map<int, mpz_class> shared_keys;
		std::map<int, mpz_class> shares;
		Matrix A(NR_SERVERS, THRESHOLD);

		c.generate_shares(N, &shares, &shared_keys, &A, &hash);

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

	cout << "Generate Partial Evals " << endl;

	std::map<int, mpz_class> partial_evals;
	for (int j = 0; j < NR_SERVERS; j++) {
		Server s = servers[j];
		mpz_class partial = tss.partial_eval(s.getJ(), s.getShares());

		partial_evals.insert(std::pair<int, mpz_class>(s.getJ(), partial));
	}

	cout << "Finish Partial Evals " << endl;

	cout << "Generate Final Evals " << endl;
	mpz_class y = tss.final_eval(partial_evals);
	cout << "Finish Final Evals " << endl;

	cout << "Generate Partial Proofs  " << endl;
	std::map<int, std::map<int, mpz_class>> r = tss.partial_proof(omegas,
			hash_Hs, matrix_as, N, phi_N);
	cout << "Finish Partial Proofs  " << endl;

	cout << "Generate Final Proof  " << endl;
	mpz_class sigma = tss.final_proof(r, pub_keys, hash_Hs, matrix_as, N,
			phi_N);
	cout << "Finish Final Proof  " << endl;

	int result = tss.verify(hash_Hs, sigma, y);
	if (result == 1) {
		std::cout << "PERFECT!!! " << std::endl;
		std::cout << "y = " << y << std::endl;
	} else {
		std::cout << "WRONG!!! " << std::endl;
	}

	return 0;
}


