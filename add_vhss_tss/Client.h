#ifndef CLIENT_H_
#define CLIENT_H_

#include <gmpxx.h>

#include "VHSSTSS.h"

class Client {
public:
	Client(int i, mpz_class secret_input, mpz_class sk, mpz_class pk,
			mpz_class random_e, VHSS_TSS tss);
	virtual ~Client();

	void generate_shares(mpz_class N, std::map<int, mpz_class> *shares,
			std::map<int, mpz_class> *shared_key, Matrix *A_i, mpz_class *hash_h);

	int getI();
private:
	int i;
	mpz_class secret_input;
	mpz_class sk;
	mpz_class pk;
	VHSS_TSS tss;
	mpz_class random_e;
};

#endif /* CLIENT_H_ */
