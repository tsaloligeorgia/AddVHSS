#ifndef CLIENT_H_
#define CLIENT_H_

#include <gmpxx.h>

#include "VHSSLHS.h"

class Client {
public:
	Client(int i, mpz_class secret_input, mpz_class random_e, VHSS_LHS vhss);
	virtual ~Client();

	void generate_shares(std::map<int, mpz_class> *shares);

	int getI();

	mpz_class getSecretInput();
private:
	int i;
	mpz_class secret_input;
	VHSS_LHS vhss;
	mpz_class random_e;
};

#endif /* CLIENT_H_ */
