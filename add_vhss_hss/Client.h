#ifndef CLIENT_H_
#define CLIENT_H_

#include <gmpxx.h>

#include "VHSSHSS.h"

class Client {
public:
	Client(int i, mpz_class secret_input, mpz_class random_e, VHSS_HSS vhss);
	virtual ~Client();

	void generate_shares(std::map<int, mpz_class> *shares, mpz_class *tau);

	int getI();
private:
	int i;
	mpz_class secret_input;
	VHSS_HSS vhss;
	mpz_class random_e;
	mpz_class tau;
};

#endif /* CLIENT_H_ */
