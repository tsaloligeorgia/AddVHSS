#ifndef SERVER_H_
#define SERVER_H_

#include <gmpxx.h>
#include <map>
#include <iostream>

#include "params.h"

class Server {
public:
	Server(int j);
	virtual ~Server();

	int getJ();

	void setShare(int i, mpz_class share);
	std::map<int, mpz_class> getShares();

	void printShares();

private:
	int j;
	std::map<int, mpz_class> shares;
};

#endif /* SERVER_H_ */
