/*
 * SecretKey.h
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#ifndef SECRETKEY_H_
#define SECRETKEY_H_

#include <gmpxx.h>
#include <iostream>


class SecretKey {
public:
	SecretKey();
	virtual ~SecretKey();

	void setPHat(mpz_class newqHat);
	void setQHat(mpz_class newQHat);

	mpz_class getP();
	mpz_class getQ();
private:
	mpz_class p_hat;
	mpz_class q_hat;
};

#endif /* SECRETKEY_H_ */
