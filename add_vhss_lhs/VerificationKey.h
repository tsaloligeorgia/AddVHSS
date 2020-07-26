/*
 * VerificationKey.h
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#ifndef VERIFICATIONKEY_H_
#define VERIFICATIONKEY_H_
#include <vector>
#include <gmpxx.h>
#include <iostream>

class VerificationKey {
public:
	VerificationKey();
	virtual ~VerificationKey();

	void setN(mpz_class N);

	void setNHat(mpz_class N_hat);

	void setG(mpz_class g);

	void setG1(mpz_class g1);

	void setHs(std::vector<mpz_class> Hs);

	mpz_class getN();

	mpz_class getNHat();

	mpz_class getG();

	mpz_class getG1();

	std::vector<mpz_class> getHS();

private:
	mpz_class N;
	mpz_class N_hat;
	mpz_class g;
	mpz_class g_1;
	std::vector<mpz_class> Hs;
};

#endif /* VERIFICATIONKEY_H_ */
