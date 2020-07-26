/*
 * Proof.h
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#ifndef PROOF_H_
#define PROOF_H_

#include <gmpxx.h>
#include <iostream>


class Proof {
public:
	Proof();
	virtual ~Proof();

	void setE(mpz_class e);
	void setSi(mpz_class si);
	void setX_i(mpz_class x_i_tilde);

	mpz_class getE();

	mpz_class getSi();

	mpz_class getX();

private:
	mpz_class e;
	mpz_class si;
	mpz_class x_i_tilde;
};

#endif /* PROOF_H_ */
