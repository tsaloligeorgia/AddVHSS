/*
 * SecretKey.cpp
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#include "SecretKey.h"

SecretKey::SecretKey() {
	// TODO Auto-generated constructor stub

}

SecretKey::~SecretKey() {
	// TODO Auto-generated destructor stub
}

void SecretKey::setPHat(mpz_class newP) {
	this->p_hat = newP;

}

void SecretKey::setQHat(mpz_class newQ) {
	this->q_hat = newQ;

}

mpz_class SecretKey::getP() {

	return this->p_hat;

}

mpz_class SecretKey::getQ() {

	return this->q_hat;

}

