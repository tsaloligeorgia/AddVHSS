/*
 * VerificationKey.cpp
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#include "VerificationKey.h"

VerificationKey::VerificationKey() {
	// TODO Auto-generated constructor stub

}

VerificationKey::~VerificationKey() {
	// TODO Auto-generated destructor stub
}

void VerificationKey::setN(mpz_class N) {
	this->N = N;

}

void VerificationKey::setNHat(mpz_class N_hat) {
	this->N_hat = N_hat;
}

void VerificationKey::setG(mpz_class g) {
	this->g = g;
}

void VerificationKey::setG1(mpz_class g1) {
	this->g_1 = g1;
}

void VerificationKey::setHs(std::vector<mpz_class> Hs) {
	this->Hs = Hs;
}

mpz_class VerificationKey::getN() {
	return this->N;
}

mpz_class VerificationKey::getNHat() {
	return this->N_hat;
}


mpz_class VerificationKey::getG() {
	return this->g;
}

mpz_class VerificationKey::getG1() {
	return this->g_1;
}

std::vector<mpz_class> VerificationKey::getHS() {
	return this->Hs;
}
