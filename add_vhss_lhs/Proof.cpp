/*
 * Proof.cpp
 *
 *  Created on: Jul 24, 2020
 *      Author: yoda
 */

#include "Proof.h"

Proof::Proof() {
	// TODO Auto-generated constructor stub

}

Proof::~Proof() {
	// TODO Auto-generated destructor stub
}

void Proof::setE(mpz_class e) {
	this->e = e;

}
void Proof::setSi(mpz_class si) {
	this->si = si;

}
void Proof::setX_i(mpz_class x_i_tilde) {
	this->x_i_tilde = x_i_tilde;

}

mpz_class Proof::getE(){
	return this->e;
}

mpz_class Proof::getSi(){
	return this->si;
}


mpz_class Proof::getX(){
	return this->x_i_tilde;
}


