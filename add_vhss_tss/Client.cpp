#include "Client.h"

Client::Client(int i, mpz_class secret_input, mpz_class sk, mpz_class pk,
		mpz_class random_e, VHSS_TSS tss) {
	this->i = i;
	this->secret_input = secret_input;
	this->sk = sk;
	this->pk = pk;
	this->tss = tss;
	this->random_e = random_e;
}

Client::~Client() {
	// TODO Auto-generated destructor stub
}

void Client::generate_shares(mpz_class N, std::map<int, mpz_class> *shares,
		std::map<int, mpz_class> *shared_key, Matrix *A_i, mpz_class *hash_h) {

	this->tss.gen_shares_threshold(this->secret_input, this->sk, this->random_e,
			this->pk, shares, shared_key, A_i, hash_h);

}

int Client::getI(){
	return this->i;
}
