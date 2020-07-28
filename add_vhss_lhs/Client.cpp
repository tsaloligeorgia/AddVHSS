#include "Client.h"

Client::Client(int i, mpz_class secret_input, mpz_class random_e, VHSS_LHS vhss) {
	this->i = i;
	this->secret_input = secret_input;
	this->vhss = vhss;
	this->random_e = random_e;
}

Client::~Client() {
	// TODO Auto-generated destructor stub
}

void Client::generate_shares(std::map<int, mpz_class> *shares) {

	this->vhss.gen_shares_threshold(this->secret_input,this->random_e, shares);

}

int Client::getI(){
	return this->i;
}

mpz_class Client::getSecretInput(){
	return this->secret_input + random_e;
}
