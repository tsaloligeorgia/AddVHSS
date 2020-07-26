#include "Client.h"

Client::Client(int i, mpz_class secret_input, mpz_class random_e, VHSS_HSS vhss) {
	this->i = i;
	this->secret_input = secret_input;
	this->vhss = vhss;
	this->random_e = random_e;
}

Client::~Client() {
	// TODO Auto-generated destructor stub
}

void Client::generate_shares(std::map<int, mpz_class> *shares, mpz_class *tau) {

	this->vhss.gen_shares(this->secret_input,this->random_e, shares, tau);

}

int Client::getI(){
	return this->i;
}
