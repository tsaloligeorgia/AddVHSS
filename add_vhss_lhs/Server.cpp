#include "Server.h"

Server::Server(int j) {
	this->j = j;

}

Server::~Server() {
	// TODO Auto-generated destructor stub
}

int Server::getJ() {
	return this->j;
}

void Server::printShares() {
	std::cout << "PRINTING SHARES " << std::endl;

	std::cout << "WHile SHARES " << std::endl;
	std::map<int, mpz_class>::iterator it = this->shares.begin();
	// Iterate over the map using Iterator till end.
	while (it != this->shares.end()) {
		std::cout << " i : " << it->first << "\t";
		std::cout << "share: " << it->second << std::endl;
		it++;
	}

}

void Server::setShare(int i, mpz_class share) {
	this->shares.insert(std::pair<int, mpz_class>(i, share));

}

std::map<int, mpz_class> Server::getShares() {

	return this->shares;
}
