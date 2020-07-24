
#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <gmpxx.h>
#include <vector>

class Matrix {
public:
	Matrix(int row, int col);
	Matrix(int row, int col, bool random);
	virtual ~Matrix();
	mpz_class getElement(int i, int j);
	void setElement(int i, int j, mpz_class value);
	void display();

	mpz_class determinant();
	Matrix adjoint();

	Matrix multiply(Matrix &m);

	std::vector<mpz_class> multiplyByVector(std::vector<mpz_class> vec);

	Matrix submatrix(int startRow, int endRow, int startColumn, int endColumn );




private:
	int row;
	int col;
	std::vector<std::vector<mpz_class> > data;

	void getCofactor(Matrix *temp, int p, int q, int n);
	mpz_class determinant_inner(Matrix A, int n);

};

#endif /* MATRIX_H_ */
