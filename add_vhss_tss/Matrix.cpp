#include "Matrix.h"

Matrix::Matrix(int _rows, int _cols) {
	this->row = _rows;
	this->col = _cols;
	data.resize(_rows);
	for (unsigned i = 0; i < data.size(); i++) {
		data[i].resize(_cols, mpz_class(0));
	}

}

Matrix::Matrix(int _rows, int _cols, bool random) {
	this->row = _rows;
	this->col = _cols;
	gmp_randclass rr(gmp_randinit_default);

	data.resize(_rows);

	for (unsigned i = 0; i < data.size(); i++) {
		data[i].resize(_cols, mpz_class(0));
	}

	if (random) {

		for(int i = 0;  i < row;i++){
			for(int j= 0; j < col;j++){
				mpz_class value = rr.get_z_bits(mpz_class(30));
				data[i][j] = value;

			}
		}

	}

}

Matrix::~Matrix() {

	/*if (col > 0) {
	 // release the memory allocated for each row
	 for (int i = 0; i < row; i++)
	 delete[] data[i];
	 }*/

	/*if (row > 0)
	 delete[] data;*/
}

mpz_class Matrix::getElement(int i, int j) {
	return data[i][j];
}

void Matrix::setElement(int i, int j, mpz_class value) {
	data[i][j] = value;
}

void Matrix::display() {

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++)
			std::cout << data[i][j] << "\t";
		std::cout << std::endl;
	}
	std::cout << "---------------------" << std::endl << std::endl;
}

void Matrix::getCofactor(Matrix *temp, int p, int q, int n) {
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int r = 0; r < n; r++) {
		for (int c = 0; c < n; c++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (r != p && c != q) {
				mpz_class value = this->data[r][c];
				///std::cout<< value << std::endl;
				temp->setElement(i, j++, value);
				//temp[i][j++] = A[row][col];
				//temp->display();

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
 n is current dimension of A[][]. */
mpz_class Matrix::determinant_inner(Matrix A, int n) {
	mpz_class D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1)
		return A.getElement(0, 0);

	Matrix temp(A.row, A.col); // To store cofactors

	mpz_class sign(1);  // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++) {
		// Getting Cofactor of A[0][f]
		A.getCofactor(&temp, 0, f, n);
		//temp.display();
		D += sign * A.data[0][f] * determinant_inner(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

Matrix Matrix::multiply(Matrix &m) {

	Matrix res(row, m.col);

	int i, j, k;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			res.setElement(i, j, 0);
			//res[i][j] = 0;
			for (k = 0; k < m.col; k++) {
				mpz_class result = this->getElement(i, k) * m.getElement(k, j);
				result += res.getElement(i, j);
				res.setElement(i, j, result);
				//res[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
	return res;

}

std::vector<mpz_class> Matrix::multiplyByVector(std::vector<mpz_class> vec) {
	//this->display();
	std::vector<mpz_class> result;
	for (int i = 0; i < row; i++) {
		result.push_back(mpz_class(0));
	}
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			result[i] += (this->data[i][j] * vec[j]);
		}
	}
	return result;
}

Matrix Matrix::submatrix(int startRow, int endRow, int startColumn,
		int endColumn) {
	Matrix m_new(endRow, endColumn);
	int j_temp = startColumn;
	int i_temp = startRow;
	int new_row = 0;
	for (int r = 0; r < endRow; r++) {
		int new_col = 0;
		j_temp = startColumn;
		for (int col = 0; col < endColumn; col++) {
			m_new.setElement(new_row, new_col, data[i_temp][j_temp]);
			//m_new->data[new_row][new_col] = this->data[i_temp][j_temp];
			new_col++;
			j_temp++;
		}
		i_temp++;
		new_row++;
	}
	return m_new;
}

mpz_class Matrix::determinant() {
	return determinant_inner(*this, row);
}

Matrix Matrix::adjoint() {

	Matrix adj(row, col);

	mpz_class sign(1);
	Matrix temp(row, col);

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			// Get cofactor of A[i][j]
			this->getCofactor(&temp, i, j, row);

			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			mpz_class value = (sign) * (determinant_inner(temp, row - 1));
			adj.setElement(j, i, value);

		}
	}
	return adj;

}

