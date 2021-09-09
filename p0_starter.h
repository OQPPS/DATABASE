//===----------------------------------------------------------------------===//
//
//
//
// p0_starter.h
//
// Identification: src/include/primer/p0_starter.h
//
//
//
//===----------------------------------------------------------------------===//

#pragma once

#include <memory>
#include <stdexcept>
#include <vector>

namespace scudb {

/**
 * The Matrix type defines a common
 * interface for matrix operations.
 */
template <typename T>
class Matrix {
 protected:
  /**
   * TODO(P0): Add implementation
   * 
   * Construct a new Matrix instance.
   * @param rows The number of rows
   * @param cols The number of columns
   *
   */
	 T **p;
  Matrix(int rows, int cols) {
	  rows_ = rows;
	  cols_ = cols;
	  for (int i = 0; i < rows_; i++){
		  for (int j = 0; j < cols_; j++){
			  p[i][j] = 0;
		  }
	  }
  }

  /** The number of rows in the matrix */
  int rows_;
  /** The number of columns in the matrix */
  int cols_;

  /**
   * TODO(P0): Allocate the array in the constructor.
   * TODO(P0): Deallocate the array in the destructor.
   * A flattened array containing the elements of the matrix.
   */
  T *linear_;

 public:
  /** @return The number of rows in the matrix */
  virtual int GetRowCount() const = 0;

  /** @return The number of columns in the matrix */
  virtual int GetColumnCount() const = 0;

  /**
   * Get the (i,j)th matrix element.
   *
   * Throw OUT_OF_RANGE if either index is out of range.
   *
   * @param i The row index
   * @param j The column index
   * @return The (i,j)th matrix element
   * @throws OUT_OF_RANGE if either index is out of range
   */
  virtual T GetElement(int i, int j) const = 0;

  /**
   * Set the (i,j)th matrix element.
   *
   * Throw OUT_OF_RANGE if either index is out of range.
   *
   * @param i The row index
   * @param j The column index
   * @param val The value to insert
   * @throws OUT_OF_RANGE if either index is out of range
   */
  virtual void SetElement(int i, int j, T val) = 0;

  /**
   * Fill the elements of the matrix from `source`.
   *
   * Throw OUT_OF_RANGE in the event that `source`
   * does not contain the required number of elements.
   *
   * @param source The source container
   * @throws OUT_OF_RANGE if `source` is incorrect size
   */
  virtual void FillFrom(const std::vector<T> &source) = 0;

  /**
   * Destroy a matrix instance.
   * TODO(P0): Add implementation
   */
  virtual ~Matrix() = default;
};

/**
 * The RowMatrix type is a concrete matrix implementation.
 * It implements the interface defined by the Matrix type.
 */
template <typename T>
class RowMatrix : public Matrix<T> {
 public:
  /**
   * TODO(P0): Add implementation
   *
   * Construct a new RowMatrix instance.
   * @param rows The number of rows
   * @param cols The number of columns
   */
  RowMatrix(int rows, int cols) : Matrix<T>(rows, cols) {
	  rows_ = rows;
	  cols_ = cols;
	  for (int i = 0; i < rows_; i++){
		  for (int j = 0; j < cols_; j++){
			  p[i][j] = 0;
		  }
	  }
  }

  /**
   * TODO(P0): Add implementation
   * @return The number of rows in the matrix
   */
  int GetRowCount() const override { return this->rows_; }

  /**
   * TODO(P0): Add implementation
   * @return The number of columns in the matrix
   */
  int GetColumnCount() const override { return this->cols_; }

  /**
   * TODO(P0): Add implementation
   *
   * Get the (i,j)th matrix element.
   *
   * Throw OUT_OF_RANGE if either index is out of range.
   *
   * @param i The row index
   * @param j The column index
   * @return The (i,j)th matrix element
   * @throws OUT_OF_RANGE if either index is out of range
   */
  T GetElement(int i, int j) const override {
	  if (i >= rows_ || i<0 || j<0 || j >= cols_) 
		  throw NotImplementedException{ "RowMatrix::GetElement() not implemented." }£»
		  return p[i][j]; 
  }

  /**
   * Set the (i,j)th matrix element.
   *
   * Throw OUT_OF_RANGE if either index is out of range.
   *
   * @param i The row index
   * @param j The column index
   * @param val The value to insert
   * @throws OUT_OF_RANGE if either index is out of range
   */
  void SetElement(int i, int j, T val) override
  {
	  if (i >= rows_ || i < 0 || j < 0 || j >= cols_)
	  {
		  throw NotImplementedException{ "RowMatrix::SetElement() not implemented." };
	  }
	  else{ p[i][j] = val; }
  }
  /**
   * TODO(P0): Add implementation
   *
   * Fill the elements of the matrix from `source`.
   *
   * Throw OUT_OF_RANGE in the event that `source`
   * does not contain the required number of elements.
   *
   * @param source The source container
   * @throws OUT_OF_RANGE if `source` is incorrect size
   */
  void FillFrom(const std::vector<T> &source) override
  {
	  if (source.size()<i*k) throw NotImplementedException{ "RowMatrix::FillFrom() not implemented." };
	  for (int i = 0; i < rows; i++)
	  {
		  for (int j = 0; j < cols; j++)
		  {
			  p[i][j] = source[i * cols + j];
		  }
	  }
  }

  /**
   * TODO(P0): Add implementation
   *
   * Destroy a RowMatrix instance.
   */
  ~RowMatrix() override
  {
	  for (int i = 0; i < rows; i++)  delete[]p[i];
  }

 private:
  /**
   * A 2D array containing the elements of the matrix in row-major format.
   *
   * TODO(P0):
   * - Allocate the array of row pointers in the constructor.
   * - Use these pointers to point to corresponding elements of the `linear` array.
   * - Don't forget to deallocate the array in the destructor.
   */
  T **data_;
};

/**
 * The RowMatrixOperations class defines operations
 * that may be performed on instances of `RowMatrix`.
 */
template <typename T>
class RowMatrixOperations {
 public:
  /**
   * Compute (`matrixA` + `matrixB`) and return the result.
   * Return `nullptr` if dimensions mismatch for input matrices.
   * @param matrixA Input matrix
   * @param matrixB Input matrix
   * @return The result of matrix addition
   */
    // TODO(P0): Add implementation
	 static unique_ptr<RowMatrix<T>> Add(const RowMatrix<T> *matrixA, const RowMatrix<T> *matrixB)
	  {
		  if (matrixA.GetColumnCount() != matrixB.GetColumnCount() || matrixA.GetRowCount() != matrixB.GetRowCount())  return unique_ptr<RowMatrix<T>>(nullptr);
		  unique_ptr<RowMatrix<T>> matrixC(new RowMatrix<T>(matrixA.GetRowCount(), matrixA.GetColumnCount()));
		  for (int i = 0; i< MatrixA.GetRowCount(); i++)
		  {
			  for (int j = 0; j< MatrixA.GetColumnCount(); j++)
			  {
				  matrixC->p[i][j] = matrixA->p[i][j] = matrixB->p[i][j];
			  }
		  }
		  return matrixC;
	  }

  /**
   * Compute the matrix multiplication (`matrixA` * `matrixB` and return the result.
   * Return `nullptr` if dimensions mismatch for input matrices.
   * @param matrixA Input matrix
   * @param matrixB Input matrix
   * @return The result of matrix multiplication
   */
    // TODO(P0): Add implementation
	 static std::unique_ptr<RowMatrix<T>> Multiply(const RowMatrix<T> *matrixA, const RowMatrix<T> *matrixB)
	  {
		  if (matrixA.GetColumnCount() != matrixB.GetRowCount())  return std::unique_ptr<RowMatrix<T>>(nullptr);
		  unique_ptr<RowMatrix<T>> matrixC(new RowMatrix<T>(matrixA.GetRowCount(), matrixA.GetColumnCount()));
		  for (int i = 0; i< Matrix.GetRowCount(); i++)
		  {
			  for (int j = 0; j< Matrix.GetColumnCount(); j++)
			  {
				  matrixC->a[i][j] = 0;
				  for (int k = 0; k< Matrix.GetColumnCount(); k++)
					  matrixC->a[i][j] += matrixA->a[i][k] * matrixB->a[k][j];
			  }
		  }
		  return matrixC;
	  }

  /**
   * Simplified General Matrix Multiply operation. Compute (`matrixA` * `matrixB` + `matrixC`).
   * Return `nullptr` if dimensions mismatch for input matrices.
   * @param matrixA Input matrix
   * @param matrixB Input matrix
   * @param matrixC Input matrix
   * @return The result of general matrix multiply
   */
    // TODO(P0): Add implementation
	 static std::unique_ptr<RowMatrix<T>> GEMM(const RowMatrix<T> *matrixA, const RowMatrix<T> *matrixB,
		 const RowMatrix<T> *matrixC)
	 {
		 if (matrixA.GetColumnCount() != matrixB.GetRowCount())  return unique_ptr<RowMatrix<T>>(nullptr);
		 if (matrixA.GetRowCount() != matrixC.GetRowCount() || matrixB.GetColumnCount() != matrixC.GetColumnCount())  return unique_ptr<RowMatrix<T>>(nullptr);
		 unique_ptr<RowMatrix<T>> matrixD(new RowMatrix<T>(matrixA.GetRowCount(), matrixA.GetColumnCount()));
		 for (int i = 0; i< Matrix.GetRowCount(); i++)
		 {
			 for (int j = 0; j< Matrix.GetColumnCount(); j++)
			 {
				 matrixD->p[i][j] = 0;
				 for (int k = 0; k< Matrix.GetColumnCount(); k++)
					 matrixD->p[i][j] += matrixA->p[i][k] * matrixB->p[k][j];
			 }
		 }
		 for (int i = 0; i< MatrixC.GetRowCount(); i++)
		 {
			 for (int j = 0; j< MatrixC.GetColumnCount(); j++)
			 {
				 matrixD->p[i][j] += matrixC->p[i][j];
			 }
		 }
		 return matrixD;
	 }
};
}  // namespace scudb
