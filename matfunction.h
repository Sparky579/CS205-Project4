#pragma once
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>
typedef struct Matrix Matrix;
bool sameSize(const Matrix *mat1, const Matrix *mat2);
float * pos(const Matrix *mat, size_t r, size_t c);
Matrix * createMatrix(size_t rows, size_t columns);
bool deleteMatrix(Matrix **mat);
bool addPosition(Matrix *mat, size_t r, size_t c, float num);
bool substractPosition(Matrix *mat, size_t r, size_t c, float num);
bool setPosition(Matrix *mat, size_t r, size_t c, float num);
bool addScalar(Matrix *mat, float num);
bool substrScalar(Matrix *mat, float num);
bool multiScalar(Matrix *mat, float num);
float minimalValue(const Matrix *mat);
float maximalValue(const Matrix *mat);
bool addMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool substractMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool copyMatrix(const Matrix *mat1, Matrix *mat2);
bool multiMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool printMatrix(const Matrix *mat);
bool matmul_plain(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved2(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved2_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved3(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved3_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_2(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_3(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat, int BLOCK_SIZE);
bool matmul_improved4_X_1_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X_2_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X_4_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X_8_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X_16_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
bool matmul_improved4_X_32_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat);
