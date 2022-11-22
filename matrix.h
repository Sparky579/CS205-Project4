#pragma once
#include "matfunction.h"
#include <immintrin.h>
#include<stdlib.h>
#include<omp.h>
typedef struct Matrix 
{
	size_t rows, columns;
	float *datas;
} Matrix;
const int MatrixSize = sizeof(Matrix);
bool sameSize(const Matrix *mat1, const Matrix *mat2)
{
    /*to test if mat1 and mat2 are same*/
    if (mat1 == NULL || mat2 == NULL) return 0;
	return mat1->rows == mat2->rows && mat1->columns == mat2->columns;
}
float * pos(const Matrix *mat, size_t r, size_t c)
{
    //return the element of mat[r][c]
    //DO NOT put a NULL matrix or it will return NULL
    if (mat == NULL) return NULL;
	if (r >= mat->rows || c >= mat->columns)
		return NULL;
	return mat->datas + (r * mat->columns + c);
}
Matrix * createMatrix(size_t rows, size_t columns)
{
	//create a matrix of r * c whose elements are all 0
	if (rows <= 0 || columns <= 0) return NULL;
	Matrix* mx = (Matrix*) malloc(MatrixSize);
	mx->datas = (float*)malloc(4 * rows * columns);
	memset(mx->datas, 0, 4 * rows * columns);
	mx->rows = rows;
	mx->columns = columns;
	return mx;
}
bool deleteMatrix(Matrix **mat)
{
    /*delete the matrix and free the memory
    return 1 if delete successfully*/
	if (*mat == NULL) return 0;
	if ((*mat)->datas != NULL) free((*mat)->datas);
	free(*mat);
	*mat = NULL;
    return 1;
}
/*
@deprecated
bool deleteMatrix(Matrix *mat)
{   
    if (mat == NULL) return 0;
	//delete a matrix
	//this version requires to set mat=NULL after using
	if (mat == NULL) return;
	if (mat->datas != NULL) free(mat->datas);
	free(mat);
    return 1;
}
*/
bool setPosition(Matrix *mat, size_t r, size_t c, float num)
{
    //set a scalar num to mat[l][r]
    if (mat == NULL) return 0;
	if (mat->rows <= r || mat->columns <= c) return 0;
	*pos(mat, r, c) = num;
	return 1;
}
bool addPosition(Matrix *mat, size_t r, size_t c, float num)
{
	//add a scalar num to mat[l][r]
    if (mat == NULL) return 0;
	if (mat->rows <= r || mat->columns <= c) return 0;
	*pos(mat, r, c) += num;
	return 1;
}
bool substrPosition(Matrix *mat, size_t r, size_t c, float num)
{
	//substract a scalar num to mat[l][r]
	return addPosition(mat, r, c, -num);
}
bool addScalar(Matrix *mat, float num)
{
    //add a scalar num to all the elements in matrix mat
    if (mat == NULL) return 0;
    for (int r = 0; r < mat->rows; r++) {
        for (int c = 0; c < mat->columns; c++) {
            (*pos(mat, r, c)) += num;
        }
    }
    return 1;
}
bool substrScalar(Matrix *mat, float num)
{
    //substract a scalar num to all the elements in matrix mat
    return addScalar(mat, -num);
}
bool multiScalar(Matrix *mat, float num)
{
    //make matrix mat (num) times
    if (mat == NULL) return 0;
    for (int r = 0; r < mat->rows; r++) {
        for (int c = 0; c < mat->columns; c++) {
            (*pos(mat, r, c)) *= num;
        }
    }
    return 1;
}
float minimalValue(const Matrix *mat)
{
    /*find the minimum number among the elements in the matrix
    DO NOT put a null matrix in this matrix, or it will return float 0*/
    if (mat == NULL) return 0;
    float ans = *(mat->datas);
    for (int r = 0; r < mat->rows; r++){
        for (int c = 0; c < mat->columns; c++) {
            if (ans > (*pos(mat, r, c))) ans = (*pos(mat, r, c));
        }
    }
    return ans;
}
float maximalValue(const Matrix *mat)
{
    /*find the maximum number among the elements in the matrix
    DO NOT put a null matrix in this matrix, or it will return float 0*/
    if (mat == NULL) return 0;
    float ans = *(mat->datas);
    for (int r = 0; r < mat->rows; r++){
        for (int c = 0; c < mat->columns; c++) {
            if (ans < (*pos(mat, r, c))) ans = (*pos(mat, r, c));
        }
    }
    return ans;
}
bool addMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make matrix ansMat equals to mat1 + mat2
    requires same size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (!(sameSize(mat1, mat2) && sameSize(mat2, ansMat))) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
	for (int i = 0; i < mat1->rows * mat1->columns; i++) {
		*(ansMat->datas + i) = *(mat1->datas + i) + *(mat2->datas + i);
	}
	return 1;
}
bool substractMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make matrix ansMat equals to mat1 + mat2
    requires same size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
	if (!(sameSize(mat1, mat2) && sameSize(mat2, ansMat))) return 0;
	for (int i = 0; i < mat1->rows * mat1->columns; i++) {
		*(ansMat->datas + i) = *(mat1->datas + i) - *(mat2->datas + i);
	}
	return 1;
}
bool copyMatrix(const Matrix *mat1, Matrix *mat2)
{
    /*copy the matrix mat1 to mat2
    requires same size*/
    if (mat1 == NULL || mat2 == NULL) return 0;
	if (!sameSize(mat1, mat2)) return 0;
	for (int i = 0; i < mat1->rows * mat1->columns; i++) {
		*(mat2->datas + i) = *(mat1->datas + i);
	}
	return 1;
}
bool multiMatrix(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make the matrix ansMat to mat1 * mat2
    /requires corresponding size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    for (int r = 0; r < ansMat->rows; r++) {
        for (int c = 0; c < ansMat->columns; c++){
            *pos(ansMat, r, c) = 0;
            for (int i = 0; i < mat1->columns; i++) {
                *pos(ansMat, r, c) += (*pos(mat1, r, i)) * (*pos(mat2, i, c));
            }
        }
    } 
    return 1;
}
inline float * pos2(const Matrix *mat, int r, int c)
{
    return mat->datas + (r * (mat->columns) + c);
}
inline float * rpos(float *dat, int r0, int c0, int c)
{
    return dat + (r0 * c + c0);
}
bool matmul_improved2(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make the matrix ansMat to mat1 * mat2
    /requires corresponding size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    int rr = ansMat->rows, cc = ansMat->columns, r2 = mat1->columns;
    for (register int i = 0; i < rr * cc; i++) {
        *(ansMat->datas + i) = 0;
    }
    #pragma omp parallel for schedule(dynamic)
    for (int r = 0; r < rr; r++) {
        for (int k = 0; k < r2; k++) {
            for (int c = 0; c < cc; c++){
                *(ansMat->datas + r * cc + c) += (*(mat1->datas + r * r2 + k)) * (*(mat2->datas + k * cc + c));
            }
        }
    } 
    return 1;
}
void randMatrix(Matrix *mat, int maxValue)
{
    for (int i = 0; i < mat->columns * mat->rows; i++) {
        mat->datas[i] = rand() % (maxValue + 1);
    }
}
void matmul_44(float *mat1, float *mat2, float *ansMat, int mat1_c, int mat2_c)
{//mat1: 
    register float d0, d1, d2, d3;
    register float k0;
    for (int i = 0; i < 4; ++i) {
        for (int k = 0; k < 4; ++k){
            k0 = *rpos(mat1, i, k, mat1_c);
            d0 = k0 * (*rpos(mat2, k, 0, mat2_c));
            d1 = k0 * (*rpos(mat2, k, 1, mat2_c));
            d2 = k0 * (*rpos(mat2, k, 2, mat2_c));
            d3 = k0 * (*rpos(mat2, k, 3, mat2_c));
            *rpos(ansMat, i, 0, mat2_c) += d0;
            *rpos(ansMat, i, 1, mat2_c) += d1;
            *rpos(ansMat, i, 2, mat2_c) += d2;
            *rpos(ansMat, i, 3, mat2_c) += d3;
        }
    }

}
bool matmul_improved3(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    int mat1_c = mat1->columns, mat1_r = mat1->rows, mat2_r = mat2->rows, mat2_c = mat2->columns;
    int mat1_cb = mat1_c >> 2, mat1_rb = mat1_r >> 2, mat2_rb = mat2_r >> 2, mat2_cb = mat2_c >> 2;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < mat1_rb; ++i) {
        for (int k = 0; k < mat1_cb; ++k) {
            for (int j = 0; j < mat2_cb; ++j) {
                matmul_44(rpos(mat1->datas, i << 2, k << 2, mat1_c), rpos(mat2->datas, k << 2, j << 2, mat2_c), 
                rpos(ansMat->datas, i << 2, j << 2, mat2_c), mat1_c, mat2_c);
            }
        }
    }
}
bool matmul_plain(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make the matrix ansMat to mat1 * mat2
    /requires corresponding size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
  //  #pragma omp parallel for num_threads(64)
    for (int r = 0; r < ansMat->rows; r++) {
        for (int c = 0; c < ansMat->columns; c++){
            *pos(ansMat, r, c) = 0;
            for (int i = 0; i < mat1->columns; i++) {
                *pos(ansMat, r, c) += (*pos(mat1, r, i)) * (*pos(mat2, i, c));
            }
        }
    } 
    return 1;
}
bool matmul_improved(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    float buffer[4];
    #pragma omp parallel for num_threads(64)
    for (int i = 0; i < mat1->rows; i++) {
        for (int j = 0; j < mat2->columns; j++) {
            int len = mat1->columns / 4;
            __m128 vsum = _mm_set1_ps(0.0f);
            for (int k = 0; k < len; ++k) {
                __m128 m1 = _mm_set_ps(*pos(mat1, i, k*4), *pos(mat1, i, k*4+1), *pos(mat1, i, k*4+2), *pos(mat1, i, k*4+3));
                __m128 m2 = _mm_set_ps(*pos(mat2, k*4, j), *pos(mat2, k*4+1, j), *pos(mat2, k*4+2, j), *pos(mat2, k*4+3, j));
                __m128 result = _mm_mul_ps(m1, m2);
                vsum = _mm_add_ps(vsum, result);
            }
            _mm_store_ps(buffer, vsum);
            for (int k = len * 4; k < mat1->columns; k++) {
                buffer[0] += (*pos(mat1, i, k)) * (*pos(mat2, k, j)); 
            }
            *pos(ansMat, i, j) = buffer[0] + buffer[1] + buffer[2] + buffer[3];
        }
    }
    return 1;
}
bool printMatrix(const Matrix *mat)
{
    /*print out the matrix*/
    if (mat == NULL) {
        printf("It is a NULL matrix!\n");
        return 0;
    }
    for (int r = 0; r < mat->rows; r++) {
        for (int c = 0; c < mat->columns; c++){
            printf("%f ", *pos(mat, r, c));
        }
        printf("\n");
    } 
    return 1;
}