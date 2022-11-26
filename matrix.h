#pragma once
#include <immintrin.h>
#include<stdlib.h>
#include<omp.h>
#include<malloc.h>
#include<string.h>
#include<stdio.h>
#include<stdbool.h>
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
inline float * rpos(float *dat, const int r0, const int c0, const int c)
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
    #pragma omp parallel for simd schedule(dynamic)
    for (int r = 0; r < rr; ++r) {
        for (int k = 0; k < r2; ++k) {
            for (int c = 0; c < cc; ++c){
                *(ansMat->datas + r * cc + c) += (*(mat1->datas + r * r2 + k)) * (*(mat2->datas + k * cc + c));
            }
        }
    } 
    return 1;
}
bool matmul_improved2_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    /*make the matrix ansMat to mat1 * mat2
    /requires corresponding size
    You shouldn't make ansMat equal to mat1 or mat2 */
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    int rr = ansMat->rows, cc = ansMat->columns, r2 = mat1->columns;
    for (int i = 0; i < rr * cc; i++) {
        *(ansMat->datas + i) = 0;
    }
    size_t cc_b = cc >> 2;
    __m128 result[cc_b + 1];
    float buf2[4];
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < rr * cc; i++)
    ansMat->datas[i] = 0.0f;
    for (int r = 0; r < rr; r++) {
        for (int k = 0; k < r2; k++) {
            for (int c = 0; c < cc_b; c++){
                const __m128 m0 = _mm_set1_ps(*(mat1->datas + r * r2 + k));
                __m128 m1 = _mm_set_ps(*(mat2->datas + k * cc + (c << 2)), *(mat2->datas + k * cc + (c << 2) + 1),
                *(mat2->datas + k * cc + (c << 2) + 2), *(mat2->datas + k * cc + (c << 2) + 3));
                m1 = _mm_mul_ps(m0, m1);
                if (k == 0)
                result[c] = m1;
                else result[c] = _mm_add_ps(result[c], m1);
            }
        }
        for (int c = 0; c < cc_b; c++) {
            _mm_store_ps(buf2, result[c]);
            for (int j = 0; j < 4; j++) {
                *(ansMat->datas + r * cc + (c << 2) + 3 - j) = buf2[j];
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
        d0 = d1 = d2 = d3 = 0;
        for (int k = 0; k < 4; ++k){
            k0 = *rpos(mat1, i, k, mat1_c);
            d0 += k0 * (*rpos(mat2, k, 0, mat2_c));
            d1 += k0 * (*rpos(mat2, k, 1, mat2_c));
            d2 += k0 * (*rpos(mat2, k, 2, mat2_c));
            d3 += k0 * (*rpos(mat2, k, 3, mat2_c));
        }
    //    printf("%d: %f %f %f %f\n",i, d0, d1, d2, d3);
        *rpos(ansMat, i, 0, mat2_c) += d0;
        *rpos(ansMat, i, 1, mat2_c) += d1;
        *rpos(ansMat, i, 2, mat2_c) += d2;
        *rpos(ansMat, i, 3, mat2_c) += d3;
    //    printf("%f %f %f %f\n", d0, d1, d2, d3);
    }

}
void matmul_44_1(float *mat1, float *mat2, float *ansMat, int mat1_c, int mat2_c)
{
    for (int i = 0; i < 4; ++i) {
        __m128 m0 = _mm_setzero_ps();
        float buffer[4];
        for (int k = 0; k < 4; ++k) {
            // k0 = _mm_set1_ps(*rpos(mat1, i, k, mat1_c));
            // k1 = _mm_set_ps(*rpos(mat2, k, 0, mat2_c), *rpos(mat2, k, 1, mat2_c), *rpos(mat2, k, 2, mat2_c), *rpos(mat2, k, 3, mat2_c));
        //    _mm_store_ps(buffer, k1);
            // printf("%f %f %f %f\n", buffer[3], buffer[2], buffer[1], buffer[0]);
            // printf("%f %f %f %f\n\n", *rpos(mat2, k, 0, mat2_c), *rpos(mat2, k, 1, mat2_c), *rpos(mat2, k, 2, mat2_c), *rpos(mat2, k, 3, mat2_c));
        //    k1 = _mm_set_ps(*rpos(mat2, k, 0, mat2_c), *rpos(mat2, k, 1, mat2_c), *rpos(mat2, k, 2, mat2_c), *rpos(mat2, k, 3, mat2_c));
            m0 = _mm_add_ps(m0, _mm_mul_ps(_mm_set1_ps(*rpos(mat1, i, k, mat1_c)), 
            _mm_set_ps(*rpos(mat2, k, 0, mat2_c), *rpos(mat2, k, 1, mat2_c), *rpos(mat2, k, 2, mat2_c), *rpos(mat2, k, 3, mat2_c))));
        }
        _mm_store_ps(buffer, m0);
        *rpos(ansMat, i, 0, mat2_c) += buffer[3];
        *rpos(ansMat, i, 1, mat2_c) += buffer[2];
        *rpos(ansMat, i, 2, mat2_c) += buffer[1];
        *rpos(ansMat, i, 3, mat2_c) += buffer[0];
    }
}
void matmul_44_11(float *mat1, float *mat2, float *ansMat, int mat1_c, int mat2_c)
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
    //    printf("%d: %f %f %f %f\n", i, *rpos(ansMat, i, 0, mat2_c), *rpos(ansMat, i, 1, mat2_c),
    //    *rpos(ansMat, i, 2, mat2_c), *rpos(ansMat, i, 3, mat2_c));
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
    for (int i = 0; i < ansMat->rows * ansMat->columns; ++i)
    ansMat->datas[i] = 0;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < mat1_rb; ++i) {
        for (int k = 0; k < mat1_cb; ++k) {
            for (int j = 0; j < mat2_cb; ++j) {
                matmul_44(rpos(mat1->datas, i << 2, k << 2, mat1_c), rpos(mat2->datas, k << 2, j << 2, mat2_c), 
                rpos(ansMat->datas, i << 2, j << 2, mat2_c), mat1_c, mat2_c);
            }
        }
    }
    for (int i = mat1_rb << 2; i < mat1_r; ++i)
        for (int k = 0; k < mat1_c; ++k) {
            for (int j = mat2_cb << 2; j < mat2_c; ++j) {
                *rpos(ansMat->datas, i, j, mat2_c) += (*rpos(mat1->datas, i, k, mat1_c)) * (*rpos(mat2->datas, k, j, mat2_c));
            }
        }
}
bool matmul_improved3_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
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
                matmul_44_1(rpos(mat1->datas, i << 2, k << 2, mat1_c), rpos(mat2->datas, k << 2, j << 2, mat2_c), 
                rpos(ansMat->datas, i << 2, j << 2, mat2_c), mat1_c, mat2_c);
            }
        }
    }
    for (int i = mat1_rb << 2; i < mat1_r; ++i)
        for (int k = 0; k < mat1_c; ++k) {
            for (int j = mat2_cb << 2; j < mat2_c; ++j) {
                *rpos(ansMat->datas, i, j, mat2_c) += (*rpos(mat1->datas, i, k, mat1_c)) * (*rpos(mat2->datas, k, j, mat2_c));
            }
        }
}
void matmul_improved4_mul44(float *mat1, float *mat2_t, float *ansMat, const size_t row1, const size_t col1, const size_t col2)
{
    
}

void matmul_improved5_mul44(float *mat1, float *mat2_t, float *ansMat, const size_t row1, const size_t col1, const size_t col2)
{
    register float f0, f1, f2, f3;
    for (int i = 0; i < 4; ++i) {
        f0 = f1 = f2 = f3 = 0.0;
        for (int k = 0; k < col1; ++k) {
            f0 = f0 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 0, k, col1));
            f1 = f1 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 1, k, col1));
            f2 = f2 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 2, k, col1));
            f3 = f3 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 3, k, col1));
        }
        *rpos(ansMat, i, 0, col2) = f0;
        *rpos(ansMat, i, 1, col2) = f1;
        *rpos(ansMat, i, 2, col2) = f2;
        *rpos(ansMat, i, 3, col2) = f3;
    }
}

void matmul_improved5_mulXX(float *mat1, float *mat2_t, float *ansMat, const size_t row1, const size_t col1, const size_t col2, const int BLOCK_SIZE)
{
    float f[BLOCK_SIZE];
    for (int i = 0; i < BLOCK_SIZE; ++i) {
        memset(f, 0, sizeof(f));
        for (int k = 0; k < col1; ++k) {
            for (int j = 0; j < BLOCK_SIZE; ++j) {
                f[j] += (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, j, k, col1));
            }
        }
        for (int j = 0; j < BLOCK_SIZE; ++j) 
            *rpos(ansMat, i, j, col2) = f[j];
    }
}

void matmul_improved5_mul88(float *mat1, float *mat2_t, float *ansMat, const size_t row1, const size_t col1, const size_t col2)
{
    register float f0, f1, f2, f3, f4, f5, f6, f7;
    for (int i = 0; i < 8; ++i) {
        f0 = f1 = f2 = f3 = f4 = f5 = f6 = f7 = 0.0;
        for (int k = 0; k < col1; ++k) {
            f0 = f0 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 0, k, col1));
            f1 = f1 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 1, k, col1));
            f2 = f2 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 2, k, col1));
            f3 = f3 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 3, k, col1));
            f4 = f4 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 4, k, col1));
            f5 = f5 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 5, k, col1));
            f6 = f6 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 6, k, col1));
            f7 = f7 + (*rpos(mat1, i, k, col1)) * (*rpos(mat2_t, 7, k, col1));
        }
        *rpos(ansMat, i, 0, col2) = f0;
        *rpos(ansMat, i, 1, col2) = f1;
        *rpos(ansMat, i, 2, col2) = f2;
        *rpos(ansMat, i, 3, col2) = f3;
        *rpos(ansMat, i, 4, col2) = f4;
        *rpos(ansMat, i, 5, col2) = f5;
        *rpos(ansMat, i, 6, col2) = f6;
        *rpos(ansMat, i, 7, col2) = f7;
    }
}

void matmul_improved5_mul44_1(float *mat1, float *mat2_t, float *ansMat, const size_t row1, const size_t col1, const size_t col2)
{
    register float f0, f1, f2, f3;
    f0 = f1 = f2 = f3 = 0.0;
    for (int k = 0; k < col1; ++k) {
        f0 = f0 + (*rpos(mat1, 0, k, col1)) * (*rpos(mat2_t, 0, k, col1));
        f1 = f1 + (*rpos(mat1, 0, k, col1)) * (*rpos(mat2_t, 1, k, col1));
        f2 = f2 + (*rpos(mat1, 0, k, col1)) * (*rpos(mat2_t, 2, k, col1));
        f3 = f3 + (*rpos(mat1, 0, k, col1)) * (*rpos(mat2_t, 3, k, col1));
    }
    *rpos(ansMat, 0, 0, col2) = f0;
    *rpos(ansMat, 0, 1, col2) = f1;
    *rpos(ansMat, 0, 2, col2) = f2;
    *rpos(ansMat, 0, 3, col2) = f3;
    f0 = f1 = f2 = f3 = 0.0;
    for (int k = 0; k < col1; ++k) {
        f0 = f0 + (*rpos(mat1, 1, k, col1)) * (*rpos(mat2_t, 0, k, col1));
        f1 = f1 + (*rpos(mat1, 1, k, col1)) * (*rpos(mat2_t, 1, k, col1));
        f2 = f2 + (*rpos(mat1, 1, k, col1)) * (*rpos(mat2_t, 2, k, col1));
        f3 = f3 + (*rpos(mat1, 1, k, col1)) * (*rpos(mat2_t, 3, k, col1));
    }
    *rpos(ansMat, 1, 0, col2) = f0;
    *rpos(ansMat, 1, 1, col2) = f1;
    *rpos(ansMat, 1, 2, col2) = f2;
    *rpos(ansMat, 1, 3, col2) = f3;
    f0 = f1 = f2 = f3 = 0.0;
    for (int k = 0; k < col1; ++k) {
        f0 = f0 + (*rpos(mat1, 2, k, col1)) * (*rpos(mat2_t, 0, k, col1));
        f1 = f1 + (*rpos(mat1, 2, k, col1)) * (*rpos(mat2_t, 1, k, col1));
        f2 = f2 + (*rpos(mat1, 2, k, col1)) * (*rpos(mat2_t, 2, k, col1));
        f3 = f3 + (*rpos(mat1, 2, k, col1)) * (*rpos(mat2_t, 3, k, col1));
    }
    *rpos(ansMat, 2, 0, col2) = f0;
    *rpos(ansMat, 2, 1, col2) = f1;
    *rpos(ansMat, 2, 2, col2) = f2;
    *rpos(ansMat, 2, 3, col2) = f3;
    f0 = f1 = f2 = f3 = 0.0;
    for (int k = 0; k < col1; ++k) {
        f0 = f0 + (*rpos(mat1, 3, k, col1)) * (*rpos(mat2_t, 0, k, col1));
        f1 = f1 + (*rpos(mat1, 3, k, col1)) * (*rpos(mat2_t, 1, k, col1));
        f2 = f2 + (*rpos(mat1, 3, k, col1)) * (*rpos(mat2_t, 2, k, col1));
        f3 = f3 + (*rpos(mat1, 3, k, col1)) * (*rpos(mat2_t, 3, k, col1));
    }
    *rpos(ansMat, 3, 0, col2) = f0;
    *rpos(ansMat, 3, 1, col2) = f1;
    *rpos(ansMat, 3, 2, col2) = f2;
    *rpos(ansMat, 3, 3, col2) = f3;
}
bool matmul_improved4(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    const Matrix *mat2_t = createMatrix(mat2->columns, mat2->rows);
    const size_t row1 = mat1->rows, column1 = mat1->columns, column2 = mat2->columns;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < column2; i++) {
        for (int j = 0; j < column1; j++) {
            *rpos(mat2_t->datas, j, i, column1) = *rpos(mat2->datas, i, j, column2);
        }
    }
    const size_t row1_b = row1 >> 2;
    const size_t column2_b = column2 >> 2;
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < row1_b; ++i) {
        for (int j = 0; j < column2_b; ++j) {
            matmul_improved5_mul44(rpos(mat1->datas, i << 2, 0, column1), rpos(mat2_t->datas, j << 2, 0, column1), rpos(ansMat->datas, i << 2, j << 2, column2), 
            row1, column1, column2);
        }
    }
    return 1;
}
bool matmul_improved4_1(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    const Matrix *mat2_t = createMatrix(mat2->columns, mat2->rows);
    const size_t row1 = mat1->rows, column1 = mat1->columns, column2 = mat2->columns;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < column2; i++) {
        for (int j = 0; j < column1; j++) {
            *rpos(mat2_t->datas, j, i, column1) = *rpos(mat2->datas, i, j, column2);
        }
    }
    const size_t row1_b = row1 >> 2;
    const size_t column2_b = column2 >> 2;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < row1_b; ++i) {
        for (int j = 0; j < column2_b; ++j) {
            matmul_improved5_mul44_1(rpos(mat1->datas, i << 2, 0, column1), rpos(mat2_t->datas, j << 2, 0, column1), rpos(ansMat->datas, i << 2, j << 2, column2), 
            row1, column1, column2);
        }
    }
    return 1;
}
bool matmul_improved4_2(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)//use 8*8 block
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    const Matrix *mat2_t = createMatrix(mat2->columns, mat2->rows);
    const size_t row1 = mat1->rows, column1 = mat1->columns, column2 = mat2->columns;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < column2; i++) {
        for (int j = 0; j < column1; j++) {
            *rpos(mat2_t->datas, j, i, column1) = *rpos(mat2->datas, i, j, column2);
        }
    }
    const size_t row1_b = row1 >> 3;
    const size_t column2_b = column2 >> 3;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < row1_b; ++i) {
        for (int j = 0; j < column2_b; ++j) {
            matmul_improved5_mul88(rpos(mat1->datas, i << 3, 0, column1), rpos(mat2_t->datas, j << 3, 0, column1), rpos(ansMat->datas, i << 3, j << 3, column2), 
            row1, column1, column2);
        }
    }
    return 1;
}
bool matmul_improved4_3(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    const Matrix *mat2_t = createMatrix(mat2->columns, mat2->rows);
    const size_t row1 = mat1->rows, column1 = mat1->columns, column2 = mat2->columns;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < column2; i++) {
        for (int j = 0; j < column1; j++) {
            *rpos(mat2_t->datas, j, i, column1) = *rpos(mat2->datas, i, j, column2);
        }
    }
    // for (int i = 0; i < row1 * column2; i++) {
    //     *(ansMat->datas + i) = 0;
    // }
    register float f0;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < row1; ++i) {
        for (int j = 0; j < column2; ++j) {
            f0 = 0;
            for (int k = 0; k < column1; ++k) {
                f0 += (*rpos(mat1->datas, i, k, column1)) * (*rpos(mat2_t->datas, j, k, column1));
            }
            *rpos(ansMat->datas, i, j, column2) = f0;
        }
    }
    // for (int i = 0; i < row1; ++i) {
    //     for (int k = 0; k < column1; ++k){
    //         for (int j = 0; j < column2; ++j)
    //             *rpos(ansMat->datas, i, j, column2) += (*rpos(mat1->datas, i, k, column1)) * (*rpos(mat2_t->datas, j, k, column1));
    //     }
    // }
    return 1;
}
bool matmul_improved4_X(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat, int BLOCK_SIZE)
{
    if (mat1 == NULL || mat2 == NULL || ansMat == NULL) return 0;
	if (mat1->columns != mat2->rows) return 0;
    if (mat1->rows != ansMat->rows || mat2->columns != ansMat->columns) return 0;
    if (mat1 == ansMat || mat2 == ansMat) return 0;
    const Matrix *mat2_t = createMatrix(mat2->columns, mat2->rows);
    const size_t row1 = mat1->rows, column1 = mat1->columns, column2 = mat2->columns;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < column2; i++) {
        for (int j = 0; j < column1; j++) {
            *rpos(mat2_t->datas, j, i, column1) = *rpos(mat2->datas, i, j, column2);
        }
    }
    const size_t row1_b = row1 / BLOCK_SIZE;
    const size_t column2_b = column2 / BLOCK_SIZE;
    #pragma omp parallel for simd schedule(dynamic)
    for (int i = 0; i < row1_b; ++i) {
        for (int j = 0; j < column2_b; ++j) {
            matmul_improved5_mulXX(rpos(mat1->datas, i * BLOCK_SIZE, 0, column1), rpos(mat2_t->datas, j * BLOCK_SIZE, 0, column1), rpos(ansMat->datas, i * BLOCK_SIZE, j * BLOCK_SIZE, column2), 
            row1, column1, column2, BLOCK_SIZE);
        }
    }
    return 1;
}
bool matmul_improved4_X_1_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 1);    
}
bool matmul_improved4_X_2_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 2);    
}
bool matmul_improved4_X_4_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 4);    
}
bool matmul_improved4_X_8_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 8);    
}
bool matmul_improved4_X_16_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 16);    
}
bool matmul_improved4_X_32_test(const Matrix *mat1, const Matrix *mat2, Matrix *ansMat) {
    matmul_improved4_X(mat1, mat2, ansMat, 32);    
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
