#include"matrix.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
#include<cblas.h>
#define debug0
#if defined __WIN32 || defined __WIN64
#include<windows.h>
#endif
#ifdef __linux__
#include<time.h>
#include<sys/time.h>
#endif
int s1[10] = {1,2,3,4,3,2,1,3,5};
int s2[10] = {3,1,2,3,1,2,3,1,3};
const int N = 512;
long long testSpeed(int size, bool (*mul)(const Matrix*, const Matrix*, Matrix*))
{
    Matrix *x1 = createMatrix(size, size), *x2 = createMatrix(size, size), *x3 = createMatrix(size, size);
    for (int i = 0; i < size * size; i++) {
        x1->datas[i] = rand();
        x2->datas[i] = rand();
    }
    #if defined __Win32 || defined __WIN64
        int begintime = clock();
        mul(x1, x2, x3);
        int endtime = clock();
        return endtime - begintime;
    #endif
    #ifdef __linux__
        struct timeval start, end;
        gettimeofday(&start, NULL);
        mul(x1, x2, x3);
        gettimeofday(&end, NULL);
        return 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    #endif
}
bool assertEqual(int size, bool (*mul1)(const Matrix*, const Matrix*, Matrix*), bool (*mul2)(const Matrix*, const Matrix*, Matrix*))
{
    Matrix *x1 = createMatrix(size, size), *x2 = createMatrix(size, size), *x3 = createMatrix(size, size), *x4 = createMatrix(size, size);
    for (int i = 0; i < size * size; i++) {
        x1->datas[i] = rand() % 10;
        x2->datas[i] = rand() % 10;
    }
    mul1(x1, x2, x3);
//    printf("\n\n\n");
    mul2(x1, x2, x4);
    bool flag = true;
    // printMatrix(x1);printf("\n");
    // printMatrix(x2);printf("\n");
    // printMatrix(x3);printf("\n");
    // printMatrix(x4);printf("\n");
    for (int i = 0; i < size * size; i ++) {
        if (x3->datas[i] != x4->datas[i]) flag = false;
    }
    return flag;
}
int main()
{

    int N = 1024;
    for (int i = N; i <= N; i <<= 1) {
        printf("%d: %lld\n", i, testSpeed(i, matmul_improved4_1));
    }
    
    float A[N * N], B[N * N], C[N * N];
    for (int i = 0; i < N * N; i++){
        for (int j = 0; j < N * N; j++) {
            A[i * N + j] = rand() %100;
            B[i * N + j] = rand() % 100;
        }
    }
    struct timeval start, end;
    gettimeofday(&start, NULL);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1, A, N, B, N, 0, C, N);
    gettimeofday(&end, NULL);
    printf("%lld\n", 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec);

    
}
