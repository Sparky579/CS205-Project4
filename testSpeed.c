#include"matrix.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
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
int testSpeed(int size, bool (*mul)(const Matrix*, const Matrix*, Matrix*))
{
    Matrix *x1 = createMatrix(size, size), *x2 = createMatrix(size, size), *x3 = createMatrix(size, size);
    for (int i = 0; i < size * size; i++) {
        x1->datas[i] = rand() % 100;
        x2->datas[i] = rand() % 100;
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
        return 1000*(end.tv_sec - start.tv_sec + end.tvusec - start.tvusec);
    #endif
}
bool assertEqual(int size, bool (*mul1)(const Matrix*, const Matrix*, Matrix*), bool (*mul2)(const Matrix*, const Matrix*, Matrix*))
{
    Matrix *x1 = createMatrix(size, size), *x2 = createMatrix(size, size), *x3 = createMatrix(size, size), *x4 = createMatrix(size, size);
    for (int i = 0; i < size * size; i++) {
        x1->datas[i] = rand() % 100;
        x2->datas[i] = rand() % 100;
    }
    mul1(x1, x2, x3);
    mul2(x1, x2, x4);
    bool flag = true;
    for (int i = 0; i < size * size; i ++) {
        if (x3->datas[i] != x4->datas[i]) flag = false;
    }
    return flag;
}
int main()
{
    if (!assertEqual(256, matmul_plain, matmul_improved3)) {
        printf("Error!");
        return 0;
     }
    Matrix *mat1, *mat2, *mat3;
    mat1 = createMatrix(8, 8);
    mat2 = createMatrix(8, 8);
    mat3 = createMatrix(8, 8);
    randMatrix(mat1, 3);
    randMatrix(mat2, 3);
    matmul_improved3(mat1, mat2, mat3);
    // matmul_44(pos2(mat1, 4, 4), pos2(mat2, 4, 4), pos2(mat3, 4, 4), 8, 8);
    printMatrix(mat1);printf("\n");
    printMatrix(mat2);printf("\n");
    printMatrix(mat3);printf("\n");
    if (!assertEqual(32, matmul_plain, matmul_improved2)) {
        printf("Error!");
        return 0;
    }
    testSpeed(256, matmul_plain);
    testSpeed(256, matmul_improved);
    int N = 1024;
    printf("%d %d %d %d\n", testSpeed(N / 4, matmul_plain), testSpeed(N, matmul_improved), testSpeed(N, matmul_improved2), testSpeed(512, matmul_improved3));
    system("pause");
    
}