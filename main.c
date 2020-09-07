//
//  main.c
//  Task_2
//
//  Created by ArtTikidji on 26/02/2020.
//  Copyright © 2020 ArtTikidji. All rights reserved.
//

#define __USE_MINGW_ANSI_STDIO

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define EPS 0.1

struct Vector_info {
    long double a;
    long double b;
    long double x0;
    size_t size_N;
};

struct Result_by_blocks{
    long double **arr;
    int blocks_count;
    int block_size;
};

long double get_num_bloks_vec (struct Result_by_blocks vec_by_block, size_t iter){
  return vec_by_block.arr[iter/vec_by_block.block_size][iter%vec_by_block.block_size];
}

// colculation of vector's values consistently
long double* calcul_n_parallel(struct Vector_info inp, long long *t0, long long *t1){
        long double* result = calloc(inp.size_N, sizeof(long double));
     long double a = inp.a;
     long double b = inp.b;
     long double x0 = inp.x0;
     size_t size_N = inp.size_N;

    result[0] = x0;


    *t0 = clock();
   for(int i = 1; i < size_N; ++i)
       result[i] = a*result[i-1] + b;
    *t1 = clock();
    return result;

}


// calculate elements values parralel
// creating result array for each thread
 struct Result_by_blocks calcul(struct Vector_info inp, long long *t0, long long *t1){

   long double a = inp.a;
   long double b = inp.b;
   long double x0 = inp.x0;
   size_t size_N = inp.size_N;

   size_t arr_size = size_N/4;

    long double **result = calloc(4, sizeof(long double*));
    for (int i = 0; i < 4; ++i)
        result[i] = calloc(arr_size, sizeof(long double));


    result[0][0] = x0;

     long double a0, a1, a2, a3, arr[4];

     long double tmp_pow = pow(a, (long double)arr_size);
     long double seq_sum = (long double)(1 - tmp_pow)/(1 - a);
    a0 = (long double)x0;
    a1 = tmp_pow*a0 + b*seq_sum;
    a2 = tmp_pow*a1 + b*seq_sum;
    a3 = tmp_pow*a2 + b*seq_sum;


    printf("\n-----debag-----\na = %Lf", a);
    printf("\nb = %Lf\narr_size = %ld\n", b, arr_size);
    printf("pow = %Lf\na = %Lf\n", tmp_pow, a);
     printf("x0 = %Lf\n", inp.x0);
    printf("rez = %Lf", seq_sum);

    arr[0] = a0;
    arr[1] = a1;
    arr[2] = a2;
    arr[3] = a3;


    *t0 = clock();
#pragma omp parralel for
    for (int i = 0; i < 4; ++i){
        long double *rez0 = result[i];
        rez0[0] = arr[i];
        for(size_t j = 1; j < arr_size; ++j)
             rez0[j] = a*rez0[j-1] + b;
    }
    *t1 = clock();
    struct Result_by_blocks rez;
    rez.arr = result;
    rez.block_size = arr_size;
    return rez;
}


// calculate elements values parralel
// it dosen't create result array for each thread, and trying to work with critical section by changing pointer
 long double* calcul_v2(struct Vector_info inp, long long *t0, long long *t1){
    long double* result = calloc(inp.size_N, sizeof(long double));
    long double a = inp.a;
    long double b = inp.b;
    long double x0 = inp.x0;
    size_t size_N = inp.size_N;

    result[0] = x0;

    long double a0, a1, a2, a3, arr[4];
    size_t arr_size = size_N/4;
    long double tmp_pow = pow(a, (long double)arr_size);
    long double seq_sum = (long double)(1 - tmp_pow)/(1 - a);
    printf("\n-----debag-----\na = %Lf", a);
    printf("\nb = %Lf\narr_size = %ld\n", b, arr_size);
    printf("pow = %Lf\na = %Lf\n", tmp_pow, a);
     printf("x0 = %Lf\n", inp.x0);
    printf("rez = %Lf", seq_sum);
    a0 = (long double)x0;
    a1 = tmp_pow*a0 + b*seq_sum;
    a2 = tmp_pow*a1 + b*seq_sum;
    a3 = tmp_pow*a2 + b*seq_sum;

    arr[0] = a0;
    arr[1] = a1;
    arr[2] = a2;
    arr[3] = a3;

    printf("\n a0, a1, a2, a3\n %Lf", a0);
    printf(" %Lf", a1);
    printf(" %Lf", a2);
    printf(" %Lf\n", a3);

    long double rez_tmp[size_N];
    *t0 = clock();
#pragma omp parralel for
    for (int i = 0; i < 4; ++i){
        //long double rez0[arr_size];
        long double *rez0 = result + arr_size*i;
        rez0[0] = arr[i];
        for(size_t j = 1; j < arr_size; ++j)
             rez0[j] = a*rez0[j-1] + b;
        //for(size_t j = 0; j < arr_size; ++j)
        //    rez_tmp[arr_size*i + j] = rez0[j];
    }
    *t1 = clock();

    //for(size_t i = 0; i < size_N; ++i)
    //    result[i] = rez_tmp[i];
    return result;
}



int main(int argc, const char * argv[]) {
    long long t0, t1, t0_np, t1_np, tv20, tv21;
      long double a, b, x0;
    size_t size_N;
    FILE *finp;
    finp = fopen("inp.txt", "r");
    fscanf(finp, "%Lf", &a);
    fscanf(finp, "%Lf", &b);
    fscanf(finp, "%Lf", &x0);
    fscanf(finp, "%ld", &size_N);
    fclose(finp);
    printf("a = %Lf, b = %Lf, x0 = %Lf, size = %ld\n", a, b, x0, size_N);

    struct Vector_info calc_inp;
    calc_inp.a = a;
    calc_inp.b = b;
    calc_inp.x0 = x0;
    calc_inp.size_N = size_N;

    struct Result_by_blocks result = calcul(calc_inp, &t0, &t1);
    long double* result_n_p= calcul_n_parallel(calc_inp, &t0_np, &t1_np);
    long double* result_v2= calcul_v2(calc_inp, &tv20, &tv21);

    FILE *output;
    output = fopen("result.txt", "w");
    for(size_t i = 0; i < size_N; ++i){
        fprintf(output, "%Lf", get_num_bloks_vec(result, i));
        if (abs(result_n_p[i] - get_num_bloks_vec(result, i)) < EPS)
            fprintf(output, " Ok ");
        else
            fprintf(output, " Err ");
        fprintf(output, "%Lf \n", result_n_p[i]);
    }
    fclose(output);
    //free(result);
    free(result_n_p);
    printf("\n==== Time of the not parralel method ====\n");
    printf("\ntime = %lld\n", t1_np-t0_np);
    printf("\n==== Time of the parralel method ====\n");
    printf("\ntime = %lld\n", t1-t0);
    printf("\n==== Time of the parralel method v2 ====\n");
    printf("\ntime = %lld\n", tv21-tv20);
    return 0;
}
