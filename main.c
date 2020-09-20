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
#include <float.h>
#include <pmmintrin.h>

#define EPS 0.1

struct Vector_info {
    float a;
    float b;
    float x0;
    size_t size_N;
};

struct Result_by_blocks{
    float **arr;
    int blocks_count;
    int block_size;
};

float get_num_bloks_vec (struct Result_by_blocks vec_by_block, size_t iter){
  return vec_by_block.arr[iter/vec_by_block.block_size][iter%vec_by_block.block_size];
}


// colculation of vector's values consistently
float* calcul_n_parallel(struct Vector_info inp, long long *t0, long long *t1){
    float* result = calloc(inp.size_N, sizeof(float));
     float a = inp.a;
     float b = inp.b;
     float x0 = inp.x0;
     size_t size_N = inp.size_N;

    result[0] = x0;


    *t0 = clock();
   for(int i = 1; i < size_N; ++i)
       result[i] = a*result[i-1] + b;
    *t1 = clock();
    return result;

}


float* calcul_n_parallel_v2(struct Vector_info inp, long long *t0, long long *t1){
    float* result = calloc(inp.size_N, sizeof(float));
     float a = inp.a;
     float b = inp.b;
     float x0 = inp.x0;
     size_t size_N = inp.size_N;

    result[0] = x0;


    *t0 = clock();
    float tmp_xi = x0;
   for(int i = 1; i < size_N; ++i){
       tmp_xi = a*tmp_xi + b;
       result[i] = tmp_xi;
   }
    *t1 = clock();
    return result;

}

// 
float* calcul_n_parallel_with_AVX(struct Vector_info inp, long long *t0, long long *t1){
    float* result = calloc(inp.size_N, sizeof(float));
     float a = inp.a;
     float b = inp.b;
     float x0 = inp.x0;
     size_t size_N = inp.size_N;

    result[0] = x0;
    result[1] = a*result[0] + b;
    result[2] = a*result[1] + b;
    result[3] = a*result[2] + b;
    
    float tmp_vec_a[4];
    float tmp_vec_b[4];
    
    for(int i = 0; i < 4; ++i){
        tmp_vec_a[i] = a;
        tmp_vec_b[i] = b;
    }
    
    __m128 a_vec = _mm_load_ps(&tmp_vec_a[0]);
    __m128 b_vec = _mm_load_ps(&tmp_vec_b[0]);
    
    __m128 sum = _mm_load_ps(&tmp_vec_b[0]);
    __m128 pow_a = _mm_mul_ps(a_vec, a_vec);
    __m128 tmp1 = _mm_mul_ps(a_vec, b_vec);
    __m128 tmp2 = sum;
    sum = _mm_add_ps(sum, tmp1);
    tmp1 = _mm_mul_ps(pow_a, a_vec); // a^3
    pow_a = tmp1;
    tmp1 = _mm_mul_ps(pow_a, b_vec);
    tmp2 = sum;
    sum = _mm_add_ps(tmp2, tmp1);
    tmp2 = _mm_mul_ps(pow_a, a_vec);// a^4
    pow_a = tmp2;
    
    *t0 = clock();
   for(int i = 4; i < size_N; i+=4){
       __m128 x_vec = _mm_load_ps(&result[i-4]);
       __m128 tmp_vec = _mm_mul_ps(pow_a, x_vec);
       x_vec = _mm_add_ps(tmp_vec, sum);
       _mm_store_ps(&result[i],x_vec);
   }
    *t1 = clock();
    return result;

}

// calculate elements values parralel
// creating result array for each thread
 struct Result_by_blocks calcul(struct Vector_info inp, long long *t0, long long *t1){

   float a = inp.a;
   float b = inp.b;
   float x0 = inp.x0;
   size_t size_N = inp.size_N;

   size_t arr_size = size_N/4;

    float **result = calloc(4, sizeof(float*));
    for (int i = 0; i < 4; ++i)
        result[i] = calloc(arr_size, sizeof(float));



     float arr_x0[4];

     float tmp_pow = pow(a, (float)arr_size);
     float seq_sum = (float)(1 - tmp_pow)/(1 - a);


    printf("\n-----debag-----\na = %f", a);
    printf("\nb = %f\narr_size = %ld\n", b, arr_size);
    printf("pow = %f\na = %f\n", tmp_pow, a);
     printf("x0 = %f\n", inp.x0);
    printf("rez = %f", seq_sum);

    arr_x0[0] = (float)x0;
    arr_x0[1] = tmp_pow*arr_x0[0] + b*seq_sum;
    arr_x0[2] = tmp_pow*arr_x0[1] + b*seq_sum;
    arr_x0[3] = tmp_pow*arr_x0[2] + b*seq_sum;


    *t0 = clock();
#pragma omp parralel for
    for (int i = 0; i < 4; ++i){
        float *rez0 = result[i];
        rez0[0] = arr_x0[i];
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
 float* calcul_v2(struct Vector_info inp, long long *t0, long long *t1){
    float* result = calloc(inp.size_N, sizeof(float));
    float a = inp.a;
    float b = inp.b;
    float x0 = inp.x0;
    size_t size_N = inp.size_N;
    size_t arr_size = size_N/4;

     float arr_x0[4];

     float tmp_pow = pow(a, (float)arr_size);
     float seq_sum = (float)(1 - tmp_pow)/(1 - a);


    printf("\n-----debag-----\na = %f", a);
    printf("\nb = %f\narr_size = %ld\n", b, arr_size);
    printf("pow = %f\na = %f\n", tmp_pow, a);
     printf("x0 = %f\n", inp.x0);
    printf("rez = %f", seq_sum);

    arr_x0[0] = (float)x0;
    arr_x0[1] = tmp_pow*arr_x0[0] + b*seq_sum;
    arr_x0[2] = tmp_pow*arr_x0[1] + b*seq_sum;
    arr_x0[3] = tmp_pow*arr_x0[2] + b*seq_sum;



    float rez_tmp[size_N];
    *t0 = clock();
#pragma omp parralel for
    for (int i = 0; i < 4; ++i){
        float *rez0 = result + arr_size*i;
        rez0[0] = arr_x0[i];
        for(size_t j = 1; j < arr_size; ++j)
             rez0[j] = a*rez0[j-1] + b;
    }
    *t1 = clock();

    return result;
}



int main(int argc, const char * argv[]) {
    long long t0, t1, t0_np, t1_np, tv20, tv21, t0_np_v2, t1_np_v2, t0_np_AVX, t1_np_AVX;
      float a, b, x0;
    size_t size_N;
    FILE *finp;
    finp = fopen("inp.txt", "r");
    fscanf(finp, "%f", &a);
    fscanf(finp, "%f", &b);
    fscanf(finp, "%f", &x0);
    fscanf(finp, "%ld", &size_N);
    fclose(finp);
    printf("a = %f, b = %f, x0 = %f, size = %ld\n", a, b, x0, size_N);

    struct Vector_info calc_inp;
    calc_inp.a = a;
    calc_inp.b = b;
    calc_inp.x0 = x0;
    calc_inp.size_N = size_N;

    struct Result_by_blocks result = calcul(calc_inp, &t0, &t1);
    float* result_n_p= calcul_n_parallel(calc_inp, &t0_np, &t1_np);
    float* result_v2= calcul_v2(calc_inp, &tv20, &tv21);
    float* result_n_p_v2= calcul_n_parallel_v2(calc_inp, &t0_np_v2, &t1_np_v2);
    float* result_n_p_AVX= calcul_n_parallel_with_AVX(calc_inp, &t0_np_AVX, &t1_np_AVX);
    FILE *output;
    output = fopen("result.txt", "w");
    for(size_t i = 0; i < size_N; ++i){
        fprintf(output, "%f", result_n_p_AVX[i]);
        if (abs(result_n_p[i] - result_n_p_AVX[i]) < EPS)
            fprintf(output, " Ok ");
        else
            fprintf(output, " Err ");
        fprintf(output, "%f \n", result_n_p[i]);
    }
    fclose(output);
    //free(result);
    free(result_n_p);
    printf("\n==== Time of the not parralel method ====\n");
    printf("\ntime = %lld\n", t1_np-t0_np);
    printf("\n==== Time of the not parralel method v2 ====\n");
    printf("\ntime = %lld\n", t1_np_v2-t0_np_v2);
        printf("\n==== Time of the not parralel method with AVX ====\n");
    printf("\ntime = %lld\n", t1_np_AVX-t0_np_AVX);
    printf("\n==== Time of the parralel method ====\n");
    printf("\ntime = %lld\n", t1-t0);
    printf("\n==== Time of the parralel method v2 ====\n");
    printf("\ntime = %lld\n", tv21-tv20);
    return 0;
}
