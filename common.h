
#ifndef __GA_COMMON_H__
#define __GA_COMMON_H__
#endif
long long tvdiff_us(struct timeval *tv1, struct timeval *tv2);
long long ga_usleep(long long interval, struct timeval *ptv);
int*** init3dArrayInt(int, int , int);
int** init2dArrayInt(int, int);
double*** init3dArrayDouble(int, int , int);
double** init2dArrayDouble(int, int);
double avg(double* , int);
double sum(double* , int);
void get_face_tid(int No_face, int No_tile_h, int No_tile_v, int tid, int* face_id, int* tile_id);
void showArrayInt(int *arr, int N);
void showArrayDouble(double *arr, int N);




