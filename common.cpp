#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#ifndef WIN32
#ifndef ANDROID
#include <execinfo.h>
#endif /* !ANDROID */
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/syscall.h>
#endif /* !WIN32 */
#ifdef ANDROID
#include <android/log.h>
#endif /* ANDROID */
#ifdef __APPLE__
#include <syslog.h>
#endif

#if !defined(WIN32) && !defined(__APPLE__) && !defined(ANDROID)
#include <X11/Xlib.h>
#endif

#include "common.h"

#include <map>
#include <list>
using namespace std;

/**
 * Compute the time difference for two \a timeval data structure, i.e.,
 * \a tv1 - \a tv2.
 *
 * @param tv1 [in] Pointer to the first \a timeval data structure.
 * @param tv2 [in] Pointer to the second \a timeval data structure.
 * @return The difference in micro seconds.
 */
long long
tvdiff_us(struct timeval *tv1, struct timeval *tv2) {
	struct timeval delta;
	delta.tv_sec = tv1->tv_sec - tv2->tv_sec;
	delta.tv_usec = tv1->tv_usec - tv2->tv_usec;
	if(delta.tv_usec < 0) {
		delta.tv_sec--;
		delta.tv_usec += 1000000;
	}
	return 1000000LL*delta.tv_sec + delta.tv_usec;
}

/**
 * Sleep and wake up at \a ptv + \a interval (micro seconds).
 *
 * @param interval [in] The expected sleeping time (in micro seconds).
 * @param ptv [in] Pointer to the baseline time.
 * @return Currently always return 0.
 *
 * This function is useful for controlling precise sleeps.
 * We usually have to process each video frame in a fixed interval.
 * Each time interval includes the processing time and the sleeping time.
 * However, the processing time could be different in each iteration, so
 * the sleeping time has to be adapted as well.
 * To achieve the goal, we have to obtain the baseline time \a ptv
 * (using \a gettimeofday function)
 * \em before the processing task and call this function \em after
 * the processing task. In this case, the \a interval is set to the total
 * length of the interval, e.g., 41667 for 24fps video.
 *
 * This function sleeps for \a interval micro seconds if the baseline
 * time is not specified.
 */
long long
ga_usleep(long long interval, struct timeval *ptv) {
	long long delta;
	struct timeval tv;
	if(ptv != NULL) {
		gettimeofday(&tv, NULL);
		delta = tvdiff_us(&tv, ptv);
		if(delta >= interval) {
			usleep(1);
			return -1;
		}
		interval -= delta;
	}
	usleep(interval);
	return 0LL;
}
int*** init3dArrayInt(int M, int N, int P){
	int*** ret;
	int i,j;
	ret = new int**[M];
	for(i=0; i < M; i++){
		ret[i] = new int*[N];
		for(j=0; j < N; j++)
			ret[i][j] = new int[P];
	}
	return ret;
}
int** init2dArrayInt(int M, int N){
	int** ret;
	int i,j;
	ret = new int*[M];
	for(i=0; i < M; i++){
		ret[i] = new int[N];
	}
	return ret;
}
double*** init3dArrayDouble(int M, int N, int P){
	double*** ret;
	int i,j;
	ret = new double**[M];
	for(i=0; i < M; i++){
		ret[i] = new double*[N];
		for(j=0; j < N; j++)
			ret[i][j] = new double[P];
	}
	return ret;
}
double** init2dArrayDouble(int M, int N){
	double** ret;
	int i,j;
	ret = new double*[M];
	for(i=0; i < M; i++){
		ret[i] = new double[N];
	}
	return ret;
}
double avg(double* a, int N){
	double avg=0;
	for(int i=0; i < N; i++)
		avg += a[i] / N;
	return avg;
}
double sum(double* a, int N){
	double sum=0;
	for(int i=0; i < N; i++)
		sum += a[i];
	return sum;
}
void get_face_tid(int No_face, int No_tile_h, int No_tile_v, int tid, int* face_id, int* tile_id){
	if(No_face == 1 || No_face == 6){
		*face_id = tid / (No_tile_h * No_tile_v);
		*tile_id = tid - (*face_id) * No_tile_h * No_tile_v;
	}
	if(No_face == 2){
		*face_id = 0;
		*tile_id = tid;
	}
}
void showArrayInt(int *arr, int N){
	int i;
	for(i=0; i < N; i++){
		printf("%d, ", arr[i]);
	}
	printf("\n");
}
void showArrayDouble(double *arr, int N){
	int i;
	for(i=0; i < N; i++){
		printf("%.2f, ", arr[i]);
	}
	printf("\n");
}