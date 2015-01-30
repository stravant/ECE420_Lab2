
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <stdint.h>

#include "Lab2IO"
#include "sys/time.h"
#include "timer.h"
#include <math.h>

/* #define N 128, defined by command line invokation */

/* Matrices and sizes. */
int W[N*N];

/* Is async mode? */
int is_async = 0;

/* Thread info, an array of threads that we use for the calculation */
pthread_t *thread_array;
int thread_count;

/* Shared function to create and join threads */
void create_and_join_threads(void *(*thread_func)(void*)) {
	/* Start all of the threads */
	thread_array = calloc(thread_count, sizeof(pthread_t));

	/* Attributes */
	pthread_attr_t thread_attr;
	pthread_attr_init(&thread_attr);
	pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

	/* Start them */
	for (int rank = 0; rank < thread_count; ++rank) {
		int res = pthread_create(&thread_array[rank], 
			&thread_attr,
			thread_func,
			(void*)(intptr_t)rank);
		if (res != 0) {
			printf("Failed to create thread: %d\n", res);
			exit(0);
		}
	}

	/* Wait for all of the threads to finish */
	for (int rank = 0; rank < thread_count; ++rank) {
		pthread_join(thread_array[rank], NULL);
	}

	free(thread_array);
}

/* Code */
void floyd_warshall_single() {
	/* Single threaded approach */
	for (int k = 0; k < N; ++k) {
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				int sum = W[i*N + k] + W[k*N + j];
				if (sum < W[i*N + j])
					W[i*N + j] = sum;
			}
		}
	}
}

int async_row_current_k[N];
pthread_mutex_t async_lock[N];
pthread_cond_t async_row_completed[N];

void *fw_async_thread(void *rank_void) {
	/* Rank -> Which rows to calculate */
	int rank = (intptr_t)rank_void;
	int blk_sz = N / thread_count;
	int i_start = blk_sz*rank;

	for (int k = 0; k < N; ++k) {
		/* Do the work for a given K */
		for (int i = i_start; i < (i_start + blk_sz); ++i) {
			for (int j = 0; j < N; ++j) {
				/*
				 * At this point:
				 * W^(k-1)(i, k)
				 *   Must be complete, since this thread is the one
				 *   calculating i in this range
				 * W^(k-1)(k, j)
				 *   Must be complete if k is in our i range.
				 *   Otherwise, may or may not be complete, in that
				 *   case, we need to wait on that row's thread's 
				 *   condition until it is up to speed.
				 */
				if (k < i_start || k >= (i_start + blk_sz)) {
					/* Need to do wait */
					pthread_mutex_lock(&async_lock[k]);
					while (async_row_current_k[k] < (k-1)) {
						pthread_cond_wait(&async_row_completed[k],
							&async_lock[k]);
					}
					pthread_mutex_unlock(&async_lock[k]);
				}
				int sum = W[i*N + k] + W[k*N + j];
				if (sum < W[i*N + j])
					W[i*N + j] = sum;
			}

			/* Fire row completed */
			pthread_mutex_lock(&async_lock[i]);
			async_row_current_k[i] = k;
			pthread_cond_broadcast(&async_row_completed[i]);
			pthread_mutex_unlock(&async_lock[i]);
		}
	}

	/* Done */
	pthread_exit(0);
}

/* Multi threaded code */
void floyd_warshall_multi_async() {
	/* Multi threaded async approach */

	/* Start stuff */
	create_and_join_threads(fw_async_thread);
}

/* Synchronous solution shared values */
int sync_barrier_count;
pthread_mutex_t sync_barrier_lock;
pthread_cond_t sync_barrier_cond;

void *fw_sync_thread(void *rank_void) {
	/* Rank -> Which rows to calculate */
	int rank = (intptr_t)rank_void;
	int blk_sz = N / thread_count;
	int i_start = blk_sz*rank;

	for (int k = 0; k < N; ++k) {
		/* Do the work for a given K */
		for (int i = i_start; i < (i_start + blk_sz); ++i) {
			for (int j = 0; j < N; ++j) {
				int sum = W[i*N + k] + W[k*N + j];
				if (sum < W[i*N + j])
					W[i*N + j] = sum;
			}
		}

		/* Barrier Synchronize, pthread_barrier_t is not provided */
		pthread_mutex_lock(&sync_barrier_lock);
		if (--sync_barrier_count == 0) {
			sync_barrier_count = thread_count;
			pthread_cond_broadcast(&sync_barrier_cond);
			pthread_mutex_unlock(&sync_barrier_lock);
		} else {
			pthread_cond_wait(&sync_barrier_cond, &sync_barrier_lock);
			pthread_mutex_unlock(&sync_barrier_lock);
		}
	}

	/* Done */
	pthread_exit(0);
}

void floyd_warshall_multi_sync() {
	/* Multi threaded barriers approach */

	/* Init the barrier */
	sync_barrier_count = thread_count;
	pthread_mutex_init(&sync_barrier_lock, NULL);
	pthread_cond_init(&sync_barrier_cond, NULL);

	/* Start stuff */
	create_and_join_threads(fw_sync_thread);
}

/* Main function */
int main(int argc, char *argv[]) {
	/* What mode? */
	if (argc == 1) {
		thread_count = 1;
		is_async = 0;
	}
	if (argc > 1) {
		if (sscanf(argv[1], "%d", &thread_count) < 1) {
			printf("Bad thread count: `%s'\n", argv[1]);
			return EXIT_FAILURE;
		}
	}
	if (argc > 2) {
		if (sscanf(argv[2], "%d", &is_async) < 1) {
			printf("Bad async status: `%s'\n", argv[2]);
			return EXIT_FAILURE;
		}
	}

	/* Load in the input */
	if (Loaddata(W, N)) {
		return EXIT_FAILURE;
	}

	/* Start time */
	double start, end;
	GET_TIME(start);

	if (thread_count == 1) {
		floyd_warshall_single();
	} else {
		if (is_async) {
			floyd_warshall_multi_async();
		} else {
			floyd_warshall_multi_sync();
		}
	}

	/* End time */
	GET_TIME(end);

	/* Print perf */
	if (is_async) {
		printf("Asynchronous");
	} else {
		printf("Synchronous");
	}
	printf(" (threads=%d): %fs\n", thread_count, (end - start));
	
	/* Save output */
	if (Savedata(W, N)) {
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}