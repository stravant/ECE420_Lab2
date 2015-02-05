
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <stdint.h>
#include <semaphore.h>
#include <memory.h>

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

/* Execution start */
double start_time;

double elapsed() {
	double end;
	GET_TIME(end);
	return 1000*(end - start_time);
}

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

/* Interesting implementation */
int (*WS)[N*N];
sem_t *WS_sem;

void dump(int k, const char *format, ...) {
	char name_buffer[20];
	snprintf(name_buffer, 20, "thlog_%d", k);
	FILE *f = fopen(name_buffer, "a");

	va_list args;
	va_start(args, format);
	vfprintf(f, format, args);
	va_end(args);

	fclose(f);
}

void dump_buf(int k, int buf) {
	char buffer[(N*(N+3))*10];
	char *output = buffer;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			output += sprintf(output, "%5d ", WS[buf][i*N + j]);
		}
		output += sprintf(output, "\n");
	}
	output += sprintf(output, "\n");
	dump(k, "%s", buffer);
}

void cleardump() {
	for (int i = 0; i < N; ++i) {
		char name_buffer[20];
		snprintf(name_buffer, 20, "thlog_%d", i);
		remove(name_buffer);
	}
}

void *fw_striped_thread(void *rank_void) {
	int rank = (intptr_t)rank_void;

	/* For each K that we are processing */
	for (int k = rank; k < N; k += thread_count) {
		/* Figure out where our data source and target are in the WS ring
		 * buffer. */
		int src = (k)   % (thread_count + 1);
		int dst = (k+1) % (thread_count + 1);
		int *src_W = WS[src];
		int *dst_W = WS[dst];

		//dump(k, "Start doing k=%d at %f [src %d -> dst %d]\n", 
		//	k, elapsed(), src, dst);

		/* Need to wait for k+1 rows to be completed in the source
		 * before we can proceed. */
		for (int dep_rows = 0; dep_rows < (k+1); ++dep_rows) {
			sem_wait(&WS_sem[src]);
		}

		/* Now loop over the slices */
		for (int slice = 0; slice < N; ++slice) {
			/* We also need to wait for the additional rows after the
			 * main row or col = k dependency */
			if (slice > k && k > 0) {			
				sem_wait(&WS_sem[src]);
			}

			//dump(k, "Doing slice s=%d at %f\n", slice, elapsed());
			//dump_buf(k, dst);

			/* Now, for each element in the slice */
			for (int i = slice+1; i < N; ++i) {
				/* Coord = (i, slice) */
				int val = src_W[i*N + k] + src_W[k*N + slice];
				if (val < src_W[i*N + slice]) {
					dst_W[i*N + slice] = val;
				} else {
					dst_W[i*N + slice] = src_W[i*N + slice];
				}
			}
			for (int j = slice; j < N; ++j) {
				/* Coord = (slice, j) */
				int val = src_W[slice*N + k] + src_W[k*N + j];
				if (val < src_W[slice*N + j]) {
					dst_W[slice*N + j] = val;
				} else {
					dst_W[slice*N + j] = src_W[slice*N + j];
				}
			}

			/* Notify the next K that a slice was completed  */
			sem_post(&WS_sem[dst]);
		}

		//dump(k, "Finished k=%d at %f:\n", k, elapsed());
		//dump_buf(k, dst);
	}

	pthread_exit(0);
}

void floyd_warshall_striped() {
	//cleardump();

	/* Set up the memory */
	WS = malloc(N*N*(thread_count+1)*sizeof(int));
	WS_sem = malloc((thread_count+1)*sizeof(sem_t));

	/* Move the data into WS initial */
	memcpy(WS[0], W, N*N*sizeof(int));

	/* Init the semaphores */
	sem_init(&WS_sem[0], 0, 1); /* Initial value = 1 */
	for (int i = 1; i < (thread_count+1); ++i) {
		sem_init(&WS_sem[i], 0, 0); /* Initial value = 0 */
	}

	/* Do the work */
	create_and_join_threads(fw_striped_thread);

	/* Move the result into W */
	int finalDst = N % (thread_count+1);
	memcpy(W, WS[finalDst], N*N*sizeof(int));
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
	GET_TIME(start_time);

	// if (thread_count == 1) {
	// 	floyd_warshall_single();
	// } else {
		if (is_async) {
			floyd_warshall_striped();
		} else {
			floyd_warshall_multi_sync();
		}
	// }

	/* End time */
	double end;
	GET_TIME(end);

	/* Print perf */
	if (is_async) {
		printf("Asynchronous");
	} else {
		printf("Synchronous");
	}
	printf(" (threads=%d): %fs\n", thread_count, (end - start_time));
	
	/* Save output */
	if (Savedata(W, N)) {
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}