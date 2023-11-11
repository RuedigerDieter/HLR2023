/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include "partdiff.h"
#include <pthread.h>
#include <semaphore.h>

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};


struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_old (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;           /* local variables for loops */
	int m1, m2;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

struct thread_arg{
	int start_index;
	int work_length;
	double* maxResiduum;
	sem_t* maxResiduum_sem;
	struct calculation_arguments* arguments;
	struct options* options;
	int *m1;
	int *m2; 
};

void *runThread(void *args)
{
	struct thread_arg *thread_args = (struct thread_arg*) args;
	int start_index = thread_args->start_index;
	int work_length = thread_args->work_length;
	struct calculation_arguments* arguments = thread_args->arguments;
	struct options* options = thread_args->options;
	
	int* m1 = thread_args->m1;
	int* m2 = thread_args->m2;
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	int row_len = (N-1);
	double** Matrix_Out = arguments->Matrix[*m1];
	double** Matrix_In  = arguments->Matrix[*m2];

	double fpisin_i = 0.0;

	if (options->inf_func == FUNC_FPISIN)
	{
		fpisin_i = fpisin * sin(0);
	}

	for(int x = start_index; x < start_index + work_length; x++) {
		int i = (x / row_len) + 1;
		int j = (x % row_len) + 1;

		if(j == N) { //if next iteration will begin at j=0

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

		}
		
		star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

		if (options->inf_func == FUNC_FPISIN)
		{
			star += fpisin_i * sin(pih * (double)j);
		}

		if (options->termination == TERM_PREC || options->term_iteration == 1)
		{
			residuum = Matrix_In[i][j] - star;
			residuum = (residuum < 0) ? -residuum : residuum;

			sem_wait(thread_args->maxResiduum_sem);
			double mr = *(thread_args->maxResiduum);
			*(thread_args->maxResiduum) = (residuum < mr) ? mr : residuum;
			sem_post(thread_args->maxResiduum_sem);
		}

		Matrix_Out[i][j] = star;

	}

	return NULL;

}

int createThreads(struct calculation_arguments* arguments, struct options* options, pthread_t** threads, struct thread_arg** thread_args)
{
    int t = options->number;
    int N = arguments->N;
    
    int M = (N-1) * (N-1);
    int L = (int) (M / t);
    int R = M - L * t;

	*threads = (pthread_t*) allocateMemory(sizeof(pthread_t) * t);
    *thread_args = (struct thread_arg*) allocateMemory(sizeof(struct thread_arg) * t);

	int poscounter = 0;

	int *m1 = (int*) allocateMemory(sizeof(int));
	int *m2 = (int*) allocateMemory(sizeof(int));
	double* maxResiduum = (double*) allocateMemory(sizeof(double));

	sem_t *maxResiduum_sem = (sem_t*) allocateMemory(sizeof(sem_t));

	if (sem_init(maxResiduum_sem, 0, 1) != 0) {
   		printf("Semaphore initialization failed.\n");
        return 1;
    }

    for(int i = 0; i < t; i++) {
		int has_remainder = t < R;
		(*thread_args)[i].start_index = poscounter;
		int work_length = L + has_remainder;
        (*thread_args)[i].work_length = work_length;
		(*thread_args)[i].m1 = m1;
		(*thread_args)[i].m2 = m2;
		(*thread_args)[i].maxResiduum = maxResiduum;
		(*thread_args)[i].maxResiduum_sem = maxResiduum_sem;
		poscounter += work_length;
		(*thread_args)[i].arguments = arguments;
		(*thread_args)[i].options = options;
    }

	return 0;
}

int calculate_new(struct options* options, struct calculation_results* results, pthread_t** threads, struct thread_arg** thread_args)
{
	int t = options->number;
	int m1 = 0;
	int m2 = 1;    
	
	while (options->term_iteration > 0)
	{
		
		*((*thread_args)[0].maxResiduum) = 0;

		for(int i = 0; i < t; i++) {
			if (pthread_create(&(*threads)[i], NULL, runThread, (void *)&(*thread_args)[i]) != 0) {
				printf("Error creating thread\n");
				return 1;
			}
		}

		printf("created threads successfully!\n");

		for(int i = 0; i < t; i++) {
			pthread_join((*threads)[i], NULL);
		}

		printf("calculation done\n");

		results->stat_iteration++;

		int ii = m1;
		m1 = m2;
		m2 = ii;
		*((*thread_args)[0].m1) = m1;
		*((*thread_args)[0].m2) = m2;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (*((*thread_args)[0].maxResiduum) < options->term_precision)
			{
				options->term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			options->term_iteration--;
		}
	}

	results->m = m2;
	
	return 0;
}

void freeThreads(pthread_t* threads, struct thread_arg* thread_args) {
	sem_destroy(thread_args[0].maxResiduum_sem);
	free(threads);
	free(thread_args[0].maxResiduum_sem);
	free(thread_args[0].maxResiduum);
	free(thread_args[0].m1);
	free(thread_args[0].m2);
    free(thread_args);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	/**
	 * fuer t threads jeweils die laenge an zu bearbeitenden Positionen
	*/
    struct thread_arg *thread_args;

    pthread_t *threads;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

    if(options.method == METH_JACOBI) {
		createThreads(&arguments, &options, &threads, &thread_args);
        gettimeofday(&start_time, NULL);
		calculate_new(&options, &results, &threads, &thread_args);
		gettimeofday(&comp_time, NULL);
	}
	else {
		gettimeofday(&start_time, NULL);
		calculate_old(&arguments, &results, &options);
		gettimeofday(&comp_time, NULL);
	}

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

    if(options.method == METH_JACOBI) {
        freeThreads(threads, thread_args);
    }
    
	freeMatrices(&arguments);

	return 0;
}
