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
#include <sys/time.h>

#include "partdiff.h"
#include <pthread.h>
#include <limits.h>

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

struct t_args 
{
	struct t_data* t_data;
	int range_start,range_end;
};
struct t_data
{
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	int lock;
	uint64_t N;
	double** Matrix_In;
	double** Matrix_Out;
	int lock_leave_count;
	
	double pih,fpisin;
	double maxResiduum;
	const struct options* options;
	const struct calculation_arguments* arguments;
	struct calculation_results* results;
};

void* t_calculate (void * args)
{
	struct t_args* t_args = (struct t_args*) args;
	struct t_data* t_data = (struct t_data*) t_args->t_data;

	const struct options* options = t_data->options;
	const struct calculation_arguments* arguments = t_data->arguments;
	struct calculation_results* results = t_data->results;

	int i = 0;
	int j = 0;           
	int m1 = 0;
	int m2 = 1;

	int const N = t_data->N;
	double const pih = t_data->pih;
	double const fpisin = t_data->fpisin;

	double residuum = 0;
	double local_maxResiduum = 0;
	double star = 0;
	int term_iteration = options->term_iteration;
	int thread_count = (int) options->number;

	double fpisin_i = 0.0;

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		/* over all rows */
		for (i = t_args->range_start; i < t_args->range_end; i++)
		{
			fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j <  N; j++)
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
					local_maxResiduum = (residuum < local_maxResiduum) ? local_maxResiduum : residuum;
				}
				Matrix_Out[i][j] = star;
			}

		}

		pthread_mutex_lock(&t_data->mutex);
		t_data->maxResiduum = (local_maxResiduum < t_data->maxResiduum) ? t_data->maxResiduum : local_maxResiduum;		

		if (t_data->lock_leave_count == 0)
		{
			pthread_cond_broadcast(&t_data->cond);

		}
		else
		{
			while (t_data->lock_leave_count > 0)
			{
				pthread_cond_wait(&t_data->cond, &t_data->mutex);
			}
		}
		
		results->stat_iteration++;
		
		if (results->stat_iteration % thread_count == 0)
		{
			t_data->lock_leave_count = thread_count;
			pthread_cond_broadcast(&t_data->cond);

		}
		else
		{
			while (results->stat_iteration % thread_count != 0)
			{
				pthread_cond_wait(&t_data->cond, &t_data->mutex);
			}
		}
		t_data->lock_leave_count--;
		pthread_mutex_unlock(&t_data->mutex);

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (t_data->maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
	pthread_exit(NULL);
}

static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
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

	int thread_count = 0;
	pthread_t* threads = NULL;

	struct t_args* t_args = NULL;
	struct t_data* t_data = NULL;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;

		thread_count = (int) options->number;
		threads = malloc(sizeof(pthread_t) * thread_count);

		t_args = malloc(sizeof(struct t_args) * thread_count);
		t_data = malloc(sizeof(struct t_data));

		t_data->N = N;
		pthread_mutex_init(&t_data->mutex, NULL);
		pthread_cond_init(&t_data->cond, NULL);	
		t_data->maxResiduum = 0;
		t_data->arguments = arguments;
		t_data->results = results;
		t_data->options = options;
		t_data->lock_leave_count = 0;

		int rpt = (N - 1)  / thread_count;
		int rst = (N - 1) % thread_count;
		int start = 1;
		for (i = 0; i < thread_count; i++)
		{
			t_args[i].t_data = t_data;
			t_args[i].range_start = start;
			t_args[i].range_end = start + rpt + (rst > 0 ? 1 : 0);
			start = t_args[i].range_end;
			rst--;
			printf("Finished setting up Thread %d\n", i);
		}
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
		if (options->method == METH_JACOBI)
		{
			t_data->pih = pih;
			t_data->fpisin = fpisin;
		}
	}

	if (options->method == METH_JACOBI)
	{
		printf("Starting Threads...\n");
		for (i = 0; i < thread_count; i++)
		{
			printf("Started Thread %d\n", i);
			pthread_create(&threads[i], NULL, t_calculate, &t_args[i]);
		}
		t_data->lock = 0;
		printf("Waiting for Threads...\n");
		for (i = 0; i < thread_count; i++)
		{
			pthread_join(threads[i], NULL);
			results->stat_precision = t_data->maxResiduum;
		}
		pthread_mutex_destroy(&t_data->mutex);
		free(t_args);
		free(t_data);
		free(threads);
		results->stat_iteration = results->stat_iteration / thread_count;
	}
	else
	{
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

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
