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

struct t_data {
	int num_threads;						// Anzahl der Threads, vllt unnötig
	pthread_t* threads;						// Array für die Threads
	int N;									// N
	const struct calculation_arguments* arguments;
	const struct calculation_results* results;
	const struct options* options;
	double** Matrix_In;
	double** Matrix_Out;

	double pih;
	double fpisin;

	int term_iteration;
};

typedef struct  {
	int position;							// Startindex pro thread
	int chunksize;							// Größe des zu berechnenden Bereichs
	int lock;								// Sperre für die Threads
	struct t_data* t_data;					// globale Variablen
}T_args;

void t_calculate(void* args) 
{
	int pos_r, pos_l, pos_u, pos_d; 		//Positionen für den Stern
	int zeile,spalte;
	int star = 0;
	int residuum = 0;
	int maxResiduum = 0;

	T_args* t_args = (T_args*) args;
	// double** Matrix_In = t_args->t_data->Matrix_In;
	// double** Matrix_Out = t_args->t_data->Matrix_Out;

	while (1)
	{
		if (!t_args->lock)
		{
			for (int i = t_args->position; i < (t_args->position + t_args->chunksize); i++)
			{
				if(i % t_args->t_data->N == 0 
				|| i % t_args->t_data->N == t_args->t_data->N - 1) 
				continue; // Ignoriere Randzellen

				pos_r = i + 1;
				pos_l = i - 1;
				pos_u = i - t_args->t_data->N;
				pos_d = i + t_args->t_data->N;

				star = .25 * (*t_args->t_data->Matrix_In[pos_r] 
							+ *t_args->t_data->Matrix_In[pos_l] 
							+ *t_args->t_data->Matrix_In[pos_u] 
							+ *t_args->t_data->Matrix_In[pos_d]);


				if (t_args->t_data->options->inf_func == FUNC_FPISIN)
				{
					zeile = i / t_args->t_data->N;
					spalte = i % t_args->t_data->N;


					star += t_args->t_data->fpisin * t_args->t_data->fpisin 
					* sin(t_args->t_data->pih * (double) zeile) 
					* sin(t_args->t_data->pih * (double) spalte);
				}

				if (t_args->t_data->options->termination == TERM_PREC || t_args->t_data->term_iteration == 1)
				{
					residuum = *t_args->t_data->Matrix_In[i] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
				}

				*t_args->t_data->Matrix_Out[i] = star;

			}
			t_args->lock = 1; 
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct t_data* t_data)
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

	T_args* t_args;
	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;

		t_data = malloc(sizeof(struct t_data));
		t_args = malloc(sizeof(T_args) * options->number);

		t_data->num_threads = (int) options->number;	// Anzahl der Threads
		t_data->arguments = arguments;			
		t_data->results = results;
		t_data->options = options;

		t_data->threads = malloc(sizeof(pthread_t) * options->number);	// Array für die Threads 
		
		t_data->N = N;
		int M = (N - 1) * (N - 1);				// Anzahl der Zellen


		int cpt = M / options->number;				// Cells pro Thread
		int cpt_rest = M % options->number;			// Rest Cells pro Thread

		int position;

		for (i = 0; i < t_data->num_threads; i++)		// setup für alle threads
		{
			t_args[i].position = position;			// setzt die Startposition für die Matrix und verteilt dabei den Rest, damit alle Zellen sequentiell zugewiesen sind.
			if (cpt_rest > 0)
			{
				t_args[i].chunksize = cpt + 1;
				cpt_rest--;
				position += cpt + 1;
			}
			else
			{
				t_args[i].chunksize = cpt;
				position += cpt;
			}
			t_args[i].lock = 1;						// Rechensperre
			pthread_create(&t_data->threads[i], NULL, t_calculate, (void*) t_args);
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

		t_data->pih = pih;
		t_data->fpisin = fpisin;
	}


	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxResiduum = 0;

		if (options->method == METH_JACOBI) // Jacobi
		{
			t_data->term_iteration = term_iteration;
			for (int i = 0; i < options->number; i++)
			{
				t_args[i].lock = 0; // Rechensperre aufheben
			}
		}
		else // Gauß-Seidel
		{
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

	if (options->method == METH_JACOBI) //freeing memory for Threading
	{
		free(t_data->threads);
		free(t_args);
	}

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
	struct t_data t_data;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options,&t_data);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
