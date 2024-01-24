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
#include <stdbool.h>
#include <mpi.h>

#include "partdiff.h"

struct calculation_arguments
{
	uint64_t N;			   /* number of spaces between lines (lines=N+1)     */
	uint64_t num_matrices; /* number of matrices                             */

	int rank;
	int world_size;
	uint64_t from;
	uint64_t to;
	uint64_t localN;

	double h;		  /* length of a space between two lines            */
	double ***Matrix; /* index matrix used for addressing M             */
	double *M;		  /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t m;
	uint64_t stat_iteration; /* number of current iteration                    */
	double stat_precision;	 /* actual precision of all slaves in iteration    */
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
static void
initVariables(struct calculation_arguments *arguments, struct calculation_results *results, struct options const *options, int world_size, int rank)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	if ((uint64_t) world_size > arguments->N - 1)
		world_size = arguments->N - 1;

	arguments->rank = rank;
	arguments->world_size = world_size;

	int n = arguments->N - 1;
	int remainder = n % world_size;
	// From und to sind inkl., halo-Zeilen sind indizes from-1 und to+1
	arguments->from = 1 + (n / world_size) * rank + (rank < remainder ? rank : remainder);
	// arguments->to ist inklusive!
	arguments->to = n / world_size * (rank + 1) + (rank + 1 < remainder ? rank + 1 : remainder);
	arguments->localN = arguments->to - arguments->from + 3; // All lines from from to to, plus 2 halo lines

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void
freeMatrices(struct calculation_arguments *arguments)
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
static void *
allocateMemory(size_t size)
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
static void
allocateMatrices(struct calculation_arguments *arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;
	uint64_t const localN = arguments->localN;

	arguments->M = allocateMemory(arguments->num_matrices * localN * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double **));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double *));

		for (j = 0; j < localN; j++)
		{
			// These are the columns of the matrix
			// We only need localN columns, but the width is N
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * localN) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void
initMatrices(struct calculation_arguments *arguments, struct options const *options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	uint64_t const localN = arguments->localN;
	double const h = arguments->h;
	double ***Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < localN; i++)
		{
			for (j = 0; j <= N; j++)
			{
				// localN rows, N+1 columns
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{

		for (i = 0; i < localN; i++) // Linke Kante
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1 + (1 - (h * (i + arguments->from - 1)));
			}
		}

		if (arguments->rank == arguments->world_size - 1) // Untere Kante
			for (i = 0; i < N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][localN - 1][i] = 1 - (h * i);
				}
			}

		for (i = 0; i < localN; i++) // Rechte Kante
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][N] = (1 - (h * (i + arguments->from - 1)));
			}
		}

		if (arguments->rank == 0) // Obere Kante
			for (i = 0; i <= N; i++)
			{
				for (j = 0; j < arguments->num_matrices; j++)
				{
					Matrix[j][0][N - i] = 1 + h * i;
				}
			}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void
calculate_Jacobi(struct calculation_arguments const *arguments, struct calculation_results *results, struct options const *options)
{
	uint64_t i, j;								/* local variables for loops */
	int m1, m2;									/* used as indices for old and new matrices */
	double star;								/* four times center value minus 4 neigh.b values */
	double residuum;							/* residuum of current iteration */
	double localMaxResiduum, globalMaxResiduum; /* maximum residuum value of a slave in iteration */

	uint64_t const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	m1 = 0;
	m2 = 1;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double **Matrix_Out = arguments->Matrix[m1];
		double **Matrix_In = arguments->Matrix[m2];

		localMaxResiduum = 0;

		if (arguments->rank < arguments->world_size) { // Idle threads shouldn't do anything

			/* over all rows */
			for (i = 1; i < arguments->localN - 1; i++)
			{
				double fpisin_i = 0.0;

				if (options->inf_func == FUNC_FPISIN)
				{
					fpisin_i = fpisin * sin(pih * (double) (i + arguments->from - 1));
				}

				/* over all columns */
				for (j = 1; j < N; j++)
				{
					star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

					if (options->inf_func == FUNC_FPISIN)
					{
						star += fpisin_i * sin(pih * (double)j);
					}

					if (options->termination == TERM_PREC || term_iteration == 1)
					{
						residuum = Matrix_In[i][j] - star;
						residuum = (residuum < 0) ? -residuum : residuum;
						localMaxResiduum = (residuum < localMaxResiduum) ? localMaxResiduum : residuum;
					}

					Matrix_Out[i][j] = star;
				}
			}

			// Send last computed row to next process
			if (arguments->rank != arguments->world_size - 1)
				MPI_Ssend(Matrix_Out[arguments->localN - 2], N + 1, MPI_DOUBLE, arguments->rank + 1, 0, MPI_COMM_WORLD);
			// Receive first halo row from previous process
			if (arguments->rank != 0)
				MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, arguments->rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// Send first computed row to previous process
			if (arguments->rank != 0)
				MPI_Ssend(Matrix_Out[1], N + 1, MPI_DOUBLE, arguments->rank - 1, 0, MPI_COMM_WORLD);
			// Receive last halo row from next process
			if (arguments->rank != arguments->world_size - 1)
				MPI_Recv(Matrix_Out[arguments->localN - 1], N + 1, MPI_DOUBLE, arguments->rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}

		// Max Reduction over maxResiduum and communicate to all processes
		MPI_Allreduce(&localMaxResiduum, &globalMaxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

		results->stat_iteration++;
		results->stat_precision = globalMaxResiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (globalMaxResiduum < options->term_precision)
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
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static void
calculate_GS(struct calculation_arguments const *arguments, struct calculation_results *results, struct options const *options)
{
	uint64_t i, j;								/* local variables for loops */
	int m1, m2;									/* used as indices for old and new matrices */
	double star;								/* four times center value minus 4 neigh.b values */
	double residuum;							/* residuum of current iteration */
	double localMaxResiduum, globalMaxResiduum; /* maximum residuum value of a slave in iteration */

	double localMaxResiduumHistory[arguments->world_size]; // Used to keep track of maxResiduums of previous iterations

	uint64_t const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;
	bool term_precision_found = false;

	MPI_Request request1, request2;

	/* initialize m1 and m2 depending on algorithm */
	m1 = 0;
	m2 = 0;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double **Matrix_Out = arguments->Matrix[m1];
		double **Matrix_In = arguments->Matrix[m2];

		localMaxResiduum = 0;

		if (arguments->rank < arguments->world_size) { // Idle threads shouldn't do anything

			/* over all rows */
			for (i = 1; i < arguments->localN - 1; i++)
			{
				double fpisin_i = 0.0;

				if (options->inf_func == FUNC_FPISIN)
				{
					fpisin_i = fpisin * sin(pih * (double) (i + arguments->from - 1));
				}

				// Receive halo rows right before they are needed (except in the first iteration)
				// Also wait for the sending of the previous iteration to complete

				// Receive first halo row from previous process
				if (i == 1 && arguments->rank != 0) {
					MPI_Recv(Matrix_In[0], N + 1, MPI_DOUBLE, arguments->rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					// If we have already sent something in the previous iteration, wait for it now
					if (results->stat_iteration > 0)
						MPI_Wait(&request1, MPI_STATUS_IGNORE);
				}
				// Receive last halo row from next process
				if (i == arguments->localN - 2 && arguments->rank != arguments->world_size - 1 && results->stat_iteration > 0) {
					MPI_Recv(Matrix_In[arguments->localN - 1], N + 1, MPI_DOUBLE, arguments->rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Wait(&request2, MPI_STATUS_IGNORE);
				}

				/* over all columns */
				for (j = 1; j < N; j++)
				{
					star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

					if (options->inf_func == FUNC_FPISIN)
					{
						star += fpisin_i * sin(pih * (double)j);
					}

					if (options->termination == TERM_PREC || term_iteration == 1)
					{
						residuum = Matrix_In[i][j] - star;
						residuum = (residuum < 0) ? -residuum : residuum;
						localMaxResiduum = (residuum < localMaxResiduum) ? localMaxResiduum : residuum;
					}

					Matrix_Out[i][j] = star;
				}

				// Send first & last rows immediately after they are calculated

				// Send first row to previous process (unless that process has terminated)
				// The previous process is always one iteration ahead of this process, and term_iteration will be set even when using term_precision
				if (i == 1 && arguments->rank != 0 && term_iteration > 1) {
					MPI_Issend(Matrix_Out[i], N + 1, MPI_DOUBLE, arguments->rank - 1, 0, MPI_COMM_WORLD, &request1);
				}
				// Send last row to next process
				if (i == arguments->localN - 2 && arguments->rank != arguments->world_size - 1) {
					MPI_Issend(Matrix_Out[arguments->localN - 2], N + 1, MPI_DOUBLE, arguments->rank + 1, 0, MPI_COMM_WORLD, &request2);
				}

			}

		}

		results->stat_iteration++;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC && !term_precision_found)
		{
			
			// Shift localMaxResiduumHistory back
			for (int i = arguments->world_size - 2; i >= 0; i--) {
				localMaxResiduumHistory[i + 1] = localMaxResiduumHistory[i];
			}
			localMaxResiduumHistory[0] = localMaxResiduum;

			// If every process has done at least one iteration, reduce over the max residuum of the last process (since it's in the lowest iteration)
			int correctedRank = arguments->rank >= arguments->world_size ? 0 : arguments->rank; // For this reduction, idle processes will be treated like the first process
			int lastProcessIteration = -arguments->world_size + correctedRank + results->stat_iteration;
			if (lastProcessIteration >= 0) {
				// What was this process's localMaxResiduum when it was in the iteration of the last process?
				double localMaxResiduumFromLastProcessIteration = localMaxResiduumHistory[arguments->world_size - correctedRank - 1];
				MPI_Allreduce(&localMaxResiduumFromLastProcessIteration, &globalMaxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				results->stat_precision = globalMaxResiduum;

				if (globalMaxResiduum < options->term_precision)
				{
					// Stop when we reach the iteration that process 0 is currently in
					term_iteration = correctedRank;
					term_precision_found = true;
				}
			}
		}
		else if (options->termination == TERM_ITER || term_precision_found)
		{
			term_iteration--;
		}
	}

	// When done, still wait for the send to the next process to finish
	if (arguments->rank < arguments->world_size - 1)
		MPI_Wait(&request2, MPI_STATUS_IGNORE);

	MPI_Allreduce(&localMaxResiduum, &globalMaxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	results->stat_precision = globalMaxResiduum;

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void
displayStatistics(struct calculation_arguments const *arguments, struct calculation_results const *results, struct options const *options)
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
	printf("Interlines:         %" PRIu64 "\n", options->interlines);
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

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static void
displayMatrix_MPI(struct calculation_arguments *arguments, struct calculation_results *results, struct options *options)
{
	int rank = arguments->rank;
	int from = arguments->from;
	int to = arguments->to;
	int size = arguments->world_size;
	int const elements = arguments->N + 1;

	int x, y;
	double **Matrix = arguments->Matrix[results->m];
	MPI_Status status;

	/* first line belongs to rank 0 */
	if (rank == 0)
		from--;

	/* last line belongs to rank size - 1 */
	if (rank + 1 == size)
		to++;

	if (rank == 0)
		printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		int line = y * (options->interlines + 1);

		if (rank == 0)
		{
			/* check whether this line belongs to rank 0 */
			if (line < from || line > to)
			{
				/* use the tag to receive the lines in the correct order
				 * the line is stored in Matrix[0], because we do not need it anymore */
				MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
			}
		}
		else
		{
			if (line >= from && line <= to)
			{
				/* if the line belongs to this process, send it to rank 0
				 * (line - from + 1) is used to calculate the correct local address */
				MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
			}
		}

		if (rank == 0)
		{
			for (x = 0; x < 9; x++)
			{
				int col = x * (options->interlines + 1);

				if (line >= from && line <= to)
				{
					/* this line belongs to rank 0 */
					printf("%7.4f", Matrix[line][col]);
				}
				else
				{
					/* this line belongs to another rank and was received above */
					printf("%7.4f", Matrix[0][col]);
				}
			}

			printf("\n");
		}
	}

	fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main(int argc, char **argv)
{
	MPI_Init(NULL, NULL);
	int world_size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	if (rank == 0)
	{
		askParams(&options, argc, argv);
	}
	// Create options struct in MPI
	int blocklengths_options[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Datatype types_options[7] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE};
	MPI_Aint offsets_options[7] = {offsetof(struct options, number), offsetof(struct options, method), offsetof(struct options, interlines), offsetof(struct options, inf_func), offsetof(struct options, termination), offsetof(struct options, term_iteration), offsetof(struct options, term_precision)};
	MPI_Datatype MPI_options_type;
	MPI_Type_create_struct(7, blocklengths_options, offsets_options, types_options, &MPI_options_type);
	MPI_Type_commit(&MPI_options_type);

	// Broadcast options struct to all processes
	MPI_Bcast(&options, 1, MPI_options_type, 0, MPI_COMM_WORLD);
	initVariables(&arguments, &results, &options, world_size, rank);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	if (options.method == METH_GAUSS_SEIDEL)
		calculate_GS(&arguments, &results, &options);
	else
		calculate_Jacobi(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	if (rank == 0)
		displayStatistics(&arguments, &results, &options);
	displayMatrix_MPI(&arguments, &results, &options);

	freeMatrices(&arguments);

	// Free MPI types
	MPI_Type_free(&MPI_options_type);

	MPI_Finalize();
	return 0;
}
