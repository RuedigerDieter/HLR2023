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

#include <mpi.h>

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


/**
 * individuelle Prozessvariablen
*/
struct proc_arguments
{
	int rank;
	int world_size;
	int lpp;
	int start_line;
};
/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
/*
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
*/
static
void
initVariablesMPI (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t lpp = (arguments->N+1) / proc_args->world_size; // lines per process
	uint64_t lpp_rest = (arguments->N+1) % proc_args->world_size; // rest lines per process

	uint64_t start_line = 0;
	for(uint64_t x = 0; x < proc_args->rank; x++) {
		start_line += lpp + (x < lpp_rest ? 1 : 0);
	}
	proc_args->start_line = start_line;

	proc_args->lpp = lpp + (proc_args->rank < lpp_rest ? 1 : 0);
	proc_args->lpp += 2; // everyone has 2 halolines...
	if(proc_args->rank == 0 || proc_args->rank == proc_args->world_size - 1)
	{
		--proc_args->lpp; // except first and last process
	}

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

static
void
allocateMatricesMPI (struct calculation_arguments* arguments, struct process_arguments* proc_args)
{
	uint64_t i, j;

	uint64_t const columnwidth = arguments->N+1;


	arguments->M = allocateMemory(arguments->num_matrices * proc_args->lpp * columnwidth * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(proc_args->lpp * sizeof(double*));

		for (j = 0; j < proc_args->lpp; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * proc_args->lpp * columnwidth) + (j * columnwidth);
		}
	}
}

// TODO: allocateMatricesMPI, Platz für Blöcke allocaten
// Definitiv diesmal 0er-Zeilen so machen wie bei sequentiell


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
		for(i = 0; i < N; i++)
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				Matrix[j][i][0] = 1 + (1 - (h * i)); // Linke Kante
				Matrix[j][N][i] = 1 - (h * i); // Untere Kante
				Matrix[j][N - i][N] = h * i; // Rechte Kante
				Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
			}
		}
	}
}

static
void
initMatricesMPI (struct calculation_arguments* arguments, struct options const* options, struct process_arguments* proc_args)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < proc_args->lpp; i++)
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
		for (j = 0; j < arguments->num_matrices; j++)
		{
			if(proc_args->rank == 0) //Nur proc0 hat obere Kante
			{
				for(i = 0; i < N; i++) {
					Matrix[j][0][N - i] = 1 + h * i; // Obere Kante
				}
			}

			if(proc_args->rank == proc_args->world_size - 1) //Nur procN hat untere Kante
			{
				for(i = 0; i < N; i++) {
					Matrix[j][proc_args->lpp - 1][i] = 1 - (h * i); // Untere Kante
				}
			}

			for(i = 0; i < N; i++){ //Jeder hat linke und rechte Kante
				if(i >= proc_args->start_line && i < proc_args->start_line + proc_args->lpp - 1){ //Wenn i in der Matrix liegt
					int ii = i - proc_args->start_line; //ii ist i in der Matrix
					Matrix[j][ii][0] = 1 + (1 - (h * ii)); // Linke Kante mit ii
					Matrix[j][N - ii][N] = h * ii; // Rechte Kante mit ii
				}
			}
		}
	}

	
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
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


// TODO: calculateMPI, MPI-Teil der Berechnung
// Unterscheidung p0, pn

/**
 * calculateMPI_GS
 * Berechnet die Matrix mit dem Gauß-Seidel-Verfahren.
*/
static void calculateMPI_GS (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct proc_arguments* proc_args)
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
	
	int rank = proc_args->rank;
	int world_size = proc_args->world_size;
	int above = (proc_args->rank == 0) ? NULL : proc_args->rank - 1;
	int below = (proc_args->rank == proc_args->world_size - 1) ? NULL : proc_args->rank + 1;
	int LAST_ITERATION = 0;

	MPI_Request request;
	MPI_Request halo_above;
	MPI_Request halo_below;
	int msg[(N + 1) +1];
	int msg_buf[(N + 1) +1];

	int term_iteration = options->term_iteration;

	if (rank >= world_size)
	{
		return;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix = arguments->Matrix[0];

		maxResiduum = 0;
		
		/* Wenn 0, prüfe ob LAST_ITERATION gesendet wurde.*/
		if (rank == 0)
		{
			MPI_Test(&request, &LAST_ITERATION, MPI_STATUS_IGNORE);
		}
		
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
				star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);

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

				Matrix[i][j] = star;
			}
		}

		if (rank == 0)
		{
			results->stat_iteration++;
			results->stat_precision = maxResiduum;
		}

		/* check for stopping calculation depending on termination method */
		// TODO MaxRes schicken, 		
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				int buf_LAST_ITERATION = 1;
				MPI_Isend(&buf_LAST_ITERATION, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
			}
			if (LAST_ITERATION)
			{
				term_iteration = 0;
			}
		
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
		*msg = Matrix[lpp - 1];
		msg[N + 1] = LAST_ITERATION;

		if (below != NULL)
		{
			MPI_Recv(Matrix[proc_args->lpp], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Wait(&halo_below, MPI_STATUS_IGNORE);
			MPI_Isend(msg, N + 1 + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, &halo_below);
		}
		if (above != NULL)
		{
			// evtl FIXME
			MPI_Recv(msg_buf, N + 1 + 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			LAST_ITERATION = msg_buf[N + 1];
			Matrix[0] = msg_buf;
			MPI_Wait(&halo_above, MPI_STATUS_IGNORE);
			MPI_Isend(Matrix[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, &halo_above);
		}
	}

	results->m = m2;
}
/**
 * calculateMPI_Jacobi
 * Berechnet die Matrix mit dem Jacobi-Verfahren.
 * Inspiriert von Patricks Gruppe.
*/
static void calculateMPI_Jacobi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct proc_arguments* proc_args)
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

	int rank = proc_args->rank;
	int world_size = proc_args->world_size;

	// TODO: LocalN

	/* initialize m1 and m2 depending on algorithm */
	m1 = 0;
	m2 = 1;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	int above = (proc_args->rank == 0) ? NULL : proc_args->rank - 1;
	int below = (proc_args->rank == proc_args->world_size - 1) ? NULL : proc_args->rank + 1;

	while (term_iteration > 0)
	{
		double **Matrix_Out = arguments->Matrix[m1];
		double **Matrix_In = arguments->Matrix[m2];

		localMaxResiduum = 0;

		if (rank < world_size) { // Idle threads shouldn't do anything

			/* over all rows */
			// local row
			for (i = 1; i < arguments-> - 1; i++)
			{
				double fpisin_i = 0.0;

				if (options->inf_func == FUNC_FPISIN)
				{
					fpisin_i = fpisin * sin(pih * (double) (i + arguments->from - 1));
				}

				/* over all columns */
				// local columns
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


			if (proc_args->rank % 2)
			{
				if (above != NULL)
				{
					// TODO check column width
					MPI_Ssend(Matrix_Out[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD);
					MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				if (below != NULL)
				{
					MPI_Ssend(Matrix_Out[proc_args->lpp - 2], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD);
					MPI_Recv(Matrix_Out[proc_args->lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			}
			else
			{
				if (below != NULL)
				{
					MPI_Recv(Matrix_Out[proc_args->lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Ssend(Matrix_Out[proc_args->lpp - 2], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD);
				}
				if (above != NULL)
				{
					MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Ssend(Matrix_Out[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD);
				}
			}

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

	struct proc_arguments proc_args;

	int rank, world_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	askParams(&options, argc, argv);

	proc_args.rank = rank;
	proc_args.world_size = world_size;
	
	if (!(world_size == 1))
	{
		// TODO inidividuelles N berechnen
		// TODO lpp berechnen
		// TODO Worldsize zum disqualifizieren von threads benutzen
		initVariablesMPI(&arguments, &results, &options);
		allocateMatricesMPI(&arguments);
		initMatricesMPI(&arguments, &options);
		
	}
	else
	{
		initVariables(&arguments, &results, &options);
		allocateMatrices(&arguments);
		initMatrices(&arguments, &options);
	}

	gettimeofday(&start_time, NULL);

	if (world_size == 1) 
	{
		calculate(&arguments, &results, &options);
	} 
	else 
	{
		if (options.method == METH_GAUSS_SEIDEL) 
		{
			calculateMPI_GS(&arguments, &results, &options, &proc_args);
		} 
		else 
		{
			calculateMPI_Jacobi(&arguments, &results, &options, &proc_args);
		}
	}
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	
	if (world_size != 1)
	{
		 // DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
		DisplayMatrix(&arguments, &results, &options, rank, world_size, proc_args.start_line, proc_args.start_line + proc_args.lpp - 1);
	}
	else
	{
		displayMatrix(&arguments, &results, &options);
	}
	freeMatrices(&arguments);

	MPI_Finalize();
	return 0;
}
