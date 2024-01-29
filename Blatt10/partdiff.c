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
//musste machen weil auf mac gabs aerger!
#ifdef __linux__
#include <malloc.h>
#endif
#include <sys/time.h>

#include "partdiff.h"

#include <mpi.h>
#include <string.h>
#include <limits.h>



/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/**
 * individuelle Prozessvariablen
*/
struct process_arguments
{
	uint64_t rank;
	uint64_t world_size;
	uint64_t lpp;
	uint64_t start_line;
};

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

static void zeroProcArgs(struct process_arguments* proc_args){
	proc_args->rank = 0;
	proc_args->world_size = 0;
	proc_args->lpp = 0;
	proc_args->start_line = 0;
}


/**
 * initVariablesMPI: Initialisiere die globalen und Prozessspezifischen Variablen. 
 * Rechnet u.a. aus, wie viele Zeilen jeder einzelne Prozess bearbeiten muss.
*/
static
void
initVariablesMPI (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args, uint64_t rank, int world_size){
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	proc_args->rank = rank;
	proc_args->world_size = world_size;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t usable_lines = ((arguments->N)+1-2);

	uint64_t lpp = (usable_lines) / proc_args->world_size; // lines per process

	/* Behandelt den Sonderfall, bei dem zu viele Prozesse für die gegebenen Interlines existieren.*/
	if (!lpp)
	{
		proc_args->world_size = usable_lines;
		proc_args->lpp = 3;
		proc_args->start_line = rank+1;

	}
	else
	{
		uint64_t lpp_pure = lpp;
		uint64_t lpp_rest = (usable_lines) % proc_args->world_size; // rest lines per process

		lpp = lpp + (proc_args->rank < lpp_rest ? 1 : 0);
		lpp += 2; // Platz für Halolines oben und unten
		
		proc_args->lpp = lpp;
		
		uint64_t start_line_rest = (rank < lpp_rest) ? rank : lpp_rest;
		uint64_t start_line = rank * lpp_pure + start_line_rest;
		start_line += 1; // Platz für Haloline oben
		
		proc_args->start_line = start_line;
	}

	// printf("[%d] Handling ll %d (lpp: %d)\n", (int) rank, (int) proc_args->start_line, (int) proc_args->lpp);
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


/**
 * allocateMatricesMPI: Reserviere Prozess-Matrizen, die später die Teilstücke der Gesamtmatrix beinhalten.
 * Jede Matrix ist lpp Zeilen (meistens chunk + 2 Halolines) "hoch" und N + 1 Spalten breit.
*/
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


/**
 * initMatricesMPI: Initialisiere die Prozess-Matrizen.
 * mit FPISIN nur mit 0,
 * ohne FPISIN zusätzlich mit Rändern.
*/
static
void
initMatricesMPI (struct calculation_arguments* arguments, struct options const* options, struct process_arguments* proc_args)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	uint64_t world_size = proc_args->world_size;
	uint64_t rank = proc_args->rank;
	uint64_t lpp = proc_args->lpp;
	uint64_t start_line = proc_args->start_line; 


	if (proc_args->rank >= proc_args->world_size)
		return;

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
			
			/* Initialisiere Untere Kante */
			if(!rank)
			{
				for(i = 0; i <= N; i++) {
					Matrix[j][0][N - i] = 1 + h * i;
				}
			}

			/* Initialisiere untere Kante */
			if(rank == world_size - 1) 
			{
				for(i = 0; i <= N; i++) {
					Matrix[j][proc_args->lpp - 1][i] = 1 - (h * i);
				}
			}

			/* Initialisiere seitliche Kanten.*/
			for (i = 0; i < lpp - 1; i++)
			{
				int global_pos = (rank > 0) ? start_line + i - 1: i;
				Matrix[j][i][0] = 1 + (1 - (h * global_pos)); // Linke Kante
			}

			for (i = lpp - 2; i > 0; i--)
			 {
				int global_pos = start_line + i - 1;
				global_pos = N - global_pos;
			 	Matrix[j][i][N] = h * global_pos; // Rechte Kante
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


static void calculateMPI_GS (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args)
{
	/* Deklariere und Initialisiere Variablen */
	/* Rechnungs-Variablen */

	int i, j;           /* local variables for loops */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */

	int const N = arguments->N;
	double const h = arguments->h;
	int term_iteration = options->term_iteration;

	double pih = 0.0;
	double fpisin = 0.0;
	double maxResiduum = 0.0;
	double localMaxResiduum = 0.0;

	double** Matrix = arguments->Matrix[0];


	/* MPI-Variablen */
	const uint64_t rank = proc_args->rank;
	const uint64_t world_size = proc_args->world_size;
	int lpp = proc_args->lpp;

	const int LAST_RANK = rank == world_size - 1;
	const int FIRST_RANK = rank == 0;

	/* Kommunikations-Variablen */
	const uint64_t invalid_rank = world_size + 1;
	uint64_t above = FIRST_RANK ? invalid_rank : proc_args->rank - 1;
	uint64_t below = LAST_RANK ? invalid_rank : proc_args->rank + 1;

	MPI_Request PREC_TERM;
	MPI_Request HALO_A, HALO_B;
	int SENT_A = 1;
	int SENT_B = 1;

	double s_buf_halo_prec[N + 1 + 2];
	double r_buf_halo_prec[N + 1 + 2];

	const double s_LAST_ITERATION = 1;
	double r_LAST_ITERATION;

	/* Breche überflüssige Prozesse ab */
	if (rank >= world_size)
	{
		printf("[%d] Ueberfluessiger Prozess, zurueck zu Main\n", (int) rank);
		return;
	}

	/* FPISIN-Berechnung */
	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
	
	/* PREC-TERM empfangen */
	if (FIRST_RANK && options->termination == TERM_PREC)
	{
		printf("Empfange PREC_TERM\n");
		MPI_Irecv(&r_LAST_ITERATION, 1, MPI_DOUBLE, LAST_RANK, 0, MPI_COMM_WORLD, &PREC_TERM);
	}

	/* Star-Berechnung */
	while (term_iteration > 0)
	{
		/**
		 *	Wenn erste Zeile:
		 *		Empfange Haloline von oben
		 *		Empfange MaxRes, Abbruch
		 *		Setze MaxResiduum zurück
		 *	Wenn letzte Zeile:
		 *		Sende Haloline nach unten
		 * 		Sende MaxRes nach unten
		 */

		 for (i = 1; i < lpp - 1; i++)
		 {
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double) (i + proc_args->start_line - 1));
			}

			if (i == 1)
			{
				if (FIRST_RANK)
					maxResiduum = 0;
				else
				{
					printf("[%d] Empfange Haloline von oben\n", (int) rank);
					/* Empfange Haloline und MaxRes von oben */
					MPI_Recv(r_buf_halo_prec, N + 1 + 2, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					r_LAST_ITERATION = r_buf_halo_prec[N + 1];
					maxResiduum = r_buf_halo_prec[N + 2];
					memcpy(Matrix[0], r_buf_halo_prec, (N + 1) * sizeof(double));
					localMaxResiduum = maxResiduum;
				}	
			}
			else if (i == lpp - 1)
			{
				MPI_Recv(Matrix[lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double) j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					localMaxResiduum = (residuum < localMaxResiduum) ? localMaxResiduum : residuum;
				}

				Matrix[i][j] = star;
			}

			if (i == 1)
			{
				if (!FIRST_RANK)
				{
					printf("[%d] Sende obere Haloline nach oben\n", (int) rank);
					if (!SENT_A)
					{
						MPI_Wait(&HALO_A, MPI_STATUS_IGNORE);
					}
					MPI_Issend(Matrix[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, &HALO_A);
					SENT_A = 0;
				}
			}
			if (i == lpp - 2)
			{
				if (LAST_RANK)
				{
					results->stat_iteration++;
					results->stat_precision = localMaxResiduum;
				}
				else
				{
					/* Baue Nachricht für Versenden nach unten */
					memcpy(s_buf_halo_prec, Matrix[lpp - 2], (N + 1) * sizeof(double));
					s_buf_halo_prec[N + 1] = r_LAST_ITERATION;
					s_buf_halo_prec[N + 2] = localMaxResiduum;

					printf("[%d] Sende untere Haloline nach unten\n", (int) rank);
					/* Sende Haloline und MaxRes nach unten */
					if (!SENT_B)
					{
						MPI_Wait(&HALO_B, MPI_STATUS_IGNORE);
					}
					MPI_Issend(s_buf_halo_prec, N + 1 + 2, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, &HALO_B);
					SENT_B = 0;
					printf("[%d] Untere Haloline nach unten gesendet\n", (int) rank);
				}
			}
		}

		if (options->termination == TERM_PREC)
		{
			if (FIRST_RANK)
			{
				printf("[%d] Prüfe ob Präzision erreicht wurde\n", (int) rank);
				int TERM_SENT = 0;
				MPI_Test(&PREC_TERM, &TERM_SENT, MPI_STATUS_IGNORE);
				if (TERM_SENT)
				{
					MPI_Wait(&PREC_TERM, MPI_STATUS_IGNORE);
					term_iteration = 0;
					r_LAST_ITERATION = 1;
				}	
			}
			else if (LAST_RANK && maxResiduum < options->term_precision)
			{
				MPI_Issend(&s_LAST_ITERATION, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &PREC_TERM);
				MPI_Wait(&PREC_TERM, MPI_STATUS_IGNORE);
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
	printf("[%d] Berechnung beendet\n", (int) rank);
	
	/* Sende Statistik an FIRST_RANK */
	if (LAST_RANK)
	{
		MPI_Ssend(&results->stat_iteration, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		MPI_Ssend(&results->stat_precision, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
	}
	else if (FIRST_RANK)
	{
		MPI_Recv(&results->stat_iteration, 1, MPI_INT, world_size - 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&results->stat_precision, 1, MPI_DOUBLE, world_size - 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

}

/**
 * calculateMPI_Jacobi
 * Berechnet die Matrix mit dem Jacobi-Verfahren.
 * Inspiriert von Patricks Gruppe.
*/
static void calculateMPI_Jacobi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args)
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


	uint64_t rank = proc_args->rank;
	uint64_t world_size = proc_args->world_size;
	uint64_t lpp = proc_args->lpp;

	const uint64_t invalid_rank = world_size + 1;

	uint64_t above = (rank == 0) ? invalid_rank : rank - 1;
	uint64_t below = (rank == world_size - 1) ? invalid_rank : rank + 1;

	if (rank >= world_size)
	{
		return;
	}

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

		/* over all rows */
		// local row
		for (i = 1; i < lpp - 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double) (i + proc_args->start_line - 1));
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
			if (above != invalid_rank)
			{
				MPI_Ssend(Matrix_Out[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD);
				MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			if (below != invalid_rank)
			{
				MPI_Ssend(Matrix_Out[proc_args->lpp - 2], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD);
				MPI_Recv(Matrix_Out[proc_args->lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		else
		{
			if (below != invalid_rank)
			{
				MPI_Recv(Matrix_Out[proc_args->lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Ssend(Matrix_Out[proc_args->lpp - 2], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD);
			}
			if (above != invalid_rank)
			{
				MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Ssend(Matrix_Out[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD);
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
		MPI_Barrier(MPI_COMM_WORLD);
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
	fflush(stdout);
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

/**
 * aus displaymatrix-mpi.c für bessere Übersichtlichkeit hierher kopiert
*/
static void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  if (rank > size)
	return;


  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
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
        MPI_Ssend(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
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



/** 
 * Displays the whole Matrix. 
 * Do not use with large Matrices.
 * Mainly used for debugging.
 * vong Sebastian
 */
// static void DisplayWholeMatrix(struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size) {
  
//   int num_columns = arguments->N + 1;
//   // Loop over all Ranks
//   for (int i = 0; i < size; i++) {
//     // Only Rank 0 should print its first line
//     if (rank == 0 && i == 0) {
//       printf("Matrix:\n");
      
//       for (int j = 0; j < num_columns; j++) {
//         printf("%7.9f ", arguments->Matrix[results->m][0][j]);
//       }

//       printf("\n");
//     }

//     // Print the main part of the matrix, without the halo lines
//     // But only, if this is the current rank
//     if (rank == i) {
//       for (int j = 1; j < (int) arguments->num_rows - 1; j++) {
//         for (int k = 0; k < num_columns; k++) {
//           printf("%7.9f ", arguments->Matrix[results->m][j][k]);
//         }
//         printf("\n");
//       }
//     }

//     // Only Rank size - 1 should print its last line
//     if (rank == size - 1 && i == size - 1) {
//       for (int j = 0; j < num_columns; j++) {
//         printf("%7.9f ", arguments->Matrix[results->m][arguments->num_rows - 1][j]);
//       }
//       printf("\n");
//     }

//     // Barrier to make sure that the output is not mixed up
//     MPI_Barrier(MPI_COMM_WORLD);
//   }
// }


/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	int rank, world_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	askParams(&options, argc, argv, rank);

	struct process_arguments proc_args;	
	
	zeroProcArgs(&proc_args);
	
	if (world_size != 1)
	{
		initVariablesMPI(&arguments, &results, &options, &proc_args, rank, world_size);
		allocateMatricesMPI(&arguments,&proc_args);
		initMatricesMPI(&arguments, &options, &proc_args);
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

	if (!rank)
	{
		displayStatistics(&arguments, &results, &options);
		printf("MPI: %d Prozesse\n", world_size);
	}
	if (world_size != 1)
	{
		int from = proc_args.start_line;
		int to = proc_args.start_line + proc_args.lpp - 3;

		DisplayMatrix(&arguments, &results, &options, rank, proc_args.world_size, from, to);
		//DisplayWholeMatrix(&arguments, &results, &options, rank, world_size);
	}
	else
	{
		displayMatrix(&arguments, &results, &options);
	}
	freeMatrices(&arguments);

	MPI_Finalize();
	return 0;
}
