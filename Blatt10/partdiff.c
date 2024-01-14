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

static
void
initVariablesMPI (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args, int rank, int world_size){
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	proc_args->rank = rank;
	proc_args->world_size = world_size;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t lpp = (arguments->N+1) / proc_args->world_size; // lines per process
	uint64_t lpp_pure = lpp;
	uint64_t lpp_rest = (arguments->N+1) % proc_args->world_size; // rest lines per process

	/* Wenn zu viele Prozesse, setze die Worldsize auf die nötige Anzahl runter	*/
	if (!lpp)
	{
		proc_args->world_size = lpp_rest;
	}

	lpp = lpp + (proc_args->rank < lpp_rest ? 1 : 0);
	lpp += 2; // Platz für Halolines oben und unten
	proc_args->lpp = lpp;

	uint64_t start_line = 0;
	start_line = proc_args->rank * lpp_pure + (proc_args->rank < lpp_rest ? proc_args->rank : lpp_rest);
	proc_args->start_line = start_line;
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

	if (proc_args->rank >= proc_args->world_size)
	{
		return;
	}
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

/**
 * calculateMPI_GS
 * Berechnet die Matrix mit dem Gauß-Seidel-Verfahren.
*/
static void calculateMPI_GS (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_arguments* proc_args)
{
	int i, j;           /* local variables for loops */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;
	double maxResiduum = 0.0;
	double localMaxResiduum = 0.0;
	
	uint64_t rank = proc_args->rank;
	uint64_t world_size = proc_args->world_size;
	const uint64_t invalid_rank = world_size + 1;
	uint64_t above = (proc_args->rank == 0) ? invalid_rank : proc_args->rank - 1;
	uint64_t below = (proc_args->rank == proc_args->world_size - 1) ? invalid_rank : proc_args->rank + 1;
	int lpp = proc_args->lpp;


	MPI_Request request;
	MPI_Request halo_above;
	int sent_above_once = 0;
	int first_iteration = 1;
	MPI_Request halo_below;
	int sent_below_once = 0;
	double msg[(N + 1) +2];
	double msg_buf[(N + 1) +2];

	//Diese Variable hat jeder Prozess
	double LAST_ITERATION = 0; 
	//Wird von pN benutzt, um zu signalisieren, dass die Präzision erreicht wurde
	//Wird von p0 benutzt, empfaengt von pN, ob Präzision erreicht wurde
	int N_to_0_PREC_REACHED = 0;

	int term_iteration = options->term_iteration;

	if (rank >= world_size)
	{
		printf("[%d] Ueberfluessig, zurueck zu Main\n", (int) rank);
		return;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}


	if(rank == 0 && options->termination == TERM_PREC){
		//Beginne im Hintergrund das Empfangen der Nachricht von pN
		MPI_Irecv(&N_to_0_PREC_REACHED, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, &request);
	}

	while (term_iteration > 0)
	{
		double** Matrix = arguments->Matrix[0];

		if (above != invalid_rank)
		{
			printf("[%d] Empfange von %d, %d\n", (int) rank, (int) above, (int) term_iteration);
			MPI_Recv(msg_buf, N + 1 + 1 + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			LAST_ITERATION = msg_buf[N + 1];
			maxResiduum = msg_buf[N + 2];
			Matrix[0] = msg_buf;
			if(sent_above_once){
				MPI_Wait(&halo_above, MPI_STATUS_IGNORE);
			}
			MPI_Isend(Matrix[1], N + 1, MPI_DOUBLE, above, 0, MPI_COMM_WORLD, &halo_above);
			if(!sent_above_once){
				sent_above_once = 1;
			}
			printf("[%d] Empfangen von %d, %d\n", (int) rank, (int) above, (int) term_iteration);
		}


		if(rank == 0){
			maxResiduum = 0;
		}else{
			localMaxResiduum = maxResiduum;
		}		

		/*
		Nur bei Praezisionsabbruch wird geprueft, ob Prozess 0 die Nachricht von Prozess N erhaelt,
		dass alles vorbei ist
		*/
		if(options->termination == TERM_PREC){
			/* Wenn 0, prüfe ob LAST_ITERATION gesendet wurde.*/
			if (rank == 0)
			{
				//Wir brauchen hier gar keinen Test, da wir ja nur eine Nachricht erwarten
				//MPI_Test(&request, &N_to_0_PREC_REACHED, MPI_STATUS_IGNORE);
				LAST_ITERATION = (double) N_to_0_PREC_REACHED;
			}
		}
		
		/* over all rows */
		for (i = 1; i < lpp - 1; i++)
		{
			if ( rank == world_size - 1 && i == lpp - 2)
				continue;
			
			/*Vor der Berechnung von Zeile N-1 (lpp-2), muss Zeile N (lpp-1) empfangen werden vom Prozess darunter*/
			if(i == lpp - 2 && below != invalid_rank && !first_iteration){
				MPI_Recv(Matrix[lpp - 1], N + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				first_iteration = 0;
			}

			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double) (i + proc_args->start_line - 1));
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
						residuum = Matrix[i][j] - star;
						residuum = (residuum < 0) ? -residuum : residuum;
						localMaxResiduum = (residuum < localMaxResiduum) ? localMaxResiduum : residuum;
					}

				Matrix[i][j] = star;
			}
		}
		if (below != invalid_rank)
		{
			// FIXME: Deadlock nach erster Iteration, alle warten auf 0
			//Nach der Berechnung von Zeile N-1, baue Nachricht die abgeschickt werden muss.
			memcpy(msg, Matrix[lpp - 2], (N + 1) * sizeof(double));
			msg[N + 1] = LAST_ITERATION;
			msg[N + 2] = maxResiduum;

			if(sent_below_once)
			{
				printf("[%d] Warte auf %d, %d\n", (int) rank, (int) below, (int) term_iteration);
				MPI_Wait(&halo_below, MPI_STATUS_IGNORE);
			}
			MPI_Isend(msg, N + 1 + 1 + 1, MPI_DOUBLE, below, 0, MPI_COMM_WORLD, &halo_below);
			if(!sent_below_once)
			{
				sent_below_once = 1;
			}
			printf("[%d] Gesendet an %d, %d\n", (int) rank, (int) below, (int) term_iteration);
		}

		// Setze eigenes maxResiduum auf localMaxResiduum
		maxResiduum = localMaxResiduum;

		if (rank == world_size - 1)
		{
			results->stat_iteration++;
			results->stat_precision = maxResiduum;
		}

		/* check for stopping calculation depending on termination method */		
		if (options->termination == TERM_PREC)
		{
			//Nur der letzte Prozess haelt das globale maxRes, deswegen schickt er die Abbruchsnachricht an Prozess 0
			if (maxResiduum < options->term_precision && rank == world_size - 1)
			{
				N_to_0_PREC_REACHED = 1;
				MPI_Isend(&N_to_0_PREC_REACHED, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
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
	}

	if (rank == world_size - 1)
	{
		MPI_Ssend(&results->stat_iteration, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Ssend(&results->stat_precision, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else if (!rank)
	{
		MPI_Recv(&results->stat_iteration, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&results->stat_precision, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	results->m = 0;
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
	}
	if (world_size != 1)
	{
		 // DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
		int from = proc_args.start_line;
		int to = proc_args.start_line + proc_args.lpp - 3;

		DisplayMatrix(&arguments, &results, &options, rank, world_size, from, to);
	}
	else
	{
		displayMatrix(&arguments, &results, &options);
	}
	freeMatrices(&arguments);

	MPI_Finalize();
	return 0;
}
