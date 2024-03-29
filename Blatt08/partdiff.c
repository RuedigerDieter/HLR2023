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

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "partdiff.h"

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

struct process_args
{
	uint64_t rank;
	uint64_t world_size;
	uint64_t starting_line;
	uint64_t working_lines;
	uint64_t working_columns;
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


/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatricesMPI (struct calculation_arguments* arguments, struct process_args* proc_args)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	uint64_t lpp = N / proc_args->world_size; // lines per process
	uint64_t lpp_rest = N % proc_args->world_size; // rest lines per process

	proc_args->working_lines = lpp + (proc_args->rank < lpp_rest ? 1 : 0);

	uint64_t starting_line = 0;
	for(uint64_t x = 0; x < proc_args->rank; x++) {
		starting_line += lpp + (x < lpp_rest ? 1 : 0);
	}
	proc_args->starting_line = starting_line;

	proc_args->working_columns = arguments->N-1;

	arguments->M = allocateMemory(arguments->num_matrices * proc_args->working_lines * proc_args->working_columns * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(proc_args->working_lines * sizeof(double*));

		for (j = 0; j <= proc_args->working_lines; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * proc_args->working_lines * proc_args->working_columns) + (j * proc_args->working_columns);
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatricesMPI (struct calculation_arguments* arguments, struct process_args* proc_args)
{
	uint64_t g, i, j; /* local variables for loops */

	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i < proc_args->working_lines; i++)
		{
			for (j = 0; j < proc_args->working_columns; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}
}


static
void
calculateMPI (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, struct process_args* proc_args)
{
	int i, j;           /* local variables for loops */
	int m1 = 0;		 	//direkt gesetzt weil jacobi 
	int m2 = 1;         /* used as indices for old and new matrices */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	double global_maxResiduum = 0;

	/**
	 * haloline_in_top:		Die Zeile der Matrix, die über dem zu berechnenden Bereich liegt.
	 * haloline_in_bottom:	Die Zeile der Matrix, die unter dem zu berechnenden Bereich liegt.
	*/
	double haloline_in_top[proc_args->working_columns];
	double haloline_in_bottom[proc_args->working_columns];
	
	// printf("Rank %d working on lines %d to %d\n", proc_args->rank, proc_args->starting_line, proc_args->starting_line + proc_args->working_lines - 1);
	#ifdef _OPENMP
	#pragma omp parallel num_threads(options->number) default(none) private(star, residuum, i, j) shared(haloline_in_top, haloline_in_bottom, arguments, options, results, term_iteration, maxResiduum, global_maxResiduum, proc_args, h,  m1, m2, pih, fpisin, ompi_mpi_double, ompi_mpi_op_max, ompi_mpi_comm_world)
	{
	#endif

		/**
		 * define fpisin_i here, since it is used in the inner loop.
		*/
		double fpisin_i = 0.0;
		while (term_iteration > 0)
		{
			double** Matrix_Out = arguments->Matrix[m1];
			double** Matrix_In  = arguments->Matrix[m2];
			

			/**
			 * haloline_out_top:	Die Zeile der Matrix, die dem Prozess davor geschickt werden muss.
			 * haloline_out_bottom:	Die Zeile der Matrix, die dem Prozess danach geschickt werden muss.
			*/	
			
			double *haloline_out_top = Matrix_In[0];
			double *haloline_out_bottom = Matrix_In[proc_args->working_lines-1];

			maxResiduum = 0;
			

			/* over all rows */
			#ifdef _OPENMP
			#pragma omp for collapse(2) reduction(max:maxResiduum)
			#endif
			for (i = 0; i < (int) proc_args->working_lines; i++)
			{
				/* over all columns */
				for (j = 0; j < (int) proc_args->working_columns; j++)
				{
					/**
					 * fpisin_i needs to be set in the inner loop, to allow collapsing the 2 for loops.
					 * outer loop is the same as when the inner loop has running variable j=0
					*/
					if(j == 0){
						fpisin_i = 0.0;
						if (options->inf_func == FUNC_FPISIN)
						{
							int translated_i = i + proc_args->starting_line + 1; //+1 for the 0 line at the beginning
							fpisin_i = fpisin * sin(pih * (double)(translated_i+1)); //fixed offset
						}
					}

					#ifdef _OPENMP
					#pragma omp critical
					{
					#endif
						if(i == 0) {
							/**
							 * Kommuniziere die obere Haloline, wobei der erste Prozess hier deadlock brechen soll.
							*/
							if(proc_args->rank == 0) 
							{
								MPI_Ssend(haloline_out_bottom, proc_args->working_columns, MPI_DOUBLE, proc_args->rank+1, 10, MPI_COMM_WORLD);
							}
							else if (proc_args->rank == proc_args->world_size -1)
							{
								MPI_Recv(haloline_in_top, proc_args->working_columns, MPI_DOUBLE, proc_args->rank-1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							else 
							{
								MPI_Recv(haloline_in_top, proc_args->working_columns, MPI_DOUBLE, proc_args->rank-1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Ssend(haloline_out_bottom, proc_args->working_columns, MPI_DOUBLE, proc_args->rank+1, 10, MPI_COMM_WORLD);
							}

							/**
							 * Kommuniziere die untere Haloline, wobei der letzte Prozess hier deadlock brechen soll.
							*/

							if(proc_args->rank == 0) 
							{
								MPI_Recv(haloline_in_bottom, proc_args->working_columns, MPI_DOUBLE, proc_args->rank + 1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							}
							else if (proc_args->rank == proc_args->world_size -1)
							{
								MPI_Ssend(haloline_out_top, proc_args->working_columns, MPI_DOUBLE, proc_args->rank - 1, 10, MPI_COMM_WORLD);
							}
							else 
							{
								MPI_Recv(haloline_in_bottom, proc_args->working_columns, MPI_DOUBLE, proc_args->rank + 1, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
								MPI_Ssend(haloline_out_top, proc_args->working_columns, MPI_DOUBLE, proc_args->rank - 1, 10, MPI_COMM_WORLD);
							}
						}
					#ifdef _OPENMP
					}
					#endif

					/**
					 * Belege die Werte für a,b,c,d mit den Werten aus der Matrix.
					 * 0 wenn Randzeilen/-spalten
					*/
					double a;
					if(proc_args->rank != 0) {
						a = (i == 0) ? haloline_in_top[j] : Matrix_In[i-1][j];
					}
					else {
						a = (i == 0) ? 0 : Matrix_In[i-1][j];
					}

					double b = (j == 0) ? 0 : Matrix_In[i][j-1];
					double c = (j == (int) proc_args->working_columns-1) ? 0 : Matrix_In[i][j+1];

					double d;
					if(proc_args->rank != proc_args->world_size -1) {
						d = (i == (int) proc_args->working_lines-1) ? haloline_in_bottom[j] : Matrix_In[i+1][j];
					}
					else {
						d = (i == (int) proc_args->working_lines-1) ? 0 : Matrix_In[i+1][j];
					}

					star = 0.25 * (a + b + c + d);

					if (options->inf_func == FUNC_FPISIN)
					{
						star += fpisin_i * sin(pih * (double)(j + 1)); //fixed offset
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

			#ifdef _OPENMP
			#pragma omp single
			{
			#endif

				MPI_Allreduce(&maxResiduum, &global_maxResiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				


				results->stat_iteration++;
				results->stat_precision = global_maxResiduum;

				/* exchange m1 and m2 */
				i = m1;
				m1 = m2;
				m2 = i;

				/* check for stopping calculation depending on termination method */
				if (options->termination == TERM_PREC)
				{
					//changed the condition to global_maxResiduum
					if (global_maxResiduum < options->term_precision)
					{
						term_iteration = 0;
					}
				}
				else if (options->termination == TERM_ITER)
				{
					term_iteration--;
				}
				
			#ifdef _OPENMP
			}
			#endif

		}

	#ifdef _OPENMP
	}
	#endif

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

static
void
displayMatrixMPI (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, struct process_args* proc_args)
{
	int x, y;

	int const interlines = options->interlines;
	double* line = malloc(sizeof(double) * proc_args->working_columns);

	if (proc_args->rank == 0) 
	{
		uint64_t line_index;
		uint64_t line_offset;	
		printf("Matrix:\n");
		for (y = 0; y < 9; y++)
		{
			line[0] = 0;
			line[proc_args->working_columns + 1] = 0;
			// line to be printed.
			line_index = y * (interlines + 1);

			/**
			 * Frage alle Prozesse nach der Zeile. Der Prozess, der sie berechnet hat, antwortet mit der Zeile.
			*/
			if ( y < 8 && (line_index >= proc_args->starting_line) && (line_index < (proc_args->starting_line + proc_args->working_lines) ))
			{
				line = arguments->Matrix[results->m][line_index - proc_args->starting_line];
			}
			else if (y < 8)
			{
				// printf("Searching for line %d\n", line_index);
				MPI_Bcast(&line_index, 1, MPI_INT, 0, MPI_COMM_WORLD);
				MPI_Recv(line, proc_args->working_columns, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}

			/**
			 * Zeige die Matrix an. 0-Fälle, da die o.g. Matrizen ohne Randzeilen arbeiten.
			*/
			for (x = 0; x < 9; x++)
			{
				if(y == 0 || y == 8 || x == 0 || x == 8) 
				{
					printf ("%7.4f\t", 0.0);
				}
				else
				{
					printf ("%7.4f\t", line[x * (interlines + 1)]);
				}
			}
			printf ("\n");
		}

		/**
		 * Terminierungssignal als unmögliche Zeile. Könnte Probleme mit Uint64_t geben, ist aber innerhalb einer realisitischen Anzahl an Zeilen
		 * unwahrscheinlich.
		*/
		int terminate = -1;
		// printf("Sending termination signal\n");
		MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD);
		fflush (stdout);
	}
	else
	{
		/**
		 * Frage alle Prozesse nach der Zeile. Der Prozess, der sie berechnet hat, antwortet mit der Zeile.
		*/
		int start = (int) proc_args->starting_line;
		int end = (int) proc_args->starting_line + (int) proc_args->working_lines;

		int line = 0;
		int termination = -1;
		while (1)
		{
			MPI_Bcast(&line, 1, MPI_INT, 0, MPI_COMM_WORLD);

			
			if (line >= start
			&& line < end)
			{
				//printf("Rank %d sending line %d\n", proc_args->rank, line);
				MPI_Ssend(arguments->Matrix[results->m][line - start], proc_args->working_columns, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			else if (line == termination)
			{
				// printf("Rank %d terminating. Current line: %d, Range: %d - %d\n", proc_args->rank, line, proc_args->starting_line, proc_args->starting_line + proc_args->working_lines);
				break;
			}
		}
	}
		
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

	int rank, world_size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	if(options.method == METH_GAUSS_SEIDEL || world_size == 1) {
		if(rank == 0) {

			allocateMatrices(&arguments);
			initMatrices(&arguments, &options);

			gettimeofday(&start_time, NULL);
			calculate(&arguments, &results, &options);
			gettimeofday(&comp_time, NULL);

			displayStatistics(&arguments, &results, &options);
			displayMatrix(&arguments, &results, &options);

			freeMatrices(&arguments);

		}
	}else {

		if(rank < (int) arguments.N-1) {

			struct process_args proc_args;
		
			proc_args.world_size = world_size <= (int) arguments.N-1 ? (uint64_t) world_size : arguments.N-1;
			proc_args.rank = rank;

			allocateMatricesMPI(&arguments, &proc_args);
			initMatricesMPI(&arguments, &proc_args);

			if(rank == 0)
			{
				gettimeofday(&start_time, NULL);
			}

			calculateMPI(&arguments, &results, &options, &proc_args);

			if(rank == 0) 
			{
				gettimeofday(&comp_time, NULL);
				displayStatistics(&arguments, &results, &options);
			}

			displayMatrixMPI(&arguments, &results, &options, &proc_args);
			// printf("Rank %d back at main\n", rank);
			freeMatrices(&arguments);

		}

	}

	MPI_Finalize();

	return 0;
}
