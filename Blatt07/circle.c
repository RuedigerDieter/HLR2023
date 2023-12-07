#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <mpi.h>
#include <string.h>


int*
init (int N, int rank)
{
	int* buf = (int*)malloc(sizeof(int) * N);
	/*
	 * Rang und Faktor um die Streuung der Zufallszahlen zu erhöhen
	 */
	srand(time(NULL) + rank * 3892);

	for (int i = 0; i < N; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	return buf;
}

int
circle (int* buf, int array_size, int world_size, int rank)
{
	int* tmp_in;
	int* tmp_out;
	int* firstbuf;
	int b_iterate = 0;
	int iterations = 0;

	int src = (rank == 0)? world_size - 1 : rank - 1;
	int dest = (rank == world_size - 1)? 0 : rank + 1;

	/**
	 * Pufferarrays zur Kommunikation
	 */
	tmp_in = malloc(sizeof(int) * array_size);
	tmp_out = malloc(sizeof(int) * array_size);
	firstbuf = NULL;

	/**
	 * Wenn es nur einen Prozess gibt, matcht der erste automatisch den letzten.
	 *  Daher kann man hier Zeit sparen. 
	*/
	if (world_size == 1)
	{
		return 0;
	}

	if (rank == world_size - 1)
	{
		firstbuf = malloc(sizeof(int) * array_size);
	}


	/**
	 * Versenden vom letzten an den ersten Prozess für die Abbruchbedingung
	*/
	if (rank == 0)
	{
		MPI_Ssend(buf, array_size, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD);
	}
	else if (rank == world_size - 1)
	{
		MPI_Recv(firstbuf, array_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}


	while(!b_iterate)
	{
		iterations++;

		memcpy(tmp_out, buf, array_size * sizeof(int));

		if (rank == 0)
		{
			MPI_Ssend(tmp_out, array_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
			MPI_Recv(tmp_in, array_size, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Recv(tmp_in, array_size, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Ssend(tmp_out, array_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
		}

		memcpy (buf, tmp_in, array_size * sizeof(int));

		/**
		 * Überprüfe Abbruchbedingung. Broadcast aktualisiert die Abbruchbedingung für alle Prozesse.
		*/
		if (rank == world_size - 1)
		{
			b_iterate = (buf[0] == firstbuf[0])? 1 : 0;

			MPI_Bcast(&b_iterate, 1, MPI_INT, world_size - 1, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Bcast(&b_iterate, 1, MPI_INT, world_size - 1, MPI_COMM_WORLD);
		}
	}

	free(tmp_in);
	free(tmp_out);
	if (rank == world_size - 1)
	{
		free(firstbuf);
		printf("\n");
		printf("Iterations: %d\n", iterations);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}

int
main (int argc, char** argv)
{
	int N;
	int rank;
	int* buf;
  	int ret;
	int* buf_ret;

	int world_size;
	int proc_array_size;
	int rest;

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}

	// Array length
	N = atoi(argv[1]);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	rest = N % world_size;
	
	proc_array_size = N / world_size;
	proc_array_size += (rest != 0)? 1 : 0;

	buf_ret = (rank == 0)? (int*)malloc(sizeof(int) * proc_array_size) : NULL;
	buf = init(proc_array_size, rank);

	/**
	 *	Alle Prozesse haben die richtigen Elemente, die Prozesse mit "Lücken" kriegen illegale Werte,
	 *	um die Kommunikation sauberer zu gestalten
	 */
	if (rest != 0 && rank >= rest)
	{
		buf[proc_array_size - 1] = 25;
	}

	if (rank == 0)
	{
		printf("\nBEFORE\n");
		
		printf("rank 0: ");
		for (int i = 0; i < proc_array_size; i++)
		{
			if (!(buf[i] == 25))
				{
					printf("%d ", buf[i]);
				}		
		}
		printf("\n");

		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(buf_ret, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: ", i);
			for (int j = 0; j < proc_array_size; j++)
			{
				if (!(buf_ret[j] == 25))
				{
					printf("%d ", buf_ret[j]);
				}
			}
			printf("\n");
		}
	}
	else
	{
		MPI_Ssend(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	ret = circle(buf, proc_array_size, world_size, rank);

	if (rank == 0)
	{
		printf("\nAFTER\n");

		printf("rank 0: ");
		for (int i = 0; i < proc_array_size; i++)
		{
			if (!(buf[i] == 25))
			{
				printf("%d ", buf[i]);
			}
		}
		printf("\n");

		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(buf_ret, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: ", i);
			for (int j = 0; j < proc_array_size; j++)
			{
				if (!(buf_ret[j] == 25))
				{
					printf("%d ", buf_ret[j]);
				}			
			}
			printf("\n");
		}
	}
	else
	{
		MPI_Ssend(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	free(buf);
	MPI_Finalize();

	return ret;
}
