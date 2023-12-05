#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


int*
init (int N)
{
	// TODO
	int* buf = (int*)malloc(sizeof(int) * N);

	srand(time(NULL));

	for (int i = 0; i < N; i++)
	{
		// Do not modify "% 25"
		buf[i] = rand() % 25;
	}

	return buf;
}

int
circle (int* buf, int array_size, int rank)
{
	int* tmp;

	tmp = (int*)malloc(sizeof(int) * array_size);

	tmp = buf;
	
	if (rank % 2)
	{
		MPI_Ssend(buf, array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Recv(buf, array_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else
	{
		MPI_Recv(buf, array_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Ssend(tmp, array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	

	return 0;
}

int
main (int argc, char** argv)
{
	int N;
	int rank;
	int* buf;
  	int ret;

	int world_size;
	int proc_array_size;

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

	proc_array_size = N / world_size;
	proc_array_size += (rank < N % world_size) ? 1 : 0;
	// TODO: Prozesse müssen wissen, wie groß ihr Array ist, bzw. Reste verteilen

	buf = init(proc_array_size);

	// TODO: zwei verschiedene verhalten, je nachdem ob man kleine oder Restgröße hat?

	if (rank == 0)
	{
		printf("\nBEFORE\n");
		// TODO: print whole array
		for (int i = 0; i < proc_array_size; i++)
		{
			printf("rank %d: %d\n", i, buf[i]);
		}
		printf("rank 0: %d\n", buf[]);
		for (int i = 1; i < world_size; i++)
		{
			// TODO: print whole array
			MPI_Recv(buf, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for
			printf("rank %d: %d\n", i, buf[i]);
		}

		
	}
	else
	{
		MPI_Send(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	ret = circle(buf);
	if (rank == 0)
	{
		printf("\nAFTER\n");

		for (int i = 1; i < world_size; i++)
		{
			MPI_Recv(buf, proc_array_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("rank %d: %d\n", i, buf[j]);
		}
	}
	else
	{
		MPI_Ssend(buf, proc_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	free(buf);
	MPI_Finalize();

	return EXIT_SUCCESS;
}
