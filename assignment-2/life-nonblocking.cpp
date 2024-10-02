#include "mpi.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

/*
 * Reads the input file line by line and stores it in a 2D matrix.
 */
void read_input_file(int **life, string const &input_file_name) {

  // Open the input file for reading.
  ifstream input_file;
  input_file.open(input_file_name);
  if (!input_file.is_open())
    perror("Input file cannot be opened");

  string line, val;
  int x, y;
  while (getline(input_file, line)) {
    stringstream ss(line);

    // Read x coordinate.
    getline(ss, val, ',');
    x = stoi(val);

    // Read y coordinate.
    getline(ss, val);
    y = stoi(val);

    // Populate the life matrix.
    life[x][y] = 1;
  }
  input_file.close();
}

/*
 * Writes out the final state of the 2D matrix to a csv file.
 */
void write_output(int **result_matrix, int X_limit, int Y_limit,
                  string const &input_name, int num_of_generations) {

  // Open the output file for writing.
  ofstream output_file;
  string input_file_name = input_name.substr(0, input_name.length() - 5);
  output_file.open(input_file_name + "." + to_string(num_of_generations) +
                   ".parallel.csv");
  if (!output_file.is_open())
    perror("Output file cannot be opened");

  // Output each live cell on a new line.
  for (int i = 0; i < X_limit; i++) {
    for (int j = 0; j < Y_limit; j++) {
      if (result_matrix[i][j] == 1) {
        output_file << i << "," << j << "\n";
      }
    }
  }
  output_file.close();
}

/*
 * Processes the life array for the specified number of iterations.
 */
void compute(int **life, int **previous_life, int X_limit, int Y_limit) {
  int neighbors = 0;

  int myrank, numpes;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numpes);

  MPI_Request req_prev, req_next;

  int prev = (myrank == 0) ? MPI_PROC_NULL : myrank - 1;
  int next = (myrank == numpes - 1) ? MPI_PROC_NULL : myrank + 1;

  MPI_Isend(&life[0], Y_limit, MPI_INT, prev, 0, MPI_COMM_WORLD, &req_prev);
  MPI_Isend(&life[X_limit - 1], Y_limit, MPI_INT, next, 0, MPI_COMM_WORLD,
            &req_next);
  MPI_Irecv(&previous_life[0], Y_limit, MPI_INT, prev, 0, MPI_COMM_WORLD,
            &req_next);
  MPI_Irecv(&previous_life[X_limit + 1], Y_limit, MPI_INT, next, 0,
            MPI_COMM_WORLD, &req_prev);

  // Update the previous_life matrix with the current life matrix state.
  for (int i = 0; i < X_limit; i++) {
    for (int j = 0; j < Y_limit; j++) {
      previous_life[i + 1][j + 1] = life[i][j];
    }
  }

  MPI_Status stat_prev, stat_next;
  MPI_Wait(&req_prev, &stat_prev);
  MPI_Wait(&req_next, &stat_next);

  // For simulating each generation, calculate the number of live
  // neighbors for each cell and then determine the state of the cell in
  // the next iteration.
  for (int i = 1; i < X_limit + 1; i++) {
    for (int j = 1; j < Y_limit + 1; j++) {
      neighbors = previous_life[i - 1][j - 1] + previous_life[i - 1][j] +
                  previous_life[i - 1][j + 1] + previous_life[i][j - 1] +
                  previous_life[i][j + 1] + previous_life[i + 1][j - 1] +
                  previous_life[i + 1][j] + previous_life[i + 1][j + 1];

      if (previous_life[i][j] == 0) {
        // A cell is born only when an unoccupied cell has 3 neighbors.
        if (neighbors == 3)
          life[i - 1][j - 1] = 1;
      } else {
        // An occupied cell survives only if it has either 2 or 3 neighbors.
        // The cell dies out of loneliness if its neighbor count is 0 or 1.
        // The cell also dies of overpopulation if its neighbor count is 4-8.
        if (neighbors != 2 && neighbors != 3) {
          life[i - 1][j - 1] = 0;
        }
      }
    }
  }
}

/**
 * The main function to execute "Game of Life" simulations on a 2D board.
 */
int main(int argc, char *argv[]) {

  if (argc != 5)
    perror("Expected arguments: ./life <input_file> <num_of_generations> "
           "<X_limit> <Y_limit>");

  float min_time, sum_time, max_time;

  int myrank, numpes;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numpes);

  string input_file_name = argv[1];
  int num_of_generations = stoi(argv[2]);
  int X_limit = stoi(argv[3]);
  int Y_limit = stoi(argv[4]);

  int **global_life = nullptr;
  if (myrank == 0) {
    global_life = new int *[X_limit];
    for (int i = 0; i < X_limit; i++) {
      global_life[i] = new int[Y_limit];
      for (int j = 0; j < Y_limit; j++) {
        global_life[i][j] = 0;
      }
    }

    read_input_file(global_life, input_file_name);
  }

  int X_limit_proc = X_limit / numpes;

  int **life = new int *[X_limit_proc];
  for (int i = 0; i < X_limit_proc; i++) {
    life[i] = new int[Y_limit];
    // TODO don't need to zero these?
    for (int j = 0; j < Y_limit; j++) {
      life[i][j] = 0;
    }
  }

  MPI_Scatter(&global_life, X_limit_proc, MPI_INT, &life, X_limit_proc, MPI_INT,
              0, MPI_COMM_WORLD);

  // Use previous_life to track the previous state of the board.
  // Pad the previous_life matrix with 0s on all four sides by setting all
  // cells in the following rows and columns to 0:
  //  1. Row 0
  //  2. Column 0
  //  3. Row X_limit+1
  //  4. Column Y_limit+1
  int **previous_life = new int *[X_limit_proc + 2];
  for (int i = 0; i < X_limit_proc + 2; i++) {
    previous_life[i] = new int[Y_limit + 2];
    for (int j = 0; j < Y_limit + 2; j++) {
      previous_life[i][j] = 0;
    }
  }

  clock_t start = clock();
  for (int numg = 0; numg < num_of_generations; numg++) {
    compute(life, previous_life, X_limit, Y_limit);
  }
  clock_t end = clock();
  float local_time = float(end - start) / CLOCKS_PER_SEC;

  MPI_Reduce(&local_time, &min_time, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_time, &sum_time, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_time, &max_time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Gather(&life, X_limit_proc, MPI_INT, &global_life, X_limit_proc, MPI_INT,
             0, MPI_COMM_WORLD);

  if (myrank == 0) {
    // For serial code: min, avg, max are the same
    cout << "TIME: Min: " << min_time << " s Avg: " << sum_time / numpes
         << " s Max: " << max_time << " s\n";

    // Write out the final state to the output file.
    write_output(global_life, X_limit, Y_limit, input_file_name,
                 num_of_generations);
  }

  for (int i = 0; i < X_limit; i++) {
    delete life[i];
  }
  for (int i = 0; i < X_limit + 2; i++) {
    delete previous_life[i];
  }
  delete[] life;
  delete[] previous_life;

  MPI_Finalize();
  return 0;
}
