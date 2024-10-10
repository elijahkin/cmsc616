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
void read_input_file(int *life, string const &input_file_name, int Y_limit) {

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
    life[x * Y_limit + y] = 1;
  }
  input_file.close();
}

/*
 * Writes out the final state of the 2D matrix to a csv file.
 */
void write_output(int *result_matrix, int X_limit, int Y_limit,
                  string const &input_name, int num_of_generations) {

  // Open the output file for writing.
  ofstream output_file;
  string input_file_name = input_name.substr(0, input_name.length() - 5);
  output_file.open(input_file_name + "." + to_string(num_of_generations) +
                   ".csv");
  if (!output_file.is_open())
    perror("Output file cannot be opened");

  // Output each live cell on a new line.
  for (int i = 0; i < X_limit; ++i) {
    for (int j = 0; j < Y_limit; ++j) {
      if (result_matrix[i * Y_limit + j] == 1) {
        output_file << i << "," << j << "\n";
      }
    }
  }
  output_file.close();
}

/*
 * The main function to execute "Game of Life" simulations on a 2D board.
 */
int main(int argc, char *argv[]) {

  if (argc != 5)
    perror("Expected arguments: ./life <input_file> <num_of_generations> "
           "<X_limit> <Y_limit>");

  double min_time, sum_time, max_time;

  int myrank, numpes;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numpes);

  string input_file_name = argv[1];
  int num_of_generations = stoi(argv[2]);
  int X_limit = stoi(argv[3]);
  int Y_limit = stoi(argv[4]);

  int *global_life = nullptr;
  if (myrank == 0) {
    global_life = new int[X_limit * Y_limit];
    read_input_file(global_life, input_file_name, Y_limit);
  }

  const int X_limit_proc = X_limit / numpes;
  const int prev = (myrank == 0) ? MPI_PROC_NULL : myrank - 1;
  const int next = (myrank == numpes - 1) ? MPI_PROC_NULL : myrank + 1;

  int *life = new int[X_limit_proc * Y_limit];
  MPI_Scatter(global_life, X_limit_proc * Y_limit, MPI_INT, life,
              X_limit_proc * Y_limit, MPI_INT, 0, MPI_COMM_WORLD);

  // Use previous_life to track the previous state of the board.
  // Pad the previous_life matrix with 0s on all four sides by setting all
  // cells in the following rows and columns to 0:
  //  1. Row 0
  //  2. Column 0
  //  3. Row X_limit+1
  //  4. Column Y_limit+1
  int *previous_life = new int[(X_limit_proc + 2) * (Y_limit + 2)];
  for (int i = 0; i < (X_limit_proc + 2) * (Y_limit + 2); ++i) {
    previous_life[i] = 0;
  }

  MPI_Request requests[4];

  double start = MPI_Wtime();
  for (int numg = 0; numg < num_of_generations; ++numg) {
    MPI_Isend(&life[0], Y_limit, MPI_INT, prev, 0, MPI_COMM_WORLD,
              &requests[0]);
    MPI_Isend(&life[(X_limit_proc - 1) * Y_limit], Y_limit, MPI_INT, next, 0,
              MPI_COMM_WORLD, &requests[1]);

    MPI_Irecv(&previous_life[1], Y_limit, MPI_INT, prev, 0, MPI_COMM_WORLD,
              &requests[2]);
    MPI_Irecv(&previous_life[(X_limit_proc + 1) * (Y_limit + 2) + 1], Y_limit,
              MPI_INT, next, 0, MPI_COMM_WORLD, &requests[3]);

    // Update the previous_life matrix with the current life matrix state.
    for (int i = 0; i < X_limit_proc; ++i) {
      for (int j = 0; j < Y_limit; ++j) {
        previous_life[(i + 1) * (Y_limit + 2) + (j + 1)] =
            life[i * Y_limit + j];
      }
    }

    // For simulating each generation, calculate the number of live
    // neighbors for each cell and then determine the state of the cell in
    // the next iteration.
    int neighbors;
    for (int i = 2; i < X_limit_proc; ++i) {
      const int prev_row = (i - 1) * (Y_limit + 2);
      const int this_row = i * (Y_limit + 2);
      const int next_row = (i + 1) * (Y_limit + 2);
      for (int j = 1; j < Y_limit + 1; ++j) {
        neighbors =
            previous_life[prev_row + (j - 1)] + previous_life[prev_row + j] +
            previous_life[prev_row + (j + 1)] +
            previous_life[this_row + (j - 1)] +
            previous_life[this_row + (j + 1)] +
            previous_life[next_row + (j - 1)] + previous_life[next_row + j] +
            previous_life[next_row + (j + 1)];

        // A cell is born only when an unoccupied cell has 3 neighbors.
        // An occupied cell survives only if it has either 2 or 3 neighbors.
        // The cell dies out of loneliness if its neighbor count is 0 or 1.
        // The cell also dies of overpopulation if its neighbor count is 4-8.
        life[(i - 1) * Y_limit + (j - 1)] =
            (neighbors == 3) ||
            (previous_life[i * (Y_limit + 2) + j] && (neighbors == 2));
      }
    }

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

    for (int i = 1; i < X_limit_proc + 1; i += X_limit_proc - 1) {
      for (int j = 1; j < Y_limit + 1; ++j) {
        neighbors = previous_life[(i - 1) * (Y_limit + 2) + (j - 1)] +
                    previous_life[(i - 1) * (Y_limit + 2) + j] +
                    previous_life[(i - 1) * (Y_limit + 2) + (j + 1)] +
                    previous_life[i * (Y_limit + 2) + (j - 1)] +
                    previous_life[i * (Y_limit + 2) + (j + 1)] +
                    previous_life[(i + 1) * (Y_limit + 2) + (j - 1)] +
                    previous_life[(i + 1) * (Y_limit + 2) + j] +
                    previous_life[(i + 1) * (Y_limit + 2) + (j + 1)];

        // A cell is born only when an unoccupied cell has 3 neighbors.
        // An occupied cell survives only if it has either 2 or 3 neighbors.
        // The cell dies out of loneliness if its neighbor count is 0 or 1.
        // The cell also dies of overpopulation if its neighbor count is 4-8.
        life[(i - 1) * Y_limit + (j - 1)] =
            (neighbors == 3) ||
            (previous_life[i * (Y_limit + 2) + j] && (neighbors == 2));
      }
    }
  }
  double end = MPI_Wtime();
  double local_time = end - start;

  MPI_Reduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_time, &sum_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Gather(life, X_limit_proc * Y_limit, MPI_INT, global_life,
             X_limit_proc * Y_limit, MPI_INT, 0, MPI_COMM_WORLD);

  if (myrank == 0) {
    // For serial code: min, avg, max are the same
    cout << "TIME: Min: " << min_time << " s Avg: " << sum_time / numpes
         << " s Max: " << max_time << " s\n";

    // Write out the final state to the output file.
    write_output(global_life, X_limit, Y_limit, input_file_name,
                 num_of_generations);
  }

  delete[] life;
  delete[] previous_life;
  if (myrank == 0) {
    delete[] global_life;
  }

  MPI_Finalize();
  return 0;
}
