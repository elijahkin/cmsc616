#include <cmath>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <vector>

/* Return the product of the vector x with every odd indexed element inverted.
   i.e. x_0 * 1/x_1 * x_2 * 1/x_3 * x_4 ...

   Example:
      input: [4, 2, 10, 4, 5]
      output: 25
*/
double productWithInverses(std::vector<double> const &x) {
  double result = 1.0;

#pragma omp parallel for reduction(* : result)
  for (int i = 0; i < x.size(); ++i) {
    result *= (i % 2) ? (1 / x[i]) : x[i];
  }

  return result;
}

int main(int argc, char **argv) {
  int N = 1024;
  int seed = 273;

  if (argc == 2) {
    N = std::stoi(argv[1]);
  }
  if (argc == 3) {
    N = std::stoi(argv[1]);
    seed = std::stoi(argv[2]);
  }

  std::vector<double> x(N);
  srand(seed);

  int points_max = N;
  int points_min = N - 1;

  for (int i = 0; i < N; i += 1) {
    x[i] = (rand() / (double)RAND_MAX) * (points_max - points_min) + points_min;
  }

  double totalTime = 0.0;
  double start = omp_get_wtime();

  double val = productWithInverses(x);
  printf("Product: %.5f\n", val);

  totalTime = omp_get_wtime() - start;
  printf("Time: %.5f\n", totalTime);
}
