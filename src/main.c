#include <assert.h>
#include <gsl/gsl_sort.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STATS_FILENAME "/fred/oz013/simulations/Tiamat/trees/forests_info-with_cumulative.hdf5"
#define FINAL_SNAP 60
#define N_RANKS 32

int main(int argc, char* argv[])
{
  // read in the unique forest_ids
  hid_t file_id = H5Fopen(STATS_FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);

  hsize_t n_ids = 0;
  {
    H5T_class_t class_id = 0;
    size_t type_size = 0;
    herr_t status = H5LTget_dataset_info(file_id, "info/forest_id", &n_ids, &class_id, &type_size);
    assert(status >= 0);
  }

  int* forest_ids = NULL;
  forest_ids = malloc(sizeof(int) * n_ids);
  assert(forest_ids != NULL);

  {
    herr_t status = H5LTread_dataset_int(file_id, "info/forest_id", forest_ids);
    assert(status >= 0);
  }

  // read in the final counts
  int* final_counts = NULL;
  final_counts = calloc(n_ids, sizeof(int));
  assert(final_counts != NULL);

  {
    char dset_name[128] = { '\0' };
    sprintf(dset_name, "snapshots/snap_%03d", FINAL_SNAP);

    herr_t status = H5LTread_dataset_int(file_id, dset_name, final_counts);
    assert(status >= 0);
  }

  // sort the final counts and store the sort indices (descending order)
  puts("Indirectly sorting forest_ids by final counts...");
  size_t* sort_ind = NULL;
  sort_ind = calloc(n_ids, sizeof(size_t));
  assert(sort_ind != NULL);
  gsl_sort_int_index(sort_ind, final_counts, 1, n_ids);
  {
    int ii = 0;
    int jj = 0;
    while (ii < jj) {
      int tmp = sort_ind[ii];
      sort_ind[ii] = sort_ind[jj];
      sort_ind[jj] = tmp;
      ii++;
      jj--;
    }
  }

  // {
  //   int* tmp = NULL;
  //   tmp = malloc(n_ids * sizeof(int));
  //   memcpy(tmp, final_counts, sizeof(int) * n_ids);
  //   assert(tmp != NULL);

  //   // TODO: reorder final_counts

  //   free(tmp);
  // }

  // use an array to store the number of halos assigned to each rank
  uint32_t rank_counts[N_RANKS] = { 0 };

  uint32_t rank_argsort_ind[N_RANKS];
  for (uint32_t ii = 0; ii < N_RANKS; ++ii) {
    rank_argsort_ind[ii] = ii;
  }

  int* snap_counts = NULL;
  snap_counts = calloc(n_ids, sizeof(int));
  assert(snap_counts != NULL);

  // Save old error handler and turn off error handling
  herr_t (*old_func)(long, void*) = NULL;
  void *old_client_data = NULL;
  hid_t estack_id = 0;
  H5Eget_auto(estack_id, &old_func, &old_client_data);
  H5Eset_auto(estack_id, NULL, NULL);

  puts("Calculating forest placements...");
  for (int snap = 0; snap < FINAL_SNAP; ++snap) {
    char dset_name[128] = { '\0' };
    sprintf(dset_name, "snapshots/snap_%03d", snap);

    herr_t status = H5LTread_dataset_int(file_id, dset_name, snap_counts);
    if (status < 0) {
      printf("%d ", snap);
      fflush(stdout);
      continue;
    }

    // loop through non-zero counts
    for (int ii = 0; ii < n_ids; ++ii) {
      int jj = sort_ind[ii];

      if ((snap_counts[jj] > 0) && (final_counts[jj] > 0)) {

        // first appearance!
        rank_counts[rank_argsort_ind[0]] += final_counts[jj];
        final_counts[jj] = 0;

        // keep the rank_argsort_inds correct
        for (int kk = 1; kk < N_RANKS; ++kk) {
          if (rank_counts[rank_argsort_ind[kk]] < rank_counts[rank_argsort_ind[kk - 1]]) {
            int tmp = rank_argsort_ind[kk];
            rank_argsort_ind[kk] = rank_argsort_ind[kk - 1];
            rank_argsort_ind[kk - 1] = tmp;
          }
        }
      }
    }
    printf("%d ", snap);
    fflush(stdout);
  }
  
  // Restore previous error handler
  H5Eset_auto(estack_id, old_func, old_client_data);

  // write the results
  FILE* fout = NULL;
  fout = fopen("results.txt", "w");
  assert(fout != NULL);
  for (int ii = 0; ii < N_RANKS - 1; ++ii) {
    fprintf(fout, "%d,", rank_counts[ii]);
  }
  fprintf(fout, "%d", rank_counts[N_RANKS - 1]);
  fclose(fout);

  free(sort_ind);
  free(snap_counts);
  free(final_counts);
  free(forest_ids);
  H5Fclose(file_id);
}
