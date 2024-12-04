#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api


#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW 50000.0
#define ALIENS_HIGH 150000.0

void usage()
{
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

double avg_power(double *data, int num)
{

  double ss = 0;
  for (int i = 0; i < num; i++)
  {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double *data, int num)
{

  double m = data[0];
  for (int i = 1; i < num; i++)
  {
    if (data[i] > m)
    {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double *data, int num)
{

  double s = 0;
  for (int i = 0; i < num; i++)
  {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double *data, int num)
{

  double dc = avg_of(data, num);

  printf("Removing DC component of %lf\n", dc);

  for (int i = 0; i < num; i++)
  {
    data[i] -= dc;
  }
}




typedef struct Inputs {
  signal * sig;
  int filter_order;
  int num_bands;
  int num_processors;
  int myid;
  //int band;
  double* band_power;
  double* filter_coeffs;
  //int ** thread_bands;
  int num_threads;
}inputs;

void* worker(void* arg) {
  
  inputs* worker_input = (inputs*)arg;
  /*
  signal* sig = worker_input->sig;
  int filter_order = worker_input->filter_order;
  int num_bands = worker_input->num_bands;
  int num_processors = worker_input->num_processors;
  int myid = worker_input->myid;
  int band = worker_input->band;
  double* band_power = worker_input->band_power;
  double* filter_coeffs = worker_input->filter_coeffs;
  */
  

  double bandwidth = (worker_input->sig->Fs) / 2 / worker_input->num_bands;

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(worker_input->myid % worker_input->num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }

  //int arr_len = (worker_input->num_bands / worker_input->num_threads + ((worker_input->num_bands % worker_input->num_threads) ? 1 : 0));


  
  for (int i = 0; i < worker_input->num_bands; i++){
    if (worker_input->myid == i % worker_input->num_threads) {
      //worker_input->band = worker_input->thread_bands[worker_input->myid][i];
      worker_input->filter_coeffs = malloc(sizeof(double) * (worker_input->filter_order + 1));
      // Make the filter
      generate_band_pass(worker_input->sig->Fs,
                          i * bandwidth + 0.0001, // keep within limits
                          (i + 1) * bandwidth - 0.0001,
                          worker_input->filter_order,
                          worker_input->filter_coeffs);
      hamming_window(worker_input->filter_order, worker_input->filter_coeffs);

      // Convolve
      convolve_and_compute_power(worker_input->sig->num_samples,
                                  worker_input->sig->data,
                                  worker_input->filter_order,
                                  worker_input->filter_coeffs,
                                  &(worker_input->band_power[i]));
    }
  }
  pthread_exit(NULL);
}

// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //

int analyze_signal(signal *sig, int filter_order, int num_bands, double *lb, double *ub, int num_processors, int num_threads)
{
  pthread_t* tid = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);

  double Fc = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data, sig->num_samples);

  double signal_power = avg_power(sig->data, sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart, THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  //double* filter_coeffs = malloc(sizeof(double) * (filter_order + 1));
  double* band_power = malloc(sizeof(double) * num_bands);

  /*
  int** thread_bands = malloc(sizeof(int*) * num_threads);


  for (int i = 0; i < num_threads; i++){
    thread_bands[i] = malloc(sizeof(int) * (num_bands / num_threads + ((num_bands % num_threads) ? 1 : 0)));
    for (int j = 0; j < (num_bands / num_threads + ((num_bands % num_threads) ? 1 : 0)); j++){
      thread_bands[i][j] = i + j * num_threads;
    }
  }
  */

  for (int thread = 0; thread < num_threads; thread++)
  {
    inputs* worker_input = malloc(sizeof(inputs));
    // set inputs
    worker_input->sig = sig;
    worker_input->filter_order = filter_order;
    worker_input->num_bands = num_bands;
    worker_input->num_processors = num_processors;
    worker_input->myid = thread;
    //worker_input->band = band;
    worker_input->band_power = band_power;
    //worker_input->thread_bands = thread_bands;
    worker_input->num_threads = num_threads;


    if (pthread_create(&(tid[worker_input->myid]), NULL, worker, (void*)worker_input)) {
      perror("Failed to start thread");
      exit(-1);
    }

    /*
    // Make the filter
    generate_band_pass(sig->Fs,
                       band * bandwidth + 0.0001, // keep within limits
                       (band + 1) * bandwidth - 0.0001,
                       filter_order,
                       filter_coeffs);
    hamming_window(filter_order, filter_coeffs);

    // Convolve
    convolve_and_compute_power(sig->num_samples,
                               sig->data,
                               filter_order,
                               filter_coeffs,
                               &(band_power[band]));
    */
  }

  // now we will join all the threads
  for (int i = 0; i < num_threads; i++) {
    if (pthread_join(tid[i], NULL)) {
      perror("join failed");
      exit(-1);
    }
  }

  unsigned long long tend = get_cycle_count();
  double end = get_seconds();

  resources rend;
  get_resources(&rend, THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power, num_bands);
  double avg_band_power = avg_of(band_power, num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++)
  {
    double band_low = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++)
    {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH))
    {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power)
      {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0)
        {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      }
      else
      {
        printf("(meh)");
      }
    }
    else
    {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n"
         "User time        %lf seconds\n"
         "System time      %lf seconds\n"
         "Page faults      %ld\n"
         "Page swaps       %ld\n"
         "Blocks of I/O    %ld\n"
         "Signals caught   %ld\n"
         "Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //
// **************************************************** //






int main(int argc, char *argv[])
{

  if (argc != 8)
  {
    usage();
    return -1;
  }

  char sig_type = toupper(argv[1][0]);
  char *sig_file = argv[2];
  double Fs = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands = atoi(argv[5]);
  
  // sets num_threads and num_processors equal to the inputs after the normal band_scan.c inputs

  int num_threads = atoi(argv[6]); // number of threads
  int num_processors = atoi(argv[7]); // numer of processors to use

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n"
         "file:     %s\n"
         "Fs:       %lf Hz\n"
         "order:    %d\n"
         "bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal *sig;
  switch (sig_type)
  {
  case 'T':
    sig = load_text_format_signal(sig_file);
    break;

  case 'B':
    sig = load_binary_format_signal(sig_file);
    break;

  case 'M':
    sig = map_binary_format_signal(sig_file);
    break;

  default:
    printf("Unknown signal type\n");
    return -1;
  }

  if (!sig)
  {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_processors, num_threads))
  {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  }
  else
  {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
