#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>


#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

typedef struct {
    signal* sig;
    int bandNum;
    int filter_order;
    double bandwidth;
    double* power;
} ThreadData;


void* worker(void* arg) {
    ThreadData* data = (ThreadData*) arg;
    
    generate_band_pass(data->sig->Fs,
                                               data->bandNum * data->bandwidth + 0.0001,
                                               (data->bandNum + 1) * data->bandwidth - 0.0001,
                                               data->filter_order, filter_coeffs);
    hamming_window(data->filter_order, filter_coeffs);

    convolve_and_compute_power(data->sig->num_samples,
                               data->sig->data,
                               data->filter_order,
                               filter_coeffs,
                               &(data->power[data->bandNum]));

    // free(filter_coeffs);
    pthread_exit(NULL);
}

void remove_dc(double* data, int num) {

    double dc = avg_of(data,num);
  
    printf("Removing DC component of %lf\n",dc);
  
    for (int i = 0; i < num; i++) {
      data[i] -= dc;
    }
  }


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {

    double Fc        = (sig->Fs) / 2;
    double bandwidth = Fc / num_bands;

    remove_dc(sig->data, sig->num_samples);

    double signal_power = avg_power(sig->data, sig->num_samples);


    int num_threads;
    int num_processors;

    double band_power[num_bands];
    // pthread_t* threadIDs[num_bands];  //add
    pthread_t* threadIDs;  //add
    ThreadData thread_data[num_bands]; //add
    threadIDs = (pthread_t*) malloc(sizeof(pthread_t) * num_threads);
    

    for (int i = 0; i < num_threads; i++) {
        ThreadData* t_data = malloc(sizeof(ThreadData));
        
        //set struct fields
        t_data->sig = sig;
        t_data->bandNum = i;
        t_data->power = band_power;

        int ret_code = pthread_create(&(threadIDs[i]), NULL, worker, (void*) t_data);
        if (ret_code != 0){
            perror("Failed thread start");
            exit(-1);
        }
    }
    
    //join all threads
    for (int i =0; i<num_threads; i++){
    int ret_code=pthread_join(threadIDs[i], NULL);

    if (ret_code !=0){
        perror("join failed");
        exit(-1);
    }
    }
}


// int main(int argc, char* argv[]) {
//     if (argc != 6) {
//         printf("Usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
//         return -1;
//     }

//     char sig_type = toupper(argv[1][0]);
//     char* sig_file = argv[2];
//     double Fs = atof(argv[3]);
//     int filter_order = atoi(argv[4]);
//     int num_bands = atoi(argv[5]);

//     signal* sig;
//     switch (sig_type) {
//         case 'T': sig = load_text_format_signal(sig_file); break;
//         case 'B': sig = load_binary_format_signal(sig_file); break;
//         case 'M': sig = map_binary_format_signal(sig_file); break;
//         default: printf("Unknown signal type\n"); return -1;
//     }

//     if (!sig) {
//         printf("Unable to load or map file\n");
//         return -1;
//     }
//     sig->Fs = Fs;

//     double start = 0, end = 0;
//     if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
//         printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
//     } else {
//         printf("no aliens\n");
//     }

//     free_signal(sig);
//     return 0;
// }