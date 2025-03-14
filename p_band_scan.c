
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#include <sched.h>   // for processor affinity
#include <unistd.h>  // unix standard apis
#include <pthread.h> // pthread api

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0


int num_threads;
int num_processors;
int num_bands;
int filter_order;
int bandwidth;

//array for each thread
pthread_t* threadIDs;

// //array holding params for each thread
// ThreadData* thread_data;

//params for each thread 
typedef struct {
    signal* sig;
    int bandNum;
    // int filter_order;
    // double bandwidth;
    double*  power; //pointer to a part in the final band_power array
} ThreadData;

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

void* worker(void* arg) {
    ThreadData* data = (ThreadData*) arg;
    int bandNum = data->bandNum;
    int blocksize = num_bands / num_threads;
    // double* filter_coeffs = (double*) malloc((data->filter_order+1)* sizeof(double));
    double filter_coeffs[filter_order+1];

    //set processor 
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(data->bandNum % num_processors, &set);
    if (sched_setaffinity(0, sizeof(set), &set) < 0)
    {
        perror("Can't set affinity");
        free(filter_coeffs);
        pthread_exit(NULL);
    }


    //figuring out chunk of band to work on based on band number
    int mystart = bandNum * blocksize;
    int myend = 0;


    if (bandNum == (num_threads - 1)){
        //additional space for end if leftover
        myend = mystart + blocksize + (num_bands % num_threads);
    } else {
        myend = mystart + blocksize;
    }


    //getting power values for each band block

    for (int i = mystart; i < myend; i++){
        generate_band_pass(data->sig->Fs, i * bandwidth + 0.0001,
                           (i + 1) * bandwidth - 0.0001,
                           filter_order, filter_coeffs);
        hamming_window(filter_order, filter_coeffs);

        convolve_and_compute_power(data->sig->num_samples,
                                   data->sig->data,
                                   filter_order,
                                   filter_coeffs,
                                   &(data->power[i]));
    }

    pthread_exit(NULL);
}



int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {

    double Fc        = (sig->Fs) / 2;
    bandwidth = Fc / num_bands;

    remove_dc(sig->data, sig->num_samples);

    double signal_power = avg_power(sig->data, sig->num_samples);

    printf("signal average power:     %lf\n", signal_power);

    //start timing
    resources rstart;
    get_resources(&rstart, THIS_PROCESS);
    double start = get_seconds();
    unsigned long long tstart = get_cycle_count();

    // double band_power[num_bands];
    // pthread_t* threadIDs[num_bands];  //add
    
    //initialize to inputed num of bands
    //holds each thread params for each thread & computed power
    // thread_data = (ThreadData*) malloc(sizeof(ThreadData) * num_bands);
    
    //pthread ids
    threadIDs = (pthread_t*) malloc(sizeof(pthread_t) * num_threads);
    double band_power[num_bands];

    // parallelize for each band
    for (int i = 0; i < num_threads; i++) {
        ThreadData* t_data = malloc(sizeof(ThreadData));
        
        //set struct fields
        t_data->sig = sig;
        t_data->bandNum = i;
        t_data->power = band_power;
        // t_data->bandwidth = bandwidth;
        // t_data->filter_order = filter_order;

        // t_data->power = band_power;

        int ret_code = pthread_create(&(threadIDs[i]), NULL, worker, (void*) t_data);
        if (ret_code != 0){
            perror("Failed thread start");
            exit(-1);
        }
    }
    
    //join all threads-- def correct
    for (int i =0; i<num_threads; i++){

        int ret_code = pthread_join(threadIDs[i], NULL);
        if (ret_code !=0){
            perror("join failed");
            exit(-1);
        }
    }

    //all results are stored in thread_data structs


    //stop timing
    unsigned long long tend = get_cycle_count();
    double end = get_seconds();

    resources rend;
    get_resources(&rend, THIS_PROCESS);

    resources rdiff;
    get_resources_diff(&rstart, &rend, &rdiff);

    //analyze band reuslts

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

    printf("Resource usages:\n\
        User time        %lf seconds\n\
        System time      %lf seconds\n\
        Page faults      %ld\n\
        Page swaps       %ld\n\
        Blocks of I/O    %ld\n\
        Signals caught   %ld\n\
        Context switches %ld\n",
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






int main(int argc, char* argv[]) {
    if (argc != 8) {
        printf("usage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
        return -1;
    }

    char sig_type = toupper(argv[1][0]);
    char* sig_file = argv[2];
    double Fs = atof(argv[3]);
    filter_order = atoi(argv[4]);
    num_bands = atoi(argv[5]);
    num_threads = atoi(argv[6]);
    num_processors = atoi(argv[7]);

    assert(Fs > 0.0);
    assert(filter_order > 0 && !(filter_order & 0x1));
    assert(num_bands > 0);

    signal* sig;
    switch (sig_type) {
        case 'T': sig = load_text_format_signal(sig_file); break;
        case 'B': sig = load_binary_format_signal(sig_file); break;
        case 'M': sig = map_binary_format_signal(sig_file); break;
        default: printf("Unknown signal type\n"); return -1;
    }

    if (!sig) {
        printf("Unable to load or map file\n");
        return -1;
    }
    sig->Fs = Fs;

    double start = 0, end = 0;
    if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
        printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
    } else {
        printf("no aliens\n");
    }

    free_signal(sig);
    return 0;
}