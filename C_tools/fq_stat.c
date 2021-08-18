#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kseq.h"
#include "stdbool.h"


KSEQ_INIT(gzFile, gzread)

typedef struct{
    unsigned long long Q;
    unsigned long long Q20;
    unsigned long long Q30;
} q_val_t;

typedef struct{
    unsigned long long base_number;
    unsigned long long read_number;
    double error_rate;
    double Q20;
    double Q30;
    double GC_ratio;
} stat_t; 

q_val_t calc_q_value(kstring_t *qual, int offset);
int count_gc(kstring_t *seq);
stat_t fastq_stat(char **files, int file_number);


int main(int argc, char **argv)
{

    stat_t sample_stat = fastq_stat(argv+1, argc-1);
    fprintf(stdout, "read_number\tbase_number\tError\tQ20\tQ30\tGC_ratio\n");
    fprintf(stdout, "%ld\t%ld\t%.4f\t%.2f\t%.2f\t%.2f\n", 
                  sample_stat.read_number,
                  sample_stat.base_number,
                  sample_stat.error_rate * 100,
                  sample_stat.Q20 * 100,
                  sample_stat.Q30 * 100,
                  sample_stat.GC_ratio * 100);


    return 0;
}


stat_t fastq_stat(char **files, int file_number)
{

    unsigned long long base_number = 0;
    unsigned long long read_number = 0;
    unsigned long long Q_val = 0 ;
    unsigned long long Q20_number = 0;
    unsigned long long Q30_number = 0;
    unsigned long long GC_base_number = 0;
    gzFile fp;
    kseq_t *seq;
    int l;
    int i;
    for (i = 0; i < file_number; ++i){
        char *file = files[i];            
        fp = gzopen(file, "r");
        seq = kseq_init(fp);
        q_val_t q_val = {0,0,0};
        while ( (l = kseq_read(seq)) >= 0){
            base_number += l;
            read_number += 1;
            q_val = calc_q_value(&seq->qual, 33);
            Q_val += q_val.Q;
            Q20_number += q_val.Q20;
            Q30_number += q_val.Q30;
            GC_base_number += count_gc(&seq->seq);

        }
        kseq_destroy(seq);
        gzclose(fp);
    
    }

    double mean_Q = (double)Q_val / base_number;
    double error_rate = pow(10.0, mean_Q /(-10));

    double Q20 = (double)Q20_number / base_number;
    double Q30 = (double)Q30_number / base_number;
    double GC_ratio = (double)GC_base_number / base_number;

    stat_t stat = {base_number, read_number, error_rate, Q20, Q30, GC_ratio};

    return stat;
}

/*
Starting in Illumina 1.8, the quality scores 
have basically returned to the use of the Sanger format (Phred+33).
*/
q_val_t calc_q_value(kstring_t *qual, int offset)
{
    size_t i;
    char *phred = qual->s;
    int c;
    q_val_t q_val = {0,0,0};
    int q = 0;
    for (i = 0; i < qual->l; ++i){
        c = phred[i];
        q = c - offset;
        q_val.Q += q;
        if (q >= 20 ){
            ++(q_val.Q20);
        }
        if ( q >= 30 ){
            ++(q_val.Q30);
        }
    }
    
    return q_val;

}

int count_gc(kstring_t *seq)
{
    int gc_number = 0;
    int i = 0;
    char c;
    for (i = 0; i < seq->l; ++i){
        c = toupper(seq->s[i]);
        if ( c == 'C' || c == 'G' ){
            ++gc_number;
        }
    }
    return gc_number;
}