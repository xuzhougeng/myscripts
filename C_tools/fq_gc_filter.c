#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "kseq.h"
#include <stdbool.h>

KSEQ_INIT(gzFile, gzread)

double calc_gc(char *seq, size_t seq_size);

int main(int argc, char *argv[])
{
    gzFile fp; 
    kseq_t *seq;
    int l ; //flag to indicate  seq length

    if (argc == 1){
        fprintf(stderr, "Usage: %s <in.fastq> min_gc max_gc\n", argv[0]);
        fprintf(stderr, "default min_gc: 35 \t max_gc: 55\n");
        return 1;
    }
    unsigned min_gc, max_gc;
    
    if (argc == 2){
        min_gc = 35;
        max_gc = 55;

    } else if (argc == 3)
    {
        min_gc = strtoul(argv[2], NULL, 10);
        max_gc = 55;
    } else{
        min_gc = strtoul(argv[2], NULL, 10);
        max_gc = strtoul(argv[3], NULL, 10);
    }

    fp = gzopen(argv[1], "r");
    seq = kseq_init(fp);
    double gc;
    while ( (l = kseq_read(seq)) >=0 ){
        gc = calc_gc(seq->seq.s, seq->seq.l);
        //fprintf(stderr, "%s, %f\n", seq->name.s, gc);
        if (gc > min_gc && gc < max_gc){
            fprintf(stdout, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s );
        }
    }

    return 0;

}

double calc_gc(char *seq, size_t seq_size){
    double gc_ratio;
    int gc_number;
    int table[5] = {0,0,0,0,0}; // A,T,C,G
    char c;
    int i;
    for (i = 0; seq[i] != '\0'; ++i){
        c = toupper(seq[i]);
        switch (c)
        {
        case 'A':
            ++table[0];
            break;
        case 'T':
            ++table[1];
            break;
        case 'C':
            ++table[2];
            break;
        case 'G':
            ++table[3];
            break;
        default:
            ++table[4];
            break;
        }
    }
    gc_number = table[2] + table[3];
        
    gc_ratio = (double)gc_number / (double)seq_size  * 100; 
    //fprintf(stderr, "%d, %d, %f\n", gc_number, seq_size, gc_ratio);
        
    return gc_ratio;
}