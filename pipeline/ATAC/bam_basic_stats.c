#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "htslib/sam.h"
#include "htslib/khash.h"


#define bam_is_unmap(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define bam_is_dup(b) (((b)->core.flag & BAM_FDUP) != 0)



long total_read = 0;
long mapped_read = 0;
long unmapped_read = 0;
long unique_read = 0;
long duplicated_read = 0;

long nuclear_dup_read = 0 ;
long nuclear_unique_read = 0 ;
long non_nuclear_dup_read = 0;
long non_nuclear_unique_read = 0;

bool is_bam(const char *fname);
int strcasecmp(const char *s1, const char *s2);
char* rename_suffix(const char *fn);

int main(int argc, char **argv)
{
    if (argc == 1){
        fprintf(stderr, "Usage: %s input.bam  [stats.txt] [organelle genome size] \n;", argv[0]);
        fprintf(stderr, "e.g %s input.bam stats.txt 20000 \n;", argv[0]);
        exit(2);
    }

    char *bam_fn = argv[1];
    if (! is_bam(bam_fn)){
        fprintf(stderr, "%s not a BAM file !\n", argv[1]);
        exit(2);
    }

    char *stats_out_file = rename_suffix( bam_fn );
    long nuclear_size = 20000;
    if (argc > 3){
        stats_out_file = argv[2];
        nuclear_size = strtol(argv[3], NULL, 10);
    } else if (argc > 2){
        stats_out_file = argv[2];
    }


    samFile *fp_in = hts_open( bam_fn , "r");
    bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
    bam1_t *aln = bam_init1();

    uint16_t flag;
    int tid_len;
    while (  sam_read1( fp_in, bamHdr, aln) > 0 ){
        //if (aln->core.tid != aln->core.mtid) continue;
        total_read += 1;
        uint32_t cur_size =  bamHdr->target_len[aln->core.tid];
        //fprintf(stderr, "%s, %ld\n", bamHdr->target_name[aln->core.tid], cur_size);
        if ( bam_is_unmap(aln) ){
            unmapped_read += 1;
        } else{
            mapped_read += 1;
            if ( bam_is_dup(aln)){
                if (cur_size > nuclear_size){
                    nuclear_dup_read += 1;
                } else{
                    non_nuclear_dup_read += 1;
                }
                duplicated_read += 1;
            } else{
                if (cur_size > nuclear_size){
                    nuclear_unique_read += 1;
                } else{
                    non_nuclear_unique_read += 1;
                }
                unique_read += 1;
            }

        }

    }

    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in);

    FILE *fp_out;
    fp_out = fopen(stats_out_file, "w");
    fprintf(fp_out, "total_read\t%ld\n", total_read);
    fprintf(fp_out, "mapped_read\t%ld\n", mapped_read);
    fprintf(fp_out, "unmapped_read\t%ld\n", unmapped_read);
    fprintf(fp_out, "unique_read\t%ld\n", unique_read);
    fprintf(fp_out, "duplicated_read\t%ld\n", duplicated_read);
    fprintf(fp_out, "nuclear_unique_read\t%ld\n", nuclear_unique_read);
    fprintf(fp_out, "nuclear_dup_read\t%ld\n", nuclear_dup_read);
    fprintf(fp_out, "non_nuclear_unique_read\t%ld\n", non_nuclear_unique_read);
    fprintf(fp_out, "non_nuclear_dup_read\t%ld\n", non_nuclear_dup_read);
    fclose(fp_out);

}

// rename the suffix
char *rename_suffix(const char *fn){
    int l = strlen(fn);
    char *temp = malloc(  (size_t)  (l+13) * sizeof(char) );

    strcpy(temp, fn);
    strcpy(temp+l-4, "_basic_stats.txt");
    return temp;
}

// determine whether is bam or not
bool is_bam(const char *fname){
    int l = strlen(fname);
    if (l >= 4 && strcasecmp(fname+l-4, ".bam")!=0){
        return false;
    }
    return true;
}

// compare string is casse-ignoring model
int strcasecmp(const char *s1, const char *s2){
    int ret;
    
    int l = strlen(s1);
    char *temp = malloc( (size_t) l * sizeof(char));
    strcpy(temp, s1);

    int i ;
    for ( i = 0 ; temp[i] != '\0'; ++i){
        temp[i] = tolower(temp[i]);
    }
    ret  = strcmp(temp, s2);
    return ret;


}