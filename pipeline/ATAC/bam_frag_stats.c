#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "htslib/sam.h"
#include "htslib/khash.h"

#define ARRSIZE 1001
long frag_size[ARRSIZE] = {0};

//hash set for qname de-duplication
KHASH_SET_INIT_STR(str)

bool is_bam(char *fname);
int strcasecmp(const char *s1, const char *s2);
char* rename_suffix(const char *fn);

int main(int argc, char *argv[]){
  //hts_verbose = 0; // suppressed htslib warnnings	

  if ( argc == 1){
    fprintf(stderr, "Usage: %s input.bam  [frag_size.txt] \n", argv[0]);
    exit(2);
  }
  char *frag_size_fname;
  if (argc == 2){
    frag_size_fname = rename_suffix(argv[1]);
  } else{
    frag_size_fname = argv[2];
  }
  
 

  // init hash set
  khash_t(str) *qname_set;
  qname_set = kh_init(str);
  khint_t k;
  int absent;


  //load BAM file
  char *bam_fn = argv[1];

  if (! is_bam(bam_fn) ){
        fprintf(stderr, "%s not a BAM file!\n", argv[1]);
        exit(2);
  } 
  fprintf(stderr, "Process %s !\n", argv[1]);

  samFile *fp_in = hts_open( bam_fn, "r" ); // open bam file
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in); // read header
  bam1_t *aln = bam_init1(); // initialize an alignment
  

  //iterate the bam file
  while ( sam_read1( fp_in, bamHdr, aln) > 0){
    // filter the alignment which not on the same tig
    if (aln->core.tid != aln->core.mtid) continue;
    char *qname = bam_get_qname(aln);
    
    // check qname whether in the set or not
    // if qname in the set, remove the qname, and the target-frag_size array +1
    k = kh_get(str, qname_set, qname);
    if (k == kh_end(qname_set)){
      k = kh_put(str, qname_set, qname, &absent);
      continue;
    } 
    kh_del(str,qname_set, k);

    long insize = abs( aln->core.isize );
    if (insize > ARRSIZE - 1){
      ++frag_size[ARRSIZE-1];
    } else{
      ++frag_size[insize];
    }
  }

  bam_destroy1(aln);
  bam_hdr_destroy(bamHdr);
  sam_close(fp_in);

  // write the array
  FILE *fp_out;
  fp_out = fopen(frag_size_fname, "w");
  long i ;
  for (i = 0; i < ARRSIZE; ++i ){
    fprintf(fp_out, "%ld\t%ld\n", i, frag_size[i] );
  }
  fclose(fp_out);


}

//replace the .bam to .txt
char* rename_suffix(const char *fn){
  int l = strlen(fn);
  char *temp = malloc( (size_t)  l * sizeof(char) );
  strcpy(temp, fn);
  strcpy(temp+l-4, ".txt");
  return temp;
}

bool is_bam(char *fname){
  int l = strlen(fname);
  if (l >= 4 && strcasecmp(fname+l-4, ".bam")!=0){
    return false;
  }
  return true;
}

//comparsion string ignoreing case
int strcasecmp(const char *s1, const char *s2){
    int ret;
    
    int l = strlen(s1);
    char *temp = malloc( (size_t) l * sizeof(char));
    strcpy(temp, s1);
  
  int i;
  for(i = 0 ; temp[i]!= '\0'; ++i){
    temp[i] = tolower(temp[i]);
  }

  ret = strcmp(temp, s2);
  return ret;
}