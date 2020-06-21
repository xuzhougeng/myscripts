#include "./kseq.h"
#include <zlib.h>
#include <string.h>
#include <stdio.h>

KSEQ_INIT(gzFile, gzread)

const char *get_filename_ext(const char *filename){
	const char *dot = strrchr(filename, '.');
	if ( !dot || dot == filename ) return "";
	return dot + 1;
}

void split_file(const char *file, const char* ext){
	gzFile fi;
	FILE *fo;
	kseq_t *seq;

	fi = gzopen(file, "r");
	seq = kseq_init(fi);
	int l;
	char *name;
	char *token;
	while ( ( l = kseq_read(seq)) >= 0) {

		name = strdup(seq->name.s);
		token = strtok(name, ":");
		strcat(token, ".");
		strcat(token, ext);
		fo = fopen(token, "a+");
		fprintf(fo, ">%s\n%s\n", seq->name.s, seq->seq.s);
		fclose(fo);

		free(name);
	
	}

	kseq_destroy(seq);
	gzclose(fi);
	return ;

}

int main(int argc, char *argv[])
{
	if (argc == 1){
		fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
		return 1;
	}
	char *file_name = argv[1];
	const char* ext = get_filename_ext(file_name);
	split_file(file_name, ext);

	return 0;
}
