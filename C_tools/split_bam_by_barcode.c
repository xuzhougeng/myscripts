#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <time.h>
#include <pthread.h>

#define NUM_RECORD 1000000
#define MAX_BARCODE_NUM 100000
#define INITIAL_CAPACITY 5000

typedef struct {
    char **barcodes;
    int count;
} Barcodes;


Barcodes load_barcodes(const char *filename) {
    FILE *file = fopen(filename, "r");
    char line[256];
    Barcodes barcodes;
    barcodes.barcodes = malloc(MAX_BARCODE_NUM * sizeof(char*));
    barcodes.count = 0;

    while (fgets(line, sizeof(line), file)) {
        int len = strlen(line);
        if (line[len - 1] == '\n') {
            line[len - 1] = '\0';
        }
        barcodes.barcodes[barcodes.count] = strdup(line);
        barcodes.count++;
    }

    fclose(file);
    return barcodes;
}

typedef struct {
    bam1_t ****records; // [barcode][chrom][record] = bam1_t*
    int **counts;
    int **capacities;
    pthread_mutex_t *locks;
} BamStorage;

BamStorage init_bam_storage(int barcode_num, int chrom_num) {
    BamStorage storage;
    storage.records = malloc(barcode_num * sizeof(bam1_t**));
    storage.counts = malloc(barcode_num * sizeof(int*));
    storage.capacities = malloc(barcode_num * sizeof(int*));
    storage.locks = malloc(barcode_num * sizeof(pthread_mutex_t));

    for (int i = 0; i < barcode_num; i++) {
        storage.records[i] = malloc(chrom_num * sizeof(bam1_t*));
        storage.counts[i] = malloc(chrom_num * sizeof(int));
        storage.capacities[i] = malloc(chrom_num * sizeof(int));
        pthread_mutex_init(&storage.locks[i], NULL);
        for (int j = 0; j < chrom_num; j++) {
            storage.records[i][j] = malloc(INITIAL_CAPACITY * sizeof(bam1_t*));
            storage.counts[i][j] = 0;
            storage.capacities[i][j] = INITIAL_CAPACITY;
        }
    }

    return storage;
}

typedef struct {
    const char *region;
    const char *input_bam;
    Barcodes barcodes;
    BamStorage *storage;
    int chrom_idx;
} thread_data_t;

void process_region(thread_data_t *data) {
    htsFile *in = sam_open(data->input_bam, "r");
    bam_hdr_t *hdr = sam_hdr_read(in);
    hts_idx_t *idx = sam_index_load(in, data->input_bam);
    if (!idx) {
        fprintf(stderr, "Error: BAM file is not indexed.\n");
        sam_close(in);
        bam_hdr_destroy(hdr);
        return;
    }

    hts_itr_t *iter = sam_itr_querys(idx, hdr, data->region);
    if (!iter) {
        fprintf(stderr, "Error: Could not parse region %s.\n", data->region);
        sam_close(in);
        bam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        return;
    }

    bam1_t *b = bam_init1();
    BamStorage *storage = data->storage;
    int chrom_idx = data->chrom_idx;

    while (sam_itr_next(in, iter, b) >= 0) {
        uint8_t *cb = bam_aux_get(b, "CB");
        if (cb != NULL) {
            char *barcode = bam_aux2Z(cb);
            for (int i = 0; i < data->barcodes.count; i++) {
                                
                if (strcmp(barcode, data->barcodes.barcodes[i]) == 0) {
                    pthread_mutex_lock(&storage->locks[i]);
                    if (storage->counts[i][chrom_idx] == storage->capacities[i][chrom_idx]) {
                        storage->capacities[i][chrom_idx] *= 2;
                        storage->records[i][chrom_idx] = realloc(storage->records[i][chrom_idx], storage->capacities[i][chrom_idx] * sizeof(bam1_t*));
                    }
                    storage->records[i][chrom_idx][storage->counts[i][chrom_idx]] = bam_dup1(b);
                    storage->counts[i][chrom_idx]++;
                    pthread_mutex_unlock(&storage->locks[i]);
                    break;
                }
            }
        }
    }

    sam_itr_destroy(iter);
    hts_idx_destroy(idx);
    sam_close(in);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
}

void *threaded_process_region(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    process_region(data);
    free(data);
    return NULL;
}

void free_bam_storage(BamStorage storage, int barcode_num, int chrom_num) {
    for (int i = 0; i < barcode_num; i++) {
        for (int j = 0; j < chrom_num; j++) {
            for (int k = 0; k < storage.counts[i][j]; k++) {
                bam_destroy1(storage.records[i][j][k]);
            }
            free(storage.records[i][j]);
        }
        free(storage.records[i]);
        free(storage.counts[i]);
        free(storage.capacities[i]);
        pthread_mutex_destroy(&storage.locks[i]);
    }

    free(storage.records);
    free(storage.counts);
    free(storage.capacities);
    free(storage.locks);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_bam> <barcode_file>\n", argv[0]);
        return 1;
    }

    Barcodes barcodes = load_barcodes(argv[2]);

    htsFile *in = sam_open(argv[1], "r");
    bam_hdr_t *hdr = sam_hdr_read(in);
    BamStorage storage = init_bam_storage(barcodes.count, hdr->n_targets);

    pthread_t threads[hdr->n_targets];

    for (int i = 0; i < hdr->n_targets; i++) {
        thread_data_t *data = malloc(sizeof(thread_data_t));
        data->region = hdr->target_name[i];
        data->input_bam = argv[1];
        data->barcodes = barcodes;
        data->storage = &storage;
        data->chrom_idx = i;

        pthread_create(&threads[i], NULL, threaded_process_region, (void *)data);
    }

    for (int i = 0; i < hdr->n_targets; i++) {
        pthread_join(threads[i], NULL);
    }

    // Do something with the storage ...

    free_bam_storage(storage, barcodes.count, hdr->n_targets);

    for (int i = 0; i < barcodes.count; i++) {
        free(barcodes.barcodes[i]);
    }
    free(barcodes.barcodes);
    bam_hdr_destroy(hdr);
    sam_close(in);

    return 0;
}