//  caggregate.c
//  TACO
//
//  Created by Balaji Pandian on 4/15/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
//

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <malloc.h>
#include <Python.h>

#include "klib/khash.h"

#define MAX_GTF_LINE_SIZE (1024)
#define MAX_TRANSCRIPT_ID_SIZE (256)

typedef struct Hash_node {
    unsigned long long* exon_start_array;
    unsigned long long* exon_end_array;
    unsigned long long exon_array_length;
    unsigned long long exon_full_array_length;
    char* chrom;
    char strand;
    bool is_ref_val;
    char* id;
    char* transcript_id;
    float expr;
} Hash_node;

KHASH_MAP_INIT_STR(aggregate_hash, Hash_node*)

bool str_in_hashtable(khash_t(aggregate_hash)* hash_table, const char* str) {
    khint_t potential_bucket_id = kh_get(aggregate_hash, hash_table, str);

    if (potential_bucket_id == kh_end(hash_table)) {
        return false;
    }

    return true;
}

void add_bucket_to_hashtable(khash_t(aggregate_hash)* hash_table, const char* transcript_string, const char* id_string, const char* seqname, const char strand, const double expr, const bool is_ref_val) {
    char* allocated_transcript_string = strdup(transcript_string);
    char* allocated_id_string = strdup(id_string);
    char* allocated_chrom_string = strdup(seqname);

    Hash_node* new_hash_node = (Hash_node*)malloc(sizeof(Hash_node));
    new_hash_node->exon_array_length = 0;
    new_hash_node->exon_full_array_length = 0;
    new_hash_node->exon_start_array = NULL;
    new_hash_node->exon_end_array = NULL;
    new_hash_node->chrom = allocated_chrom_string;
    new_hash_node->strand = strand;
    new_hash_node->is_ref_val = is_ref_val;
    new_hash_node->transcript_id = allocated_transcript_string;
    new_hash_node->id = allocated_id_string;
    new_hash_node->expr = expr;

    int ret = -1;
    khint_t bucket_id = kh_put(aggregate_hash, hash_table, allocated_transcript_string, &ret);
    if (ret != 1) {
        printf("CAGGREGATE.C ERROR %d: Could not add new node to hashtable\n", ret);
    }

    kh_val(hash_table, bucket_id) = new_hash_node;
}

int append_exon_pair_to_bucket(khash_t(aggregate_hash)* hash_table, const char* transcript_id, const unsigned long long exon_start, const unsigned long long exon_end) {
    khint_t potential_bucket_id = kh_get(aggregate_hash, hash_table, transcript_id);
    if (potential_bucket_id == kh_end(hash_table)) {
        return 1;
    }

    Hash_node* node = (Hash_node*)kh_val(hash_table, potential_bucket_id);

    if (node->exon_array_length >= node->exon_full_array_length) {
        // need to resize exon array (x2 size)
        if (node->exon_full_array_length == 0) {
            node->exon_start_array = malloc(8 * sizeof(unsigned long long));
            node->exon_end_array = malloc(8 * sizeof(unsigned long long));
            node->exon_full_array_length = 8;
        } else {
            node->exon_start_array = realloc(node->exon_start_array, node->exon_full_array_length * 2 * sizeof(unsigned long long));
            node->exon_end_array = realloc(node->exon_end_array, node->exon_full_array_length * 2 * sizeof(unsigned long long));
            node->exon_full_array_length *= 2;
        }
        
        if ((node->exon_start_array == NULL) || (node->exon_end_array == NULL)){
            return 1;
        }
    }

    node->exon_start_array[node->exon_array_length] = exon_start;
    node->exon_end_array[node->exon_array_length] = exon_end;
    node->exon_array_length++;

    return 0;
}

void free_all_hashtable_nodes(khash_t(aggregate_hash)* hash_table) {
    Hash_node* node;
    kh_foreach_value(hash_table, node, {
        free(node->transcript_id);
        free(node->id);
        free(node->exon_start_array);
        free(node->exon_end_array);
        free(node->chrom);
        free(node);
    });

    kh_destroy(aggregate_hash, hash_table);
}

static PyObject* py_caggregate(PyObject* self, PyObject* args) {
    char* gtf_file;
    char* sample_id;
    char* gtf_expr_attr;
    char* is_ref;
    char* bed_filename;
    char* filtered_bed_filename;

    if (!PyArg_ParseTuple(args, "ssssss", &gtf_file, &sample_id, &gtf_expr_attr, &is_ref, &bed_filename, &filtered_bed_filename)) {
        return NULL;
    }

    size_t bufflen = MAX_GTF_LINE_SIZE;
    char* buffer = (char*)malloc(bufflen * sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    FILE* gtf_file_handler = fopen(gtf_file, "r");
    if (gtf_file_handler == NULL) {
        return NULL;
    }

    FILE* bed_file_handler = fopen(bed_filename, "w");
    if (bed_file_handler == NULL) {
        return NULL;
    }

    FILE* filtered_bed_file_handler = fopen(filtered_bed_filename, "w");
    if (filtered_bed_file_handler == NULL) {
        return NULL;
    }

    khash_t(aggregate_hash)* hash_table;
    hash_table = kh_init(aggregate_hash);

    unsigned long long cur_transcript_id = 1;
    double total_expr = 0.0;

    while (getline(&buffer, &bufflen, gtf_file_handler) != -1) {
        char* saveptr = NULL; 
        const char* seqname = strtok_r(buffer, "\t", &saveptr);
        const char* source = strtok_r(NULL, "\t", &saveptr);
        const char* feature = strtok_r(NULL, "\t", &saveptr);
        const char* start = strtok_r(NULL, "\t", &saveptr);
        const char* end = strtok_r(NULL, "\t", &saveptr);
        const char* score = strtok_r(NULL, "\t", &saveptr);
        const char* strand = strtok_r(NULL, "\t", &saveptr);
        const char* frame = strtok_r(NULL, "\t", &saveptr);
        const char* attribute = strtok_r(NULL, "\t", &saveptr);

        (void)frame;
        (void)score;
        (void)source;

        // Get Transcript ID from attributes
        char* transcript_id_start_index = strstr(attribute, "transcript_id \"");
        transcript_id_start_index += 15;
        char* end_of_transcript_id = strchr(transcript_id_start_index, '"');
        char transcript_id[MAX_TRANSCRIPT_ID_SIZE];
        strncpy(transcript_id, transcript_id_start_index, end_of_transcript_id - transcript_id_start_index);
        transcript_id[end_of_transcript_id - transcript_id_start_index] = '\0';

        if (strcmp(feature, "transcript") == 0) {
            if (str_in_hashtable(hash_table, transcript_id)) {
                printf("Transcript %s duplicate detected\n", transcript_id);
                return NULL;
            }

            char new_transcript_id[MAX_TRANSCRIPT_ID_SIZE];
            sprintf(new_transcript_id, "%s.%llu", sample_id, cur_transcript_id);
            cur_transcript_id++;

            double expr = 0.0;
            bool is_ref_val;
            if (strcmp(is_ref, "True") == 0) {
                is_ref_val = true;
                expr = 0.0;
            } else {
                is_ref_val = false;
                if (strstr(attribute, gtf_expr_attr) == NULL) {
                    printf("GTF expression attribute %s not found\n", gtf_expr_attr);
                    return NULL;
                }

                // Get GTF_EXPR_ATTR (FPKM, for example) from attributes
                char* expr_start_index = strstr(attribute, gtf_expr_attr);
                // Add 2 for space then quotation mark
                expr_start_index += strlen(gtf_expr_attr) + 2;
                char* end_of_expr_index = strchr(expr_start_index, '"');
                char expr_string[64];
                strncpy(expr_string, expr_start_index, end_of_expr_index - expr_start_index);
                expr_string[end_of_expr_index - expr_start_index] = '\0';

                expr = atof(expr_string);
                total_expr += expr;
            }

            add_bucket_to_hashtable(hash_table, transcript_id, new_transcript_id, seqname, *strand, expr, is_ref_val);
        } else if (strcmp(feature, "exon") == 0) {
            if (str_in_hashtable(hash_table, transcript_id) == false) {
                printf("Transcript %s exon feature appeared in gtf file prior to transcript feature\n", transcript_id);
                return NULL;
            }

            // The -1 is to convert from a 1-based (inclusive) to 0-based (exclusive) intervals
            if (append_exon_pair_to_bucket(hash_table, transcript_id, atof(start) - 1, atof(end)) != 0) {
                return NULL;
            }
        }
    }

    unsigned long long nlength = 0;
    unsigned long long nexpr = 0;
    unsigned long long nsplice = 0;

    Hash_node* node;
    kh_foreach_value(hash_table, node, {
        if (strcmp(node->id, "1.2") != 0) {
            continue;
        }
        printf("%s\n", node->id);
        printf("%f\n", node->expr);
        for (unsigned long long i = 0; i < node->exon_array_length; i++) {
            printf("S: %llu\n", node->exon_start_array[i]);
            printf("E: %llu\n", node->exon_end_array[i]);
        }
        printf("\n");
        return NULL;
    });

    fclose(gtf_file_handler);
    fclose(bed_file_handler);
    fclose(filtered_bed_file_handler);
    free_all_hashtable_nodes(hash_table);
    free(buffer);

    // Return value -- returning "NULL" is an error
    return Py_BuildValue("i", 0);
}

// Mapping between python and c function names
static PyMethodDef CAggregateMethods[] = {
    {"caggregate", py_caggregate, METH_VARARGS},
    {NULL, NULL}
};

// Module initialisation routine
PyMODINIT_FUNC
initcaggregate(void)
{
    (void) Py_InitModule("caggregate", CAggregateMethods);
}
