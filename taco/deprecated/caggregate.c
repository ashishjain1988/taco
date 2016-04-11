//
//  caggregate.c
//  TACO
//
//  Created by Balaji Pandian on 2/23/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
//  Need to Relay Safety Checks on sizes of strings, etc.

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <Python.h>

#include "klib/khash.h"

#define MAX_GTF_LINE_SIZE (1024)
#define MAX_TRANSCRIPT_ID_SIZE (256)

typedef struct Hash_node {
    uint32_t num_exons;
    uint32_t length;
    char* id;
    char* transcript_id;
} Hash_node;

KHASH_MAP_INIT_STR(aggregate_hash, Hash_node*)

bool str_in_hashtable(khash_t(aggregate_hash)* hash_table, const char* str) {
    khint_t potential_bucket_id = kh_get(aggregate_hash, hash_table, str);

    if (potential_bucket_id == kh_end(hash_table)) {
        return false;
    }

    return true;
}

void add_bucket_to_hashtable(khash_t(aggregate_hash)* hash_table, const char* transcript_string, const char* id_string) {
    char* allocated_transcript_string = strdup(transcript_string);
    char* allocated_id_string = strdup(id_string);

    Hash_node* new_hash_node = (Hash_node*)malloc(sizeof(Hash_node));
    new_hash_node->num_exons = 0;
    new_hash_node->length = 0;
    new_hash_node->transcript_id = allocated_transcript_string;
    new_hash_node->id = allocated_id_string;

    int ret = -1;
    khint_t bucket_id = kh_put(aggregate_hash, hash_table, allocated_transcript_string, &ret);
    if (ret != 1) {
        printf("CAGGREGATE.C ERROR %d: Could not add new node to hashtable\n", ret);
    }

    kh_val(hash_table, bucket_id) = new_hash_node;
}

char* increment_node_and_return_id(khash_t(aggregate_hash)* hash_table, const char* transcript_id, const uint32_t length_increment) {
    khint_t potential_bucket_id = kh_get(aggregate_hash, hash_table, transcript_id);
    if (potential_bucket_id == kh_end(hash_table)) {
        return NULL;
    }

    Hash_node* node = (Hash_node*)kh_val(hash_table, potential_bucket_id);
    node->num_exons += 1;
    node->length += length_increment;

    return node->id;
}

void free_all_hashtable_nodes(khash_t(aggregate_hash)* hash_table) {
    Hash_node* node;
    kh_foreach_value(hash_table, node, {
        free(node->transcript_id);
        free(node->id);
        free(node);
    });

    kh_destroy(aggregate_hash, hash_table);
}

static PyObject* py_caggregate(PyObject* self, PyObject* args) {
    char* gtf_file;
    char* sample_id;
    char* gtf_expr_attr;
    char* output_file;
    char* stats_file;
    char* is_ref;

    if (!PyArg_ParseTuple(args, "ssssss", &gtf_file, &sample_id, &gtf_expr_attr, &output_file, &stats_file, &is_ref)) {
        return NULL;
    }

    size_t bufflen = MAX_GTF_LINE_SIZE;
    char* buffer = (char*)malloc(bufflen * sizeof(char));
    if (buffer == NULL) {
        return NULL;
    }

    unsigned long long cur_t_id = 1;

    FILE* gtf_file_handler = fopen(gtf_file, "r");
    if (gtf_file_handler == NULL) {
        return NULL;
    }

    FILE* output_file_handler = fopen(output_file, "a");
    if (output_file_handler == NULL) {
        return NULL;
    }

    khash_t(aggregate_hash)* hash_table;
    hash_table = kh_init(aggregate_hash);
    //kh_resize(aggregate_hash, hash_table, 1024);

    while (getline(&buffer, &bufflen, gtf_file_handler) != -1) {
        const char* seqname = strtok(buffer, "\t");
        const char* source = strtok(NULL, "\t");
        const char* feature = strtok(NULL, "\t");
        const char* start = strtok(NULL, "\t");
        const char* end = strtok(NULL, "\t");
        const char* score = strtok(NULL, "\t");
        const char* strand = strtok(NULL, "\t");
        const char* frame = strtok(NULL, "\t");
        const char* attribute = strtok(NULL, "\t");

        // Get Transcript ID from attributes
        char* transcript_id_start_index = strstr(attribute, "transcript_id \"");
        transcript_id_start_index += 15;
        char* end_of_transcript_id = strchr(transcript_id_start_index, '"');
        char transcript_id[MAX_TRANSCRIPT_ID_SIZE];
        strncpy(transcript_id, transcript_id_start_index, end_of_transcript_id - transcript_id_start_index);
        transcript_id[end_of_transcript_id - transcript_id_start_index] = '\0';

        if (strcmp(feature, "transcript") == 0) {
            char id_string[MAX_TRANSCRIPT_ID_SIZE];

            // Search for transcript_id in hashtable
            if (str_in_hashtable(hash_table, transcript_id)) {
                printf("GTF %s transcript_id %s is not unique", gtf_file, transcript_id);
                return NULL;
            } else {
                sprintf(id_string, "%s.T%llu", sample_id, cur_t_id);
                cur_t_id++;

                add_bucket_to_hashtable(hash_table, transcript_id, id_string);
            }

            if (strcmp(is_ref, "True") == 0) {
                fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texpr \"0.0\"; transcript_id \"%s\"; ref \"1\"; sample_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, id_string, sample_id);
            } else {

                // Get GTF_EXPR_ATTR (FPKM, for example) from attributes
                char* expr_start_index = strstr(attribute, gtf_expr_attr);

                // 2 for the space, then quotation mark
                expr_start_index += strlen(gtf_expr_attr) + 2;
                char* end_of_expr_index = strchr(expr_start_index, '"');
                char expr[64];
                strncpy(expr, expr_start_index, end_of_expr_index - expr_start_index);
                expr[end_of_expr_index - expr_start_index] = '\0';

                fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texpr \"%s\"; transcript_id \"%s\"; ref \"0\"; sample_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, expr, id_string, sample_id);
            }

        } else if (strcmp(feature, "exon") == 0) {
            char* id_string  = increment_node_and_return_id(hash_table, transcript_id, (uint32_t) (strtoul(end, NULL, 0) - strtoul(start, NULL, 0)));
            fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, id_string);
        }
    }

    fclose(gtf_file_handler);
    fclose(output_file_handler);
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
