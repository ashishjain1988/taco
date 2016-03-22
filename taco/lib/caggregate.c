//
//  caggregate.c
//  TACO
//
//  Created by Balaji Pandian on 2/23/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
//  Need to Relay Safety Checks on sizes of strings, etc.
//  Need to implement locking on the thread checks.

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <pthread.h>

#include <Python.h>

// Industry-standard FNV-32bit prime and offset
#define FNV_32_PRIME ((uint32_t)0x01000193)
#define FNV_32_OFFSET (2166136261U)

// Hash table size is 2^22 * sizeof(Hash_node) ~100mb
#define HASH_TABLE_SIZE (4194304)
#define MAX_GTF_LINE_SIZE (1024)
#define MAX_TRANSCRIPT_ID_SIZE (256)

// Maximum sie of filepaths including file
#define MAX_FILEPATH_LENGTH (512)

typedef struct Hash_node {
    uint32_t num_exons;
    uint32_t length;
    char* id;
    char* transcript_id;
    struct Hash_node* next;
} Hash_node;

typedef struct File_status_t {
    char raw_gtf_filename[MAX_FILEPATH_LENGTH];
    char aggregated_gtf_filename[MAX_FILEPATH_LENGTH];
    char sorted_gtf_filename[MAX_FILEPATH_LENGTH];
    char sample_id[16];
    char* gtf_expr_attr;
    char* is_ref;
    volatile bool aggregate_started;
    volatile bool aggregate_finished;
    volatile bool sort_started;
    volatile bool sort_finished;
} File_status_t;

typedef struct Thread_arg_struct_t {
    pthread_t pthread;
    volatile bool in_progress;
    File_status_t* file_status;
} Thread_arg_struct_t;


/* This is an industry-standard implementation of the FNV-32bit hash */
uint32_t fnv_32_str(char *str)
{
    unsigned char *s = (unsigned char *)str;
    uint32_t hval = FNV_32_OFFSET;

    while (*s) {
        hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
        hval ^= (uint32_t)*s++;
    }

    return hval % HASH_TABLE_SIZE;
}

bool str_in_hashtable(Hash_node* hash_table_root,char* str) {
    Hash_node bucket = hash_table_root[fnv_32_str(str)];

    if (bucket.transcript_id == NULL) {
        return false;
    }

    if (strcmp(bucket.transcript_id, str) == 0) {
        return true;
    }

    while (bucket.next != NULL) {
        if (strcmp(bucket.transcript_id, str) == 0) {
            return true;
        }

        bucket = *bucket.next;
    }

    return false;
}

void add_bucket_to_hashtable(Hash_node* hash_table_root,
                             char* transcript_string,
                             char* id_string) {
    char* allocated_transcript_string = strdup(transcript_string);
    char* allocated_id_string = strdup(id_string);
    uint32_t hash_value = fnv_32_str(transcript_string);

    if (hash_table_root[hash_value].transcript_id == '\0') {
        // add to initial bucket
        hash_table_root[hash_value].transcript_id = allocated_transcript_string;
        hash_table_root[hash_value].id = allocated_id_string;
        hash_table_root[hash_value].num_exons = 0;
        hash_table_root[hash_value].length = 0;
        hash_table_root[hash_value].next = NULL;
    } else {
        // add to linked list
        Hash_node* bucket = &hash_table_root[hash_value];
        while (bucket->next != NULL) {
            bucket = bucket->next;
        }

        Hash_node* new_hash_node = (Hash_node*)malloc(sizeof(Hash_node));
        new_hash_node->num_exons = 0;
        new_hash_node->length = 0;
        new_hash_node->next = NULL;
        new_hash_node->transcript_id = allocated_transcript_string;
        new_hash_node->id = allocated_id_string;

        bucket->next = new_hash_node;
    }
}

char* increment_node_and_return_id(Hash_node* hash_table_root, char* transcript_id, uint32_t length_increment) {
    Hash_node bucket = hash_table_root[fnv_32_str(transcript_id)];

    if (strcmp(bucket.transcript_id, transcript_id) == 0) {
        // Have the correct bucket.
    } else {
        while (strcmp(bucket.transcript_id, transcript_id) != 0) {
            bucket = *(bucket.next);
        }
    }

    bucket.num_exons += 1;
    bucket.length += length_increment;

    return bucket.id;
}

void free_all_hashtable_nodes(Hash_node* hash_table_root) {
    uint32_t i;
    for (i = 0; i < HASH_TABLE_SIZE; i++) {
        Hash_node* temp = &hash_table_root[i];

        if (temp->transcript_id != NULL) {
            free(temp->transcript_id);
            free(temp->id);
        }

        // Jump past first bucket -- allocated all at once in beginning
        temp = temp->next;

        while (temp != NULL) {
            Hash_node* temp_next = temp->next;
            free(temp->transcript_id);
            free(temp->id);
            free(temp);
            temp = temp_next;
        }
    }

    free(hash_table_root);
}

int aggregate_function(char* gtf_file, char* sample_id, char* gtf_expr_attr, char* output_file, char* is_ref) {
    size_t bufflen = MAX_GTF_LINE_SIZE;
    char* buffer = (char*)malloc(bufflen * sizeof(char));
    if (buffer == NULL) {
        return 1;
    }

    unsigned long long cur_t_id = 1;

    FILE* gtf_file_handler = fopen(gtf_file, "r");
    if (gtf_file_handler == NULL) {
        return 1;
    }

    FILE* output_file_handler = fopen(output_file, "a");
    if (output_file_handler == NULL) {
        return 1;
    }

    Hash_node* t_dict = (Hash_node*)calloc(HASH_TABLE_SIZE, sizeof(Hash_node));
    if (t_dict == NULL) {
        return 1;
    }

    while (getline(&buffer, &bufflen, gtf_file_handler) != -1) {
        char* saveptr;
        const char* seqname = strtok_r(buffer, "\t", &saveptr);
        const char* source = strtok_r(NULL, "\t", &saveptr);
        const char* feature = strtok_r(NULL, "\t", &saveptr);
        const char* start = strtok_r(NULL, "\t", &saveptr);
        const char* end = strtok_r(NULL, "\t", &saveptr);
        const char* score = strtok_r(NULL, "\t", &saveptr);
        const char* strand = strtok_r(NULL, "\t", &saveptr);
        const char* frame = strtok_r(NULL, "\t", &saveptr);
        const char* attribute = strtok_r(NULL, "\t", &saveptr);

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
            if (str_in_hashtable(t_dict, transcript_id)) {
                printf("GTF %s transcript_id %s is not unique", gtf_file, transcript_id);
                return 1;
            } else {
                sprintf(id_string, "%s.T%llu", sample_id, cur_t_id);
                cur_t_id++;

                add_bucket_to_hashtable(t_dict, transcript_id, id_string);
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
            char* id_string  = increment_node_and_return_id(t_dict, transcript_id, (uint32_t) (strtoul(end, NULL, 0) - strtoul(start, NULL, 0)));
            fprintf(output_file_handler, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id \"%s\";\n", seqname, source, feature, start, end, score, strand, frame, id_string);
        }
    }

    fclose(gtf_file_handler);
    fclose(output_file_handler);
    free_all_hashtable_nodes(t_dict);
    free(buffer);

    return 0;
}

void* aggregate_thread(void* packed_args) {
    Thread_arg_struct_t* args = (Thread_arg_struct_t*)packed_args;
    File_status_t* file_status = args->file_status;
    if (aggregate_function(file_status->raw_gtf_filename, file_status->sample_id, file_status->gtf_expr_attr, file_status->aggregated_gtf_filename, file_status->is_ref) != 0) {
        printf("Aggregate Function Failed\n");
    }

    args->in_progress = false;
    file_status->aggregate_finished = true;
    return NULL;
}

static PyObject* py_caggregate(PyObject* self, PyObject* args) {
    PyObject* listObj;
    PyObject* strObj;
    char* gtf_expr_attr;
    char* tmp_dir;
    char* is_ref;
    int cpu_count;

    if (!PyArg_ParseTuple(args, "O!sssb", &PyList_Type, &listObj, &gtf_expr_attr, &tmp_dir, &is_ref, &cpu_count)) {
        return NULL;
    }

    int aggregate_cpu_count = 4;
    int sort_cpu_count = 6;

    int num_files = (int)PyList_Size(listObj);
    if (num_files <= 0) {
        return NULL;
    }

    cpu_count--; // Account for this main thread.
    if (cpu_count < 1) {
        cpu_count = 1;
    }

    File_status_t* file_status = (File_status_t*)calloc(num_files, sizeof(File_status_t));
    Thread_arg_struct_t* aggregate_thread_status = (Thread_arg_struct_t*)calloc(aggregate_cpu_count, sizeof(Thread_arg_struct_t));
    Thread_arg_struct_t* sort_thread_status = (Thread_arg_struct_t*)calloc(sort_cpu_count, sizeof(Thread_arg_struct_t));

    for (uint32_t i = 0; i < num_files; i++) {
        strObj = PyList_GetItem(listObj, i);
        char* input_filename = PyString_AsString(strObj);

        char sample_id[16];
        sprintf(sample_id, "%u", i);

        char output_file[strlen(tmp_dir) + 32];
        sprintf(output_file, "%s/transcripts.%u.gtf", tmp_dir, i);

        char sorted_file[strlen(tmp_dir) + 64];
        sprintf(sorted_file, "%s/transcripts.%u.sorted.gtf", tmp_dir, i);

        strncpy(file_status[i].raw_gtf_filename, input_filename, MAX_FILEPATH_LENGTH);
        strncpy(file_status[i].aggregated_gtf_filename, output_file, MAX_FILEPATH_LENGTH);
        strncpy(file_status[i].sorted_gtf_filename, sorted_file, MAX_FILEPATH_LENGTH);
        strncpy(file_status[i].sample_id, sample_id, 16);
        file_status[i].gtf_expr_attr = gtf_expr_attr;
        file_status[i].is_ref = is_ref;
    }

    while (true) {

        for (uint32_t thread_id = 0; thread_id < aggregate_cpu_count; thread_id++) {
            if (aggregate_thread_status[thread_id].in_progress == false) {
                aggregate_thread_status[thread_id].in_progress = true;

                bool new_thread_started = false;
                for (uint32_t i = 0; i < num_files; i++) {
                    if (file_status[i].aggregate_started == false) {
                        file_status[i].aggregate_started = true;
                        new_thread_started = true;

                        aggregate_thread_status[thread_id].file_status = &file_status[i];
                        if (pthread_create((pthread_t*)&aggregate_thread_status[thread_id].pthread, NULL, aggregate_thread, (void*)&aggregate_thread_status[thread_id])) {
                            fprintf(stderr, "Error creating thread\n");
                            return NULL;
                        }
                        break;
                    }
                }

                if (new_thread_started == false) {
                    aggregate_thread_status[thread_id].in_progress = false;
                }
            }
        }

        
        for (uint32_t thread_id = 0; thread_id < sort_cpu_count; thread_id++) {
            if (sort_thread_status[thread_id].in_progress == false) {
                sort_thread_status[thread_id].in_progress = true;

                bool new_thread_started = false;
                for (uint32_t i = 0; i < num_files; i++) {
                    if ((file_status[i].aggregate_finished) && (file_status[i].sort_started == false)) {
                        file_status[i].sort_started = true;
                        new_thread_started = true;

                        sort_thread_status[thread_id].file_status = &file_status[i];
                        char sort_command[2 * strlen(file_status[i].aggregated_gtf_filename) + strlen(file_status[i].sorted_gtf_filename) + 128];
                        sprintf(sort_command, "( env LC_ALL=C sort -k1,1 -k4,4n -k3,3r -o %s %s ; rm %s ) &", file_status[i].sorted_gtf_filename, file_status[i].aggregated_gtf_filename, file_status[i].aggregated_gtf_filename);

                        system(sort_command);
                        break;
                    }
                }

                if (new_thread_started == false) {
                    sort_thread_status[thread_id].in_progress = false;
                }
            }
        }
        
        for (uint32_t thread_id = 0; thread_id < sort_cpu_count; thread_id++) {
            if (sort_thread_status[thread_id].in_progress == true) {
                if ((sort_thread_status[thread_id].file_status->sort_started == true) && (access((char*)sort_thread_status[thread_id].file_status->aggregated_gtf_filename, F_OK) == -1)) {
                    // original file doesn't exist anymore... the thread is done.
                    sort_thread_status[thread_id].file_status->sort_finished = true;
                    sort_thread_status[thread_id].in_progress = false;
                }
            }

        }
        
        bool all_files_completed = true;
        for (uint32_t i = 0; i < num_files; i++) {
            if (file_status[i].sort_finished == false) {
                all_files_completed = false;
                break;
            }
        }
        if (all_files_completed) {
            break;
        }

        sleep(2);
    }

    // Final Merge
    memset(sort_thread_status, 0, sort_cpu_count * sizeof(Thread_arg_struct_t));
    int num_layer_1_merge_files = num_files / sort_cpu_count;

    // The thread_id indexing goes to -1 because the last sort thread needs to pick up the remainder
    for (uint32_t thread_id = 0; thread_id < sort_cpu_count - 1; thread_id++) {

        char* merge_command = (char*)calloc(num_layer_1_merge_files + 1, sizeof(file_status[0].sorted_gtf_filename) + 16);
        sprintf(merge_command, "( env LC_ALL=C sort -T %s -m -k1,1 -k4,4n -k3,3r -o %s/merge.%u.gtf ", tmp_dir, tmp_dir, thread_id);
        unsigned long merge_command_length = strlen(merge_command);
        for (uint32_t i = (thread_id * num_layer_1_merge_files); i < ((thread_id * num_layer_1_merge_files) + num_layer_1_merge_files); i++) {
            merge_command_length += sprintf(merge_command + merge_command_length, "%s ", file_status[i].sorted_gtf_filename);
        }
        merge_command_length += sprintf(merge_command + merge_command_length, "; rm ");
        for (uint32_t i = (thread_id * num_layer_1_merge_files); i < ((thread_id * num_layer_1_merge_files) + num_layer_1_merge_files); i++) {
            merge_command_length += sprintf(merge_command + merge_command_length, "%s ", file_status[i].sorted_gtf_filename);
        }
        merge_command_length += sprintf(merge_command + merge_command_length, ") &");
        system(merge_command);
    }

    // Get the last merge to pick up any remainders
    char* merge_command = (char*)calloc(num_layer_1_merge_files * 2, sizeof(file_status[0].sorted_gtf_filename) + 16);
    sprintf(merge_command, "( env LC_ALL=C sort -T %s -m -k1,1 -k4,4n -k3,3r -o %s/merge.%u.gtf ", tmp_dir, tmp_dir, sort_cpu_count - 1);
    unsigned long merge_command_length = strlen(merge_command);
    for (uint32_t i = ((sort_cpu_count - 1) * num_layer_1_merge_files); i < num_files; i++) {
        merge_command_length += sprintf(merge_command + merge_command_length, "%s ", file_status[i].sorted_gtf_filename);
    }
    merge_command_length += sprintf(merge_command + merge_command_length, "; rm ");
    for (uint32_t i = ((sort_cpu_count - 1) * num_layer_1_merge_files); i < num_files; i++) {
        merge_command_length += sprintf(merge_command + merge_command_length, "%s ", file_status[i].sorted_gtf_filename);
    }
    merge_command_length += sprintf(merge_command + merge_command_length, ") &");
    system(merge_command);

    // Wait for 1st layer merge to finish
    while (true) {
        bool layer_1_merge_completed = true;
        for (uint32_t i = 0; i < num_files; i++) {
            if (access((char*)file_status[i].sorted_gtf_filename, F_OK) == -1) {
                // File does not exist... Means that thread is done
            } else {
                layer_1_merge_completed = false;
                break;
            }

        }

        if (layer_1_merge_completed) {
            break;
        } else {
            sleep(2);
        }
    }

    // Final layer 2 merge.
    sprintf(merge_command, "env LC_ALL=C sort -T %s -m -k1,1 -k4,4n -k3,3r -o %s/transfrags.gtf ", tmp_dir, tmp_dir);
    merge_command_length = strlen(merge_command);
    for (uint32_t i = 0; i < sort_cpu_count; i++) {
        merge_command_length += sprintf(merge_command + merge_command_length, "%s/merge.%u.gtf ", tmp_dir, i);
    }
    merge_command_length += sprintf(merge_command + merge_command_length, "; rm ");
    for (uint32_t i = 0; i < sort_cpu_count; i++) {
        merge_command_length += sprintf(merge_command + merge_command_length, "%s/merge.%u.gtf ", tmp_dir, i);
    }
    system(merge_command);

    free(file_status);
    free(aggregate_thread_status);
    free(sort_thread_status);

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
