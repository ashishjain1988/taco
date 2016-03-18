//
//  csort.c
//  TACO
//
//  Created by Balaji Pandian on 3/15/16.
//  Copyright © 2016 Balaji Pandian. All rights reserved.
//
//  TODO:
// 

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <pthread.h>
#include <Python.h>

#define MAX_GTF_LINE_SIZE (1024)

typedef struct Arg_struct {
    pthread_t pthread;
    volatile bool in_progress;
    char* filename_a;
    char* filename_b;
    char output_filename[512];
} Arg_struct;

inline int gtf_string_compare(const char* restrict line_a, const char* restrict line_b) {

    char feature_a[32];
    char feature_b[32];
    char chr_a[4];
    char chr_b[4];
    uint32_t start_a;
    uint32_t start_b;
    uint32_t end_a;
    uint32_t end_b;

    sscanf(line_a, "chr%s\t%*s\t%s\t%u\t%u\t%*s", chr_a, feature_a, &start_a, &end_a);
    sscanf(line_b, "chr%s\t%*s\t%s\t%u\t%u\t%*s", chr_b, feature_b, &start_b, &end_b);

    // accurate "chr" sort
    /*
    if (((chr_a[0] == 'x') || (chr_a[0] == 'X')) && ((chr_b[0] == 'y') || (chr_b[0] == 'Y'))) {
        return -1;
    } else if (((chr_a[0] == 'y') || (chr_a[0] == 'Y')) && ((chr_b[0] == 'x') || (chr_b[0] == 'X'))) {
        return 1;
    } else if ((chr_a[1] == '\0') && (chr_b[1] == '\0')) {
        if (chr_a[0] < chr_b[0]) {
            return -1;
        } else if (chr_a[0] > chr_b[0]) {
            return 1;
        }
    } else if ((chr_a[1] != '\0') && (chr_b[1] != '\0')) {
        if (chr_a[0] < chr_b[0]) {
            return -1;
        } else if (chr_a[0] > chr_b[0]) {
            return 1;
        } else if (chr_a[1] < chr_b[1]) {
            return -1;
        } else if (chr_a[1] > chr_b[1]) {
            return 1;
        }
    } else if ((chr_a[1] == '\0') && (chr_b[1] != '\0')) {
        return -1;
    } else if ((chr_a[1] != '\0') && (chr_b[1] == '\0')) {
        return 1;
    }*/

    // (intentionally) Incorrect "chr" sort
    if (chr_a[0] < chr_b[0]) {
        return -1;
    } else if (chr_a[0] > chr_b[0]) {
        return 1;
    } else {
        if (chr_a[1] < chr_b[1]) {
            return -1;
        } else if (chr_a[1] > chr_b[1]) {
            return 1;
        }
    }

    // BP Start Sort
    if (start_a < start_b) {
        return -1;
    } else if (start_a > start_b) {
        return 1;
    }

    // Feature sort
    if (feature_a[0] == 't') {
        return -1;
    } else if (feature_b[0] == 't') {
        return 1;
    }

    // BP End sort
    if (end_a < end_b) {
        return -1;
    } else if (end_a > end_b) {
        return 1;
    }

    return 1;
}

void merge_two_files_and_delete_originals(const char* file_name_a, const char* file_name_b, const char* output_filename) {

    if ((access(file_name_a, F_OK) == -1) && (access(file_name_b, F_OK) != -1)) {
        // File B exists, File A doesn't.
        rename(file_name_b, output_filename);
        return;
    } else if ((access(file_name_a, F_OK) != -1) && (access(file_name_b, F_OK) == -1)) {
        // File A exists, File B doesn't.
        rename(file_name_a, output_filename);
        return;
    }

    FILE* file_a = fopen(file_name_a, "r");
    FILE* file_b = fopen(file_name_b, "r");
    FILE* output_file = fopen(output_filename, "w");

    char line_a[MAX_GTF_LINE_SIZE];
    char line_b[MAX_GTF_LINE_SIZE];

    bool file_a_has_data = (fgets(line_a, MAX_GTF_LINE_SIZE, file_a) != NULL);
    bool file_b_has_data = (fgets(line_b, MAX_GTF_LINE_SIZE, file_b) != NULL);

    while ((file_a_has_data) && (file_b_has_data)) {
        if (gtf_string_compare(line_a, line_b) < 0) {
            fputs(line_a, output_file);
            file_a_has_data = (fgets(line_a, MAX_GTF_LINE_SIZE, file_a) != NULL);
        } else {
            fputs(line_b, output_file);
            file_b_has_data = (fgets(line_b, MAX_GTF_LINE_SIZE, file_b) != NULL);

        }
    }

    // Dump rest of other file -- one file is empty
    if ((file_a_has_data) && (file_b_has_data == false)) {
        while (file_a_has_data) {
            fputs(line_a, output_file);
            file_a_has_data = (fgets(line_a, MAX_GTF_LINE_SIZE, file_a) != NULL);
        }
    } else if ((file_a_has_data == false) && (file_b_has_data)) {
        while (file_b_has_data) {
            fputs(line_b, output_file);
            file_b_has_data = (fgets(line_b, MAX_GTF_LINE_SIZE, file_b) != NULL);
        }
    }

    fclose(file_a);
    fclose(file_b);
    fclose(output_file);

    remove(file_name_a);
    remove(file_name_b);
}

void* merge_thread(void* packed_args) {
    Arg_struct* args = (Arg_struct*)packed_args;
    char file_name_a[strlen(args->filename_a) + sizeof(char)];
    char file_name_b[strlen(args->filename_b) + sizeof(char)];
    char output_filename[strlen(args->output_filename) + sizeof(char)];

    strcpy(file_name_a, args->filename_a);
    strcpy(file_name_b, args->filename_b);
    strcpy(output_filename, args->output_filename);

    merge_two_files_and_delete_originals(file_name_a, file_name_b, output_filename);
    args->in_progress = false;
    return NULL;
}

static PyObject* py_csort(PyObject* self, PyObject* args) {
    int num_files;

    PyObject* listObj;
    PyObject* strObj;
    char* tmp_dir;
    int cpu_count;

    if (!PyArg_ParseTuple(args, "O!sb", &PyList_Type, &listObj, &tmp_dir, &cpu_count)) {
        return NULL;
    }

    num_files = (int)PyList_Size(listObj);
    if (num_files < 0) {
        return NULL;
    }

    cpu_count--; // Account for this main thread.
    if (cpu_count < 1) {
        cpu_count = 1;
    }

    volatile Arg_struct thread_status[cpu_count];
    for (int i = 0; i < cpu_count; i++) {
        thread_status[i].in_progress = false;
    }

    char* filename_array[num_files];
    for (uint32_t i = 0; i < num_files; i++) {
        strObj = PyList_GetItem(listObj, i);
        char* line = PyString_AsString(strObj);

        char sorted_filename[strlen(line) + 16];
        sprintf(sorted_filename, "%s.sorted", line);

        char sort_command[2*strlen(line) + strlen(sorted_filename) + 128];
        sprintf(sort_command, "( sort -k1,1 -k4,4n -k3,3r %s > %s ; rm %s ) &", line, sorted_filename, line);

        while (true) {
            bool created_thread = false;
            for (int thread_id = 0; thread_id < cpu_count; thread_id++) {
                if (thread_status[thread_id].in_progress == false) {
                    thread_status[thread_id].in_progress = true;
                    
                    strcpy((char*)thread_status[thread_id].output_filename, line); 

                    system(sort_command);

                    created_thread = true;
                    break;
                }
            }

            if (created_thread) {
                break;
            } else {
                // wait for other threads to finish
                for (int thread_id = 0; thread_id < cpu_count; thread_id++) {
                    if (access((char*)thread_status[thread_id].output_filename, F_OK) == -1) {
                        // original file doesn't exist anymore... the thread is done.
                        thread_status[thread_id].in_progress = false;
                    }
                }
            }

            sleep(1);
        }

        filename_array[i] = malloc(strlen(sorted_filename) + sizeof(char)); // Extra byte for null-terminator
        strcpy(filename_array[i], sorted_filename);
    }

    // Wait for all threads to finish
    while (true){
        bool all_threads_clear = true;
        for (int thread_id = 0; thread_id < cpu_count; thread_id++) {
            if (access((char*)thread_status[thread_id].output_filename, F_OK) == -1) {
                thread_status[thread_id].in_progress = false;
            } else {
                all_threads_clear = false;
            }
        }

        if (all_threads_clear) {
            break;
        }
    }

    // Reset threads
    for (int i = 0; i < cpu_count; i++) {
        thread_status[i].in_progress = false;
    }

    // Merge files in a binary pattern
    uint32_t files_to_merge = num_files;
    while (files_to_merge > 1) {
        if ((files_to_merge % 2) == 0) {
            // Even number of files to merge.
            files_to_merge /= 2;
        } else {
            // Odd nuumber of files to merge
            files_to_merge /= 2;
            files_to_merge += 1;
        }

        // Merge two files in a multi-threaded approach 
        for (uint32_t i = 0; i < files_to_merge; i++) {
            while (true) {
                bool created_thread = false;
                for (int thread_id = 0; thread_id < cpu_count; thread_id++) {
                    if (thread_status[thread_id].in_progress == false) {
                        thread_status[thread_id].in_progress = true;

                        thread_status[thread_id].filename_a = filename_array[i*2];
                        if (((i*2) + 1) >= num_files) {
                            thread_status[thread_id].filename_b[0] = '\0';
                        } else {
                            thread_status[thread_id].filename_b = filename_array[(i*2) + 1];
                        }

                        sprintf((char*)thread_status[thread_id].output_filename, "%s/merge.%u.gtf.sorted", tmp_dir, i);

                        if(pthread_create((pthread_t*)&thread_status[thread_id].pthread, NULL, merge_thread, (void*)&thread_status[thread_id])) {
                            fprintf(stderr, "Error creating thread\n");
                            return NULL;
                        }

                        created_thread = true;
                        break;
                    }
                }

                if (created_thread) {
                    break;
                } else {
                    // wait for other threads to finish.
                }

                sleep(1);
            }
        }

        // make sure all threads have joined back.
        for (int thread_id = 0; thread_id < cpu_count; thread_id++) {
            pthread_join(thread_status[thread_id].pthread, NULL);
        }

        // Rename merged files
        for (uint32_t i = 0; i < files_to_merge; i++) {
            char old_merge_filename[strlen(tmp_dir) + 64];
            sprintf(old_merge_filename, "%s/merge.%u.gtf.sorted", tmp_dir, i);

            char new_merge_filename[strlen(tmp_dir) + 64];
            sprintf(new_merge_filename, "%s/transcripts.%u.gtf.sorted", tmp_dir, i);
            rename(old_merge_filename, new_merge_filename);
        }
    }

    // final rename
    char old_finished_merged_filename[strlen(tmp_dir) + 64];
    sprintf(old_finished_merged_filename, "%s/transcripts.0.gtf.sorted", tmp_dir);

    char new_finished_merged_filename[strlen(tmp_dir) + 64];
    sprintf(new_finished_merged_filename, "%s/transcripts.sorted.gtf", tmp_dir);
    rename(old_finished_merged_filename, new_finished_merged_filename);

    // Free array of files
    for (uint32_t i = 0; i < num_files; i++) {
        free(filename_array[i]);
    }

    return Py_BuildValue("i", 0);
}

// Mapping between python and c function names
static PyMethodDef CSortMethods[] = {
    {"csort", py_csort, METH_VARARGS},
    {NULL, NULL}
};

// Module initialisation routine
PyMODINIT_FUNC
initcsort(void)
{
    (void) Py_InitModule("csort", CSortMethods);
}
