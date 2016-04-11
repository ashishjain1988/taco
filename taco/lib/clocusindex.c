//
// TACO
// clocusindex.c
//
#include <Python.h>

#define BED_FIELD_SIZE_MAX 256
#define BED_DELIM "\t\r\n"

#define LOCUS_FMT_STRING "L%zu\t%s\t%zu\t%zu\t%" PRIu64 "\t%zu\n"


int
bed_window_overlap(const char *seqname1, size_t start1, size_t end1,
                   const char *seqname2, size_t start2, size_t end2) {
    if (strcmp(seqname1, seqname2) != 0)
        return 0;
    return ((start1 < end2) && (start2 < end1));
}


int
bed_index_loci(char *bed_file, char *output_file) {
    char *line = NULL;
    char *linep;
    size_t linecap = 0;
    ssize_t linelen;

    FILE *infp = fopen(bed_file, "r");
    FILE *outfp = fopen(output_file, "w");
    char chrom[BED_FIELD_SIZE_MAX] = {0};
    size_t start = 0;
    size_t end = 0;
    size_t numlines = 0;
    size_t locus_id = 1;
    uint64_t filepos = 0;
    char *this_chrom;
    size_t this_start = 0;
    size_t this_end = 0;

    // get first line to initialize window
    linelen = getline(&line, &linecap, infp);
    if (linelen == 0) {
        fprintf(stderr, "build_locus_index(): empty input file");
        fclose(infp);
        fclose(outfp);
        return -1;
    }
    // split into fields
    linep = line;
    this_chrom = strsep(&linep, BED_DELIM);
    this_start = atoi(strsep(&linep, BED_DELIM));
    this_end = atoi(strsep(&linep, BED_DELIM));
    // initialize window
    strncpy(chrom, this_chrom, BED_FIELD_SIZE_MAX);
    chrom[sizeof(chrom) - 1] = '\0';
    start = this_start;
    end = this_end;
    numlines = 1;

    while ((linelen = getline(&line, &linecap, infp)) > 0) {
        // split into fields
        linep = line;
        this_chrom = strsep(&linep, BED_DELIM);
        this_start = atoi(strsep(&linep, BED_DELIM));
        this_end = atoi(strsep(&linep, BED_DELIM));

        // check if feature is outside current window
        if (!bed_window_overlap(chrom, start, end, this_chrom, this_start,
                                this_end)) {
            // print current locus
            fprintf(outfp, LOCUS_FMT_STRING, locus_id, chrom, start, end,
                    filepos, numlines);
            // reset window
            strncpy(chrom, this_chrom, BED_FIELD_SIZE_MAX);
            chrom[sizeof(chrom) -1] = '\0';
            start = this_start;
            end = this_end;
            filepos = ((uint64_t) ftello(infp)) - linelen;
            numlines = 1;
            locus_id += 1;
        } else {
            // expand window
            if (this_end > end)
                end = this_end;
            numlines += 1;
        }
    }

    // print final window
    if (numlines > 0) {
        // print current locus
        fprintf(outfp, LOCUS_FMT_STRING, locus_id, chrom, start, end,
                filepos, numlines);
    }

    if (line != NULL)
        free(line);
    fclose(infp);
    fclose(outfp);
    return 0;
}


static PyObject *
py_bed_index_loci(PyObject *self, PyObject *args) {
    char *bed_file;
    char *output_file;
    int ret;

    if (!PyArg_ParseTuple(args, "ss", &bed_file, &output_file)) {
        return NULL;
    }
    ret = bed_index_loci(bed_file, output_file);
    return Py_BuildValue("i", ret);
}


static PyMethodDef CLocusIndexMethods[] = {
    {"bed_index_loci", py_bed_index_loci, METH_VARARGS},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initclocusindex(void)
{
    (void) Py_InitModule("clocusindex", CLocusIndexMethods);
}
