#include <zlib.h>
#include <stdio.h>
#include <float.h>
#include <dirent.h>
#include <string.h>
#include "kseq.h"

// reference: http://lh3lh3.users.sourceforge.net/parsefastq.shtml
// Tao Zhu 2015-08-30

KSEQ_INIT(gzFile, gzread)

// extension
// https://stackoverflow.com/questions/3035225/getting-file-extension-in-c-language
const char *getExt (const char *fspec) {
    char *e = strrchr (fspec, '.');
    if (e == NULL)
        e = ""; // fast method, could also use &(fspec[strlen(fspec)]).
    return e;
}

// main
int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    int l;
    
    struct dirent *pDirent;
    DIR *pDir;
    
    int parsed_files = 0;
    unsigned long int max_size = 0;
    unsigned long int reads = 0;
    unsigned long int yield = 0;
    unsigned int min_size = 65535;
    char fastq[] = ".fastq";
    
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <dir>\n", argv[0]);
        return 1;
    }
    
    pDir = opendir(argv[1]);
    if (pDir == NULL) {
        printf ("Cannot open directory '%s'\n", argv[1]);
        return 1;
    }

    // change pwd
    chdir(argv[1]);
    
    // reads files
    while ((pDirent = readdir(pDir)) != NULL) {
        
        // check if its .fastq or .gz
        if (strcmp(fastq, getExt(pDirent->d_name)) == 0)
        {
            // printf("reading file : %s\n", pDirent->d_name);
            fp = gzopen(pDirent->d_name, "r");
            
            seq = kseq_init(fp);
            while ((l = kseq_read(seq)) >= 0)
            {
                ++reads, yield += seq->seq.l;
                
                // smallest and largest
                if (min_size > seq->seq.l) {
                    min_size = seq->seq.l;
                }
                if (max_size < seq->seq.l) {
                    max_size = seq->seq.l;
                }
            }
            kseq_destroy(seq);
            gzclose(fp);
            
            printf("reads:%ld\n", reads);
            ++parsed_files;
        }
    }

    closedir(pDir);
    
    printf("fastQ files found:%d\nreads:%ld\tbases:%ld\tmax:%ld\tmin:%ld\tavg:%lf\n", parsed_files, reads, yield, max_size, min_size, (float)yield/reads);
    
    return 0;
}
