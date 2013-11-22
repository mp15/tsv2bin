/*
 * Genotype TSV to GTypeX Binary Format converter
 *
 * Copyright (c) Genome Research Limited 2012, 2013.
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Author: Martin Pollard
 *         <mp15@sanger.ac.uk>
 * Created: 2012-09-21
 */

// Add this for getline(), etc
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>

/*
 String lengths are inclusive of terminating NULL
 A) a dictionary of chromosome names chromDict:
 dict len [int] 
 Dictionary magic number "DICT" (char[4])
 Number of chromosomes N (int)
 N x chromosome strings:
 string length len (int)
 chromosome name s (char[len])
 
 B) number of sites n (int)
 
 C) array of chromosomes chrArray (int[n])
 
 D) array of positions posArray (int[n])
 
 E) individual genotype data (read until EOF):
 length of sample name len (int)
 same name (char[len])
 genotype list gList (int[n])
 
 chrArray gives the chromosomes as keys to the dictionary chromDict, and posArray[i] gives the co-ordinates of site i on the chromosome given by chrArray[i]. The entries to gList use the conversion code {0:AA,1:AC,2:AG,3:AT,4:CC,5:CG,6:CT,7:GG,8:GT,9:TT}, with any other number taken to be "NN".
*/

int alleleEncode(const char allele_a, const char allele_b) {
    switch (allele_a)
    {
        case 'A':
            switch (allele_b)
            {
                case 'A':
                    return 0;
                case 'C':
                    return 1;
                case 'G':
                    return 2;
                case 'T':
                    return 3;
                case 'N':
                default:
                    return 10;
            }
        case 'C':
            switch (allele_b)
        {
            case 'A':
                return 1;
            case 'C':
                return 4;
            case 'G':
                return 5;
            case 'T':
                return 6;
            case 'N':
            default:
                return 10;
        }
        case 'G':
            switch (allele_b)
        {
            case 'A':
                return 2;
            case 'C':
                return 5;
            case 'G':
                return 7;
            case 'T':
                return 8;
            case 'N':
            default:
                return 10;
        }
        case 'T':
            switch (allele_b)
        {
            case 'A':
                return 3;
            case 'C':
                return 6;
            case 'G':
                return 8;
            case 'T':
                return 9;
            case 'N':
            default:
                return 10;
        }
        case 'N':
        default:
            return 10;
    }
}

/*** These three structures define the binary format for GLF genotypes ***/

struct genotype_data
{
    int32_t sample_name_len;
    char* sample_name; // [sample_name_len]
    int32_t* gList; // [num_sites]
    struct genotype_data* next; // not in file
};

struct fixed_chrom_string
{
    int32_t chrom_name_len;
    char* chrom_name;
    struct fixed_chrom_string* next;
};

struct binary_format
{
    // In file: const char magic[4];
    int32_t num_chromosomes;
    struct fixed_chrom_string* chrom_strings; // [num_chromosomes]
    int32_t num_sites;
    int32_t* chrArray; // [num_sites]
    int32_t* posArray; // [num_sites]
    int32_t num_samples; // not in file
    struct genotype_data* sample; // [num_samples]
};

/*** This structure holds SNP metadata for processing the file ***/
struct snp_data {
    int32_t chr;
    int32_t pos;
    char* snpid;
    char ref;
    char alt;
    struct snp_data* next;
};

const char* cunning_convert(int input, char ref)
{
    const char* n[] = {"./.","./.","./.","./.","./.","./.","./.","./.","./.","./."};
    const char* a[] = {"0/0","0/1","0/1","0/1","1/1","./.","./.","1/1","./.","1/1"};
    const char* t[] = {"1/1","./.","./.","0/1","1/1","./.","0/1","1/1","0/1","0/0"};
    const char* g[] = {"1/1","./.","0/1","./.","1/1","0/1","./.","0/0","0/1","1/1"};
    const char* c[] = {"1/1","0/1","./.","./.","0/0","0/1","0/1","1/1","./.","1/1"};
    const char** table = n;

    switch (ref) {
        case 'A':
	table = a;
        break;
        case 'T':
        table = t;
        break;
        case 'G':
        table = g;
        break;
        case 'C':
        table = c;
        break;
	default:
	table = n;
    }
    return table[input];
}

/* Writes binary genotype data to specified FILE */
void write_vcf(struct binary_format* bin_data, FILE* output_file, char** snplist, char* snpref, char* snpalt) {
    // Write header
    fprintf(output_file, "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    // Write sample names
    for (int s = 0; s < bin_data->num_samples; ++s) {
        fprintf(output_file,"\t%s", bin_data->sample[s].sample_name);
    }
    // Write end of header
    fprintf(output_file, "\n");
    // Write variant data
    for (int i = 0; i < bin_data->num_sites; ++i) {
        fprintf(output_file, "%s\t%d\t%s\t%c\t%c\t.\t.\t.", bin_data->chrom_strings[bin_data->chrArray[i]].chrom_name, bin_data->posArray[i], snplist[i], snpref[i], snpalt[i]);
        for (int j = 0; j < bin_data->num_samples; ++j) {
            fprintf(output_file, "\t%s", cunning_convert(bin_data->sample[j].gList[i], snpref[i]));
        }
	fprintf(output_file, "\n");
    }
}

/* Writes binary genotype data to specified FILE */
void write_binary_gt(struct binary_format* bin_data, FILE* output_file) {
    // Write header
    int32_t dict_len = strlen("DICT")+1;
    fwrite(&dict_len, sizeof(int32_t), 1, output_file);
    fwrite("DICT", dict_len, 1, output_file);
    fwrite(&bin_data->num_chromosomes, sizeof(bin_data->num_chromosomes), 1, output_file);
    for (int32_t i = 0; i < bin_data->num_chromosomes; ++i) {
        fwrite(&bin_data->chrom_strings[i].chrom_name_len, sizeof(bin_data->chrom_strings[i].chrom_name_len), 1, output_file);
        fwrite(bin_data->chrom_strings[i].chrom_name, bin_data->chrom_strings[i].chrom_name_len, 1, output_file);
    }
    fwrite(&bin_data->num_sites, sizeof(bin_data->num_sites), 1, output_file);
    fwrite(bin_data->chrArray, sizeof(int32_t), bin_data->num_sites, output_file);
    fwrite(bin_data->posArray, sizeof(int32_t), bin_data->num_sites, output_file);
    for (int32_t j = 0; j < bin_data->num_samples; ++j) {
        fwrite(&bin_data->sample[j].sample_name_len, sizeof(int32_t), 1, output_file);
        fwrite(bin_data->sample[j].sample_name, bin_data->sample[j].sample_name_len, 1, output_file);
        fwrite(bin_data->sample[j].gList, sizeof(int32_t), bin_data->num_sites, output_file); 
    }
}

int addOrFindChrom(const char* token, struct fixed_chrom_string** first, int* num_chromosomes)
{
    int arraypos = 0;
    int found = 0;
    struct fixed_chrom_string* prev = NULL;
    struct fixed_chrom_string* chrom = *first;
    while (chrom != NULL)
    {
        if (strcmp(chrom->chrom_name, token) == 0)
        {
            found = 1;
            break;
        }
        prev = chrom;
        chrom = chrom->next;
        arraypos++;
    }
    if (found == 0)
    {
        (*num_chromosomes)++;
        chrom = malloc(sizeof(struct fixed_chrom_string));
        chrom->next = NULL;
        chrom->chrom_name = strdup(token);
        chrom->chrom_name_len = strlen(token)+1;
        // add item to linked list
        if ( prev == NULL ) {*first = chrom;} else { prev->next = chrom; }
    }
    return arraypos;
}

// Does a very crude chomp on the newline
void chomp(char* i)
{
    while (*i != '\0' && *i != '\n') { i++; }
    *i = '\0';
}

bool compareSNPLT(struct snp_data* a, struct snp_data* b)
{
    if (a->chr > b->chr)
    {
        return true;
    }
    else
    {
        if ((a->chr == b->chr) && (a->pos > b->pos))
            return true;
        else
            return false;
    }
}

// inserts a snp_data item into the linked list ensuring the list stays sorted
void insertSNPSorted(struct snp_data** snplist, struct snp_data* insert)
{
    // First item to be inserted or less than first item?
    if (*snplist == NULL || compareSNPLT(*snplist, insert))
    {
        insert->next = *snplist;
        *snplist = insert;
        return;
    }
    
    struct snp_data* prev = *snplist;
    struct snp_data* curr = (*snplist)->next;
    while (curr != NULL)
    {
        if (compareSNPLT(curr, insert))
        {
            prev->next = insert;
            insert->next = curr;
            return;
        }
        // next
        prev = curr;
        curr = prev->next;
    }
    prev->next = insert;
    return;
}

/*
 * Parses a map file containing the coordinates of the SNPs used in our input file.
 * Information from the map is then used to populate the chromosome dictionary and
 * header information in the binary structure.
 */
int parse_map(struct binary_format* bin_data, FILE* map_file, char*** snplist, char**reflist, char**altlist) {
    char* line = NULL;
    size_t length = 0;
    struct fixed_chrom_string* chrom_list = NULL;
    struct snp_data* first = NULL;
    int snpcount = 0;
    
    bin_data->num_chromosomes = 0;
    
    while (feof(map_file) == 0 && ferror(map_file) == 0) {
        if (getline(&line, &length, map_file) == -1) break;
        
        chomp( line );
        
        // Is it a comment line?
        if ( line[0] == '#' ) { continue; }
        
        int rowpos = 0;
        char* field = line;
        char* token = strsep(&field, "\t");
        int chr = -1;
        int pos = -1;
        char* snpid = NULL;
        char ref = '\0';
        char alt = '\0';
        while (token != NULL && rowpos < 4) {
            switch (rowpos) {
                case 0: // chromosome
                    chr = addOrFindChrom(token, &chrom_list, &bin_data->num_chromosomes);
                    break;
                case 1: // rsid
                    snpid = strdup(token);
                    break;
                case 2:
                    break;
                case 3: // coordinates within chromosome
                    pos = atoi(token);
                    break;
                case 4:
                    ref = token[0];
                    break;
                case 5:
                    alt = token[1];
                    break;
            }
            ++rowpos;
            token = strsep(&field, "\t");
        }
        if (chr != -1 && pos != -1 && snpid != NULL)
        {
            struct snp_data* snp = (struct snp_data*)malloc(sizeof(struct snp_data));
            snp->snpid = snpid;
            snp->pos = pos;
            snp->chr = chr;
            snp->ref = ref;
            snp->alt = alt;
            snp->next = NULL;
            insertSNPSorted(&first, snp);
            snpcount++;
            //if (prev == NULL) { first = prev = snp; } else { prev->next = snp; prev = snp; }
        }
    }
    free(line);

    // Turn linked list to array
    struct fixed_chrom_string* iter = chrom_list;
    bin_data->chrom_strings = calloc(bin_data->num_chromosomes, sizeof(struct fixed_chrom_string));
    int k = 0;
    while (iter != NULL) {
        memcpy( &bin_data->chrom_strings[k], iter, sizeof(struct fixed_chrom_string));
        ++k;
        struct fixed_chrom_string *old = iter;
        iter = iter->next;
        free(old);
    }

    if (ferror(map_file) != 0) {
        printf("Error reading map file\r\n");
        return -1;
    }
    
    // Allocate the sites data based on map
    bin_data->num_sites = snpcount;
    bin_data->chrArray = (int32_t*)calloc(bin_data->num_sites, sizeof(int32_t));
    bin_data->posArray = (int32_t*)calloc(bin_data->num_sites, sizeof(int32_t));
    char** snpNames = (char**)calloc(bin_data->num_sites+1, sizeof(char*)); // +1 so we can use a NULL pointer as terminator
    char* snpRef = (char*)calloc(bin_data->num_sites, sizeof(char));
    char* snpAlt = (char*)calloc(bin_data->num_sites, sizeof(char));
    struct snp_data* curr = first;
    for (int i = 0; i < snpcount; ++i) {
        bin_data->chrArray[i] = curr->chr;
        bin_data->posArray[i] = curr->pos;
        snpNames[i] = curr->snpid;
        snpRef[i] = curr->ref;
        snpAlt[i] = curr->alt;
        struct snp_data* prev = curr;
        curr = curr->next;
        free(prev);
    }
    
    *snplist = snpNames;
    *reflist = snpRef;
    *altlist = snpAlt;

    return 0;
}

int getChrPosArray(char** snplist, const char* snpid){
    char** iter = snplist;
    int counter = 0;
    while (*iter != NULL)
    {
        if (strcmp(snpid, *iter) == 0)
            return counter;
        iter++;
        counter++;
    }
    abort();
}

int parse_input(struct binary_format* bin_data, FILE* input_file, char** snplist) {
    // Parse the header
    char* line = NULL;
    size_t length = 0;
    getline(&line, &length, input_file);
    
    chomp(line);
    
    // Count tags
    int field_count = 1;
    char* lp = line;
    char* fp = line;
    while(*lp != '\0') {
        if (*lp == '\t') field_count++;
        lp++;
        if (field_count == 8) {fp = lp;} // Save the position of the first SNP txt
    }
    if (field_count < 8 || *fp == '\0') {
        printf("Not enough fields in header.\r\n");
        return -1;
    }
    fp++; // Move on from the tab

    // Create indirection list
    int* inDirection = (int32_t*)calloc(bin_data->num_sites, sizeof(int32_t));
    char* snpid = strsep(&fp, "\t");
    for (int i = 0; i < bin_data->num_sites; ++i) {
        if (snpid == NULL)
        {
            printf("invalid input file header\r\n");
            return -1;
        }
        inDirection[i] = getChrPosArray(snplist, snpid);
        snpid = strsep(&fp, "\t");
    }
    
    // Now parse the body
    struct genotype_data* first = NULL;
    struct genotype_data* prev = NULL;
    bin_data->num_samples = 0;

    while (feof(input_file) == 0 && ferror(input_file) == 0) {
        if (getline(&line, &length, input_file) == -1) break;
        // Do a crude chomp on the newline
        chomp(line);
        // Is it a comment line?
        if ( line[0] != 'S' ) { continue; }
           
        int rowpos = 0;
        char* field = line;
        char* token = strsep(&field, "\t");
        char* sangerid = NULL;
        while (token != NULL && rowpos < 8) {
            switch (rowpos) {
                case 0: // rowtype
                    break;
                case 2: // sanger id
                    sangerid = token;
                    break;
                    // GNDN fields
                case 1: // extname
                case 3: // plexchip
                case 4: // wellid
                case 5: // plate
                case 6: // plate name
                case 7: // well
                    break;
                default:
                    break;
            }
            ++rowpos;
            token = strsep(&field, "\t");
        }
        if (rowpos != 8 || sangerid == NULL) { continue; }
        
        // Now create sample record
        struct genotype_data* sample = (struct genotype_data*)malloc(sizeof(struct genotype_data));
        sample->sample_name_len = strlen(sangerid)+1;
        sample->sample_name = strdup(sangerid);
        sample->gList = calloc(bin_data->num_sites, sizeof(int32_t));
        sample->next = NULL;
        rowpos = 0;
        // finished parsing the fixed data, now we're on to the SNP genotypes
        while (token != NULL && rowpos < bin_data->num_sites) {
            sample->gList[inDirection[rowpos]] = alleleEncode(token[0], token[1]);
            ++rowpos;
            token = strsep(&field, "\t");
        }
        // add item to linked list
        if ( prev == NULL ) {first = prev = sample;} else { prev->next = sample; prev = sample; }
        bin_data->num_samples++;
    }
    free(line);
    free(inDirection);
    // Turn linked list to array
    struct genotype_data* iter = first;
    bin_data->sample = calloc(bin_data->num_samples, sizeof(struct genotype_data));
    int k = 0;
    while (iter != NULL) {
        memcpy( &bin_data->sample[k], iter, sizeof(struct genotype_data));
        ++k;
        struct genotype_data *old = iter;
        iter = iter->next;
        free(old);
    }
    
    if (ferror(input_file) != 0) {
        printf("Error reading input file\r\n");
        return -1;
    }
    return 0;
}


int main (int argc, char *argv[])
{
    // parse inputs this should probably be in subroutine
    // Options
    char* file_to_parse = NULL;
    char* file_read = NULL;
    char* file_map = NULL;
    char* output_type = NULL;
    
    // parse inputs
    const char* opts = "i:o:m:t:";
    const struct option longopts[] = {{NULL, 0, NULL, 0},};
    int option_index = 0;

    int c;
    while ((c = getopt_long(argc, argv, opts, longopts, &option_index)) != -1) {
        switch (c) {
            case 'i':
                file_to_parse = optarg;
                break;
            case 'm':
                file_map = optarg;
                break;
            case 'o':
                file_read = optarg;
                break;
            case 't':
                output_type = optarg;
                break;
            default:
                fprintf(stderr, "\r\n");
        }
    }
    
    // validate input
    if (file_to_parse == NULL) {
        fprintf(stderr, "You must specify a file to parse\r\n");
        return -1;
    }
    if (file_map == NULL) {
        fprintf(stderr,"You must specify a map file.\r\n");
        return -1;
    }
    
    // Now open file
    FILE* input_file = fopen(file_to_parse, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Cannot open input file for reading: %s\r\n", file_to_parse);
        return -1;
    }

    // Now open file
    FILE* map_file = fopen(file_map, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Cannot open map file for reading: %s\r\n", file_to_parse);
        return -1;
    }

    // Then open output file
    if (file_read == NULL) { file_read = "default.gtypex"; }
    FILE* output_file = fopen(file_read, "w");
    if (output_file == NULL) {
        fprintf(stderr, "Cannot open output file for writing: %s\r\n", file_read);
        return -1;
    }

    // convert input data to output data
    struct binary_format* bin_data = (struct binary_format*) malloc(sizeof(struct binary_format));

    // Parse the map file
    char** snplist;
    char* snpref;
    char* snpalt;
    if (parse_map(bin_data, map_file, &snplist, &snpref, &snpalt) != 0)
    {
        return -1;
    }
    fclose(map_file);

    // parse the gt input file
    if (parse_input(bin_data, input_file, snplist) != 0)
    {
        return -1;
    }
    fclose(input_file);

    // Write records to output file
    if (!strcmp(output_type,"vcf")) {
        write_vcf(bin_data, output_file, snplist, snpref, snpalt);
    } else {
        write_binary_gt(bin_data, output_file);
    }
    
    // Clean up for valgrind
    //for (char **iter = snplist; iter != NULL; )
    //{
        //if (*iter != NULL) free(*iter);
    //    iter++;
    //}
    free (snplist);
    for (int x = 0; x< bin_data->num_chromosomes; x++) {
        free(bin_data->chrom_strings[x].chrom_name );
    }
    free(bin_data->chrom_strings);
    for (int y = 0; y< bin_data->num_samples; y++) {
        free(bin_data->sample[y].gList);
        free(bin_data->sample[y].sample_name);
    }
    free(bin_data->chrArray);
    free(bin_data->posArray);
    free(bin_data->sample);
    free(bin_data);
    
    fclose(output_file);

    return(0);
}

