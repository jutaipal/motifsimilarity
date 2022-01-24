//
//  motifsimilarity.c
//  motifsimilarity calculates similarity between two PWMs
//  by comparing maximum scores for all possible kmers
//  using cosine similarity
//  https://en.wikipedia.org/wiki/Cosine_similarity
//  Created by Jussi Taipale on 22/01/2022.
//  for parallel, compile: clang -latomic -fopenmp -lm -O3 -o motif main.c
//  for one core: clang -lm -O3 -o motif main.c
//  use: ./motifsimilarity a.pfm b.pfm gapped 12
//  for gapped 12-mer, should take approx 30 seconds.
//  PWM motifs must be formatted using tab separated columns,
//  with four rows in alphabetical order (A, C, G, T)


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/* GLOBAL VARIABLES */
double pseudocount = 0.0001;
short int max_width_of_pwm = 84; /* Nlength * 2 + 2 */
short int contacts = 0;
char *VERSION = "spacek40 v0.200 DEC 19 2018";
short int Nlength = 20;
short int rna = 0;

/* STRUCTURES */

/* ORIENTED MATCH */
struct oriented_match {short int position; signed short int strand; double score; short int id;};
short int oriented_match_init (struct oriented_match *i)
{
    (*i).position = 0;
    (*i).strand = 0;
    (*i).score = 0;
    (*i).id = 0;
    return(0);
}

/* COUNT PWM */
struct count_pwm {char *name; short int width; long int max_counts; double **incidence;};

short int count_pwm_clear (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
for (counter = 0; counter < 5 + contacts * 12; counter++)
{
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}

short int count_pwm_init (struct count_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).name = malloc(1000);
strcpy ((*i).name, name);
(*i).width = width;
(*i).max_counts = initial_value;
(*i).incidence = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
for (counter = 0; counter < 5 + contacts * 12; counter++)
{
(*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).incidence[counter][counter2] = initial_value;
}
return(0);
}

short int count_pwm_free (struct count_pwm *i)
{
    short int counter;
    free((*i).name);
    for (counter = 0; counter < 5; counter++) free((*i).incidence[counter]);
    free((*i).incidence);
    return(0);
}


struct pairwise_correlation {short int first_base; short int second_base; double delta_ic; short int max_dinucleotide; short int min_dinucleotide; double min_fold_change; double max_fold_change;};

/* NORMALIZED PWM */
struct normalized_pwm {char *name; char *seed; short int width; long int max_counts; double *information_content; short int *original_position; double *position_score; long int *total_counts_for_column; double **fraction; struct pairwise_correlation *pairwise_correlation; short int negative_values_allowed; struct oriented_match match;};
short int normalized_pwm_init (struct normalized_pwm *i, char *name, short int width, double initial_value)
{
short int maximum_width = max_width_of_pwm;
short int counter;
short int counter2;
(*i).negative_values_allowed = 0;
(*i).pairwise_correlation = malloc(sizeof(struct pairwise_correlation) * 10 + 5);
for (counter = 0; counter < 10; counter++)
{
(*i).pairwise_correlation[counter].first_base = 0;
(*i).pairwise_correlation[counter].delta_ic = 0;
(*i).pairwise_correlation[counter].second_base = 0;
(*i).pairwise_correlation[counter].max_dinucleotide = 0;
(*i).pairwise_correlation[counter].min_dinucleotide = 0;
(*i).pairwise_correlation[counter].max_fold_change = 0;
(*i).pairwise_correlation[counter].min_fold_change = 0;
}
(*i).name = malloc(100);
strcpy ((*i).name, name);
(*i).seed = malloc(1000);
strcpy ((*i).seed, "UNKNOWN");
(*i).width = width;
(*i).max_counts = initial_value;
(*i).fraction = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
(*i).information_content = malloc(sizeof(double) * maximum_width + 5);
(*i).position_score = malloc(sizeof(double) * maximum_width + 5);
(*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
(*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

for (counter = 0; counter < 5 + contacts * 12; counter++)
{
(*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
for (counter2 = 0; counter2 < maximum_width; counter2++) (*i).fraction[counter][counter2] = initial_value;
}
for (counter2 = 0; counter2 < maximum_width; counter2++)
{
(*i).information_content[counter2] = 0;
(*i).position_score[counter2] = 0;
(*i).original_position[counter2] = counter2;
(*i).total_counts_for_column[counter2] = 0;
}
(*i).match.position = 0;
(*i).match.strand = 0;
(*i).match.score = 0;
return(0);
}
short int normalized_pwm_free (struct normalized_pwm *i)
{
short int counter;
free((*i).name);
free((*i).pairwise_correlation);
free((*i).information_content);
free((*i).position_score);
free((*i).total_counts_for_column);
for (counter = 0; counter < 5; counter++) free((*i).fraction[counter]);
free((*i).fraction);
return(0);
}


/* SUBROUTINES */

/* SUBROUTINE THAT RENORMALIZES NORMALIZED PWM (ROWS IN EACH COLUMN ADD TO 1) */
short int Normalize_pwm (struct normalized_pwm *n)
{
    short int counter;
    short int position;
    double total_nucleotides = 0;
    double normalized_value = 0;
    for (position = 0; position < (*n).width; position++)
    {
        for (counter = 0, total_nucleotides = 0; counter < 4; counter++)
        {
        if ((*n).fraction[counter][position] > 0) total_nucleotides += (*n).fraction[counter][position];
        else if ((*n).negative_values_allowed == 1) total_nucleotides += -(*n).fraction[counter][position];
        }
        
        for (counter = 0; counter < 4; counter++)
        {
        normalized_value = ((double) (*n).fraction[counter][position]) / total_nucleotides;
        if ((normalized_value < 0) && ((*n).negative_values_allowed == 0)) normalized_value = 0;
        (*n).fraction[counter][position] = normalized_value;
        }
    }
    return (0);
}

/* SUBROUTINE THAT LOADS A PWM AND NORMALIZES IT */
short int Load_pwm (struct normalized_pwm *p, char *filename, short int normalize)
{
long int counter;
char text1 = '\0';
short int line = 0;
short int pwm_position = 0;
char *current_string;
current_string = malloc(200);
FILE *pwmfile;
if ((pwmfile = fopen(filename, "r")) == (void *)0) {printf("\nNo File: %s", filename); exit (2);}
for(line = 0; line <= 3 + contacts * 12;)
{
    for(counter = 0; counter < 30; counter++)
    {
        text1 = getc(pwmfile);
        if (text1 == EOF || text1 == '\n' || text1 == '\t')
        {
        current_string[counter] = '\0';
        if (counter > 0 && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3' || current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7' || current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-'))
        {
        (*p).fraction[line][pwm_position] = atof(current_string);
        /* printf("\n%f", (*p).fraction[line][pwm_position]); */
        pwm_position++;
        }
        if (text1 == '\n' || text1 == EOF) {(*p).width = pwm_position; line++; pwm_position = 0;}
        break;
        }
        current_string[counter]=text1;
        /* printf ("%c", text1); */
    }
}
free (current_string);
if (normalize == 1) Normalize_pwm(p);
// for (line = 0; line < 4 + contacts * 12; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
if (text1 == EOF && line != 3) return(1);
else return (0);
}

/* SUBROUTINE THAT CONVERTS A PWM TO A LOG RATIO PWM (DEFINES SPECIFICITY OF EACH NUCLEOTIDE AT EACH POSITION) */
short int Log_ratio_pwm (struct normalized_pwm *p)
{
short int line = 0;
short int pwm_position = 0;
for (pwm_position = 0; pwm_position < (*p).width; pwm_position++)
{
for(line = 0; line <= 3; line++)
{
(*p).fraction[line][pwm_position] = log10(((*p).fraction[line][pwm_position] + pseudocount) / (pseudocount + 1 - (*p).fraction[line][pwm_position]));
}
}
for (line = 0; line < 4; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
return(0);
}

/* SUBROUTINE THAT ADDS ZERO INFORMATION CONTENT (EQUAL) FLANKS TO A PWM */
short int Add_flanks(struct normalized_pwm *p, short int length, double fill_value)
{
short int nucleotide;
short int pwm_position = (*p).width + 2 * length;
for ( ; pwm_position >= (*p).width + length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = fill_value;
for ( ; pwm_position >= length; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = (*p).fraction[nucleotide][pwm_position-length];
for ( ; pwm_position >= 0; pwm_position--) for (nucleotide = 0; nucleotide < 4; nucleotide++) (*p).fraction[nucleotide][pwm_position] = fill_value;
(*p).width = (*p).width + 2 * length;
return(0);
}

/* faster SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR KMER IN A WIDER PWM  */
double fastKmerscore (double **pwm, signed short int pwm_width, long long int kmer_sequence_value, signed short int kmer_length)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    short int counter;
    signed short int pwm_position;
    signed short int current_position;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    
    for (counter = 0; counter < kmer_length; counter++)
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
    }
    /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length);
     printf("\nREVERSE_SEQ: "); Kmerprint(reverse_complement_kmer_sequence_value, kmer_length); */
    
    for (pwm_position = 0; pwm_position < pwm_width - kmer_length; pwm_position++)
    {
        
        current_position = pwm_width-1-pwm_position;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value; kmer_position < kmer_length; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2)
        {
            nucleotide = (current_kmer_sequence_value & 3);
            //           nucleotide = ((nucleotide & 3) << 2) + ((nucleotide & 12) >> 2);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            //         rev_nucleotide = ((rev_nucleotide & 3) << 2) + ((rev_nucleotide & 12) >> 2);
            // printf("\nDinucleotide %c%c at kmer_position %i and pwm_position %i, score addition %.2f", forward[(nucleotide & 12) >> 2], forward[nucleotide & 3], kmer_position, pwm_position, (*a).fraction[nucleotide][(*a).width-1-pwm_position-kmer_position]);
            score += pwm[nucleotide][current_position-kmer_position];
            revscore += pwm[rev_nucleotide][current_position-kmer_position];
        }
        if (revscore > score && rna == 0) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            best_score = score;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
    
    return (best_score);
}

/* SUBROUTINE THAT FINDS BEST SCORE AND POSITION FOR A GAPPED KMER IN A WIDER PWM  */
double fastgappedKmerscore (double **pwm, signed short int pwm_width, long long int kmer_sequence_value, long long int kmer2_sequence_value, signed short int kmer_length, signed short int gap)
{
    signed short int nucleotide;
    signed short int rev_nucleotide;
    signed short int kmer_position;
    double score = 0;
    double revscore = 0;
    double best_score = -100;
    short int counter;
    signed short int pwm_position;
    signed short int current_position;
    long long int reverse_complement_kmer_sequence_value = 0;
    long long int current_reverse_complement_kmer_sequence_value;
    long long int current_kmer_sequence_value = kmer_sequence_value;
    long long int reverse_complement_kmer2_sequence_value = 0;
    long long int current_reverse_complement_kmer2_sequence_value;
    long long int current_kmer2_sequence_value = kmer2_sequence_value;
    
    /* GENERATES REVERSE COMPLEMENTS */
    for (counter = 0; counter < kmer_length; counter++)
    {
        reverse_complement_kmer_sequence_value <<= 2;
        reverse_complement_kmer2_sequence_value <<= 2;
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer2_sequence_value & 3));
        reverse_complement_kmer2_sequence_value |= (3 ^ (current_kmer_sequence_value & 3));
        current_kmer_sequence_value >>= 2;
        current_kmer2_sequence_value >>= 2;
    }
    
    for (pwm_position = 0; pwm_position < pwm_width - kmer_length * 2 - gap; pwm_position++)
    {
        /* printf("\nFORWARD_SEQ: "); Kmerprint(kmer_sequence_value, kmer_length); */
        
        current_position = pwm_width-1-pwm_position;
        for(kmer_position = 0, score= 0, revscore = 0, current_kmer_sequence_value = kmer_sequence_value, current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value, current_kmer2_sequence_value = kmer2_sequence_value, current_reverse_complement_kmer2_sequence_value = reverse_complement_kmer2_sequence_value; kmer_position < kmer_length; kmer_position++, current_kmer_sequence_value >>= 2, current_reverse_complement_kmer_sequence_value >>= 2, current_kmer2_sequence_value >>= 2, current_reverse_complement_kmer2_sequence_value >>= 2)
        {
            /* KMER1 */
            nucleotide = (current_kmer_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer_sequence_value & 3);
            score += pwm[nucleotide][current_position-kmer_position];
            revscore += pwm[rev_nucleotide][current_position-kmer_position];
            /* KMER2 */
            nucleotide = (current_kmer2_sequence_value & 3);
            rev_nucleotide = (current_reverse_complement_kmer2_sequence_value & 3);
            score += pwm[nucleotide][current_position-kmer_position-gap-kmer_length];
            revscore += pwm[rev_nucleotide][current_position-kmer_position-gap-kmer_length];
        }
        if (revscore > score && rna == 0) score = revscore;
        /* printf ("\nScore %f at position %i", score, (*p).width-1-pwm_position-kmer_position); */
        if (score >= best_score) /* SCORES PWM AT THIS POSITION */
        {
            best_score = score;
            /* printf ("\nMATCH at %i: score %f higher than cut-off %f", match.position[number_of_matches], score, cut_off); */
        }
        /* else printf ("\nNo match at %i: %i is not %i", Nlength - current_position - query_sequence_length, query_sequence_ULL, test_sequence_ULL & query_sequence_ULL); */
    }
    
    return (best_score);
}

/* MAIN PROGRAM */

int main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");
    
    long double t0 = time((void *)0);
    long double t1;
    
    char *searchstring;
    searchstring = malloc(1000);
    strcpy(searchstring, "-");
    
    char *tempstring;
    tempstring = malloc(1000);
    tempstring[0] = '\0';
    
    long double sum_x_squared = 0;
    long double sum_y_squared = 0;
    long double sum_xy = 0;
    
    short int firstnoncommandposition = 0;
    
    short int gapped = 0;
    short int gappedkmerlength = 4;
    short int too_long_kmer = 8;
    short int gap;
    long double onetwo;
    long double twoone;
    long double gapped_correlation;
    long double firstscore;
    long double secondscore;
    double **firstpwm;
    double **secondpwm;
    signed short int firstlength;
    signed short int secondlength;
    signed short int secondpwmwidth;
    signed short int firstpwmwidth;
    signed short int max_pwm_width;
    long long int kmervalue;
    long long int kmervalue2;

    /* PWM STRUCTURE */
    struct normalized_pwm qp;
    normalized_pwm_init(&qp, "empty", Nlength * 2, 0);
    
    /* PWM STRUCTURE2 */
    struct normalized_pwm mfp;
    normalized_pwm_init(&mfp, "empty", Nlength * 2, 0);
    
    
    strcpy(searchstring, argv[1 + firstnoncommandposition]);
    printf("\n\nNORMALIZED MATRIX 1:\t%s", searchstring);
    Load_pwm (&qp, searchstring, 0);
    firstlength = qp.width;
    Normalize_pwm(&qp);
    printf("\n\nLOG MATRIX 1:\t%s", searchstring);
    Log_ratio_pwm(&qp);
    
    strcpy(tempstring, argv[2 + firstnoncommandposition]);
    printf("\n\nNORMALIZED MATRIX 2:\t%s", tempstring);
    Load_pwm (&mfp, tempstring, 0);
    secondlength = mfp.width;
    Normalize_pwm(&mfp);
    printf("\n\nLOG MATRIX 2:\t%s", tempstring);
    Log_ratio_pwm(&mfp);

    if(firstlength > secondlength) too_long_kmer = firstlength;
    else too_long_kmer = secondlength;
    max_pwm_width = too_long_kmer;
    
    if (argc >= 4 + firstnoncommandposition)
    {
    if (strcmp(argv[3 + firstnoncommandposition], "gapped") == 0) gapped = 1;
    else firstlength = atoi(argv[3 + firstnoncommandposition]);
    }
    if (argc >= 5 + firstnoncommandposition)
    {
    secondlength = atoi(argv[4 + firstnoncommandposition]);
    secondlength = secondlength - secondlength % 2;
    gappedkmerlength = secondlength / 2;
    }
    if (gapped == 1)
    {
    firstlength = secondlength / 2;
    if (max_pwm_width < gappedkmerlength * 2) gappedkmerlength = max_pwm_width / 2;
    }

    if(firstlength > secondlength) too_long_kmer = firstlength;
    else too_long_kmer = secondlength;
    
    Add_flanks(&qp, too_long_kmer-2, -0.60206);
    Add_flanks(&mfp, too_long_kmer-2, -0.60206);

    firstpwm = qp.fraction;
    firstpwmwidth = qp.width;
    
    secondpwm = mfp.fraction;
    secondpwmwidth = mfp.width;
    
    long long int last_kmer_value;
    
    sum_x_squared = 0;
    sum_y_squared = 0;
    sum_xy = 0;
    
    last_kmer_value = pow(4,firstlength);
    #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore)
    for(kmervalue = last_kmer_value; kmervalue > 0; kmervalue--)
    {
    firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, firstlength);
    secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, firstlength);
    sum_x_squared += pow(10,firstscore+firstscore);
    sum_y_squared += pow(10,secondscore+secondscore);
    sum_xy += pow(10,firstscore+secondscore);
   /*     if(firstscore != secondscore) different++;
        else same++; */
    }
    onetwo = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));

    /* printf("\ndifferent %li, same %li\n", different, same); */
           
    if(firstlength != secondlength)
    {
    sum_x_squared = 0;
    sum_y_squared = 0;
    sum_xy = 0;
    
    last_kmer_value = pow(4,secondlength);
    #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore)
    for(kmervalue = 0; kmervalue < last_kmer_value; kmervalue++)
    {
    firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, secondlength);
    secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, secondlength);
    sum_x_squared += pow(10,firstscore+firstscore);
    sum_y_squared += pow(10,secondscore+secondscore);
    sum_xy += pow(10,firstscore+secondscore);
    }
    twoone = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));
    }
    else twoone = onetwo;
    

    /* GAPPED kmers */
    if (gapped == 1)
    {
    
    sum_x_squared = 0;
    sum_y_squared = 0;
    sum_xy = 0;
    
    last_kmer_value = pow(4,gappedkmerlength);
    
#pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,kmervalue2, gap, firstscore,secondscore)
    for(gap = max_pwm_width - 2 * gappedkmerlength; gap >= 0; gap--)
    {
    for(kmervalue = last_kmer_value; kmervalue > 0; kmervalue--)
    {
    for(kmervalue2 = last_kmer_value; kmervalue2 > 0; kmervalue2--)
    {

        firstscore = fastgappedKmerscore(firstpwm, firstpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap);
        secondscore = fastgappedKmerscore(secondpwm, secondpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap);
        sum_x_squared += pow(10,firstscore+firstscore);
        sum_y_squared += pow(10,secondscore+secondscore);
        sum_xy += pow(10,firstscore+secondscore);
    }
    }
    }
    gapped_correlation = sum_xy/(sqrt(sum_y_squared)*sqrt(sum_x_squared));
        
    printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n%s to %s (gapped %i-mer)\t%Lg\n", searchstring, tempstring, firstlength, onetwo, tempstring, searchstring, secondlength, twoone, searchstring, tempstring, gappedkmerlength * 2, gapped_correlation);
    }
    else
    {
    printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n", searchstring, tempstring, firstlength, onetwo, tempstring, searchstring, secondlength, twoone);
    }
    
    t1 = time((void *)0);
    printf ("\nTime: %ld seconds\n", (long) (t1 - t0));
    exit(0);

    return 0;
}
