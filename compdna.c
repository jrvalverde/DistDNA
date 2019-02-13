/**
 *
 * compdna - Composition vector analysis for DNA
 *
 *	Calculate composition vector for DNA sequences
 *	TEST VERSION
 *
 *	Usage: compdna options
 *			-i infile
 *			-o outfile
 *			-k k-word size
 *			-a		(toggle on ambiguity)
 *
 * (C) Jose R. Valverde, EMBnet/CNB, CSIC. 2009
 *	jrvalverde@cnb.csic.es
 *
 * License: GNU GPL
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * @package 	distdna
 * @author  	Jose R Valverde <jrvalverde@cnb.csic.es>
 * @copyright	EMBnet/CNB
 * @license 	c/gpl.txt
 *
 * $Id: compdna.c,v 1.1 2009/10/05 13:32:21 root Exp root $
 *
 * $Log: compdna.c,v $
 * Revision 1.1  2009/10/05 13:32:21  root
 * Initial revision
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <getopt.h>

#define DEBUG	1

/*
 * Alphabet definition: 
 * 
 * all valid nucleotide codes as numbers
 */
enum {
    A,
    T,
    G,
    C,
    R,				/* Purine,          A or G */
    Y,				/* Pyrimidine       T or C */
    M,				/* Amino            A or C */
    K,				/* Keto             G or T */
    S,				/* Strong (3bond)   C or G */
    W,				/* Weak (2 bond)    A or T */
    B,				/* not A            T, C, G */
    D,				/* not C            A, T, G */
    H,				/* not G            A, T, C */
    V,				/* not T            A, C, G */
    N				/* Any              A, T, C, G */
} iupac_codes;

/* Valid nucleotide codes as characters */
char *nt = "ATGCRYMKSWBDHVN";

/* relative frequencies of contributing residues in each code */
/*	or.. information payload of each code */
double relfreq[] = {
    1.,				/* A */
    1.,				/* T */
    1.,				/* C */
    1.,				/* G */
    0.5,			/* R */
    0.5,			/* Y */
    0.5,			/* M */
    0.5,			/* K */
    0.5,			/* S */
    0.5,			/* W */
    1. / 3.,			/* B */
    1. / 3.,			/* D */
    1. / 3.,			/* H */
    1. / 3.,			/* V */
    0.25			/* N */
};

/*
 * forward definitions
 */
void usage();
void error(char *msg);
void read_sequence(char *infile, char **seq, int *slen);
long int ipow(int b, unsigned int exp);
void compute_freqs012(char *s, int slen,
		double *f, double *f1, double *f2, int ksize);
void compute_kstring_frequencies(char *s, int ssize, double *fk, 
		double *fk_1, double *fk_2, int k);
void
add_hyperkstr(double *h, int k, double *t, double prob, int pos);
void print_kfreqs(char *outfile, double *fcv, int k);
void compute_corrected_CV(double *a, 
			  double *f, double *f_1, double *f_2, 
			  int ksize);
double compute_k_distance(double *a, double *b, int k);



/*
 * CompDNA
 *
 *	Compute DNA composition in k-words (frequency
 * of each k-word on a DNA strand)
 */
int main(int argc, char *argv[])
{
    int opt, i, ambig, correct;
    extern char *optarg;
    extern int optind, optopt;
    int ksize, slen, rlen;
    char *infile, *outfile, *s, *reference, *r;
    double *fcv, *fcv_1, *fcv_2, *a, *b, d;

    /* initialize values */
    opt = ambig = 0;
    ksize = slen = 0;
    infile = reference = outfile = s = NULL;
    fcv = fcv_1 = fcv_2 = NULL;
    d = 0.0;

    /* parse arguments */
    while ((opt = getopt(argc, argv, "i:o:k:ac1:2:")) != -1) {
    	switch (opt ) {
	    case '1':
	    case 'i': infile = optarg; break;
    	    case '2': reference = optarg; break;
	    case 'o': outfile = optarg; break;
	    case 'k': ksize = atoi(optarg); break;
	    case 'a': ambig = 1; break;
	    case 'c': correct = 1; break;
	    case ':': fprintf(stderr, "Option -%c requires a filename\n",
	    	      optopt);
		      exit(1);
		      break;
	    default: usage();
		     exit(1);
		     break;
	}
    }
    if (ksize < 3) {
    	fprintf(stderr, "K-size must be >= 3\n");
	exit(1);
    }


    /* open sequence */
    read_sequence(infile, &s, &slen);

    if (reference == NULL) {
    	/* simply output sequence composition */
	/* allocate frequency vectors */
	fcv = calloc(sizeof(double), ipow(4, ksize));
	fcv_1 = calloc(sizeof(double), ipow(4, (ksize - 1)));
	fcv_2 = calloc(sizeof(double), ipow(4, (ksize - 2)));
	if ((fcv == NULL) || (fcv_1 == NULL) || (fcv_2 == NULL)) {
	    error("Not enough memory");
	}

	/* compute frequencies */
	if (ambig)
    	    compute_kstring_frequencies(s, slen, fcv, fcv_1, fcv_2, ksize);
	else
    	    compute_freqs012(s, slen, fcv, fcv_1, fcv_2, ksize);

	/* print out results for checking */
	print_kfreqs(outfile, fcv, ksize);
	print_kfreqs(outfile, fcv_1, ksize-1);
	print_kfreqs(outfile, fcv_2, ksize-2);

	a = calloc(sizeof(double), ipow(4, ksize));
	compute_corrected_CV(a, fcv, fcv_1, fcv_2, ksize); 
	print_kfreqs(outfile, a, ksize);
    }
    else {
    	/* compare two sequences (compute distance matrix) */
	read_sequence(reference, &r, &rlen);
	/* allocate frequency vectors */
	/* we'll reuse them for both sequences */
	fcv = calloc(sizeof(double), ipow(4, ksize));
	fcv_1 = calloc(sizeof(double), ipow(4, (ksize - 1)));
	fcv_2 = calloc(sizeof(double), ipow(4, (ksize - 2)));
	if ((fcv == NULL) || (fcv_1 == NULL) || (fcv_2 == NULL)) {
	    error("Not enough memory");
	}

	/* compute frequencies */
	if (ambig)
    	    compute_kstring_frequencies(s, slen, fcv, fcv_1, fcv_2, ksize);
	else
    	    compute_freqs012(s, slen, fcv, fcv_1, fcv_2, ksize);

    	/* compute corrected vector */
	a = calloc(sizeof(double), ipow(4, ksize));
	compute_corrected_CV(a, fcv, fcv_1, fcv_2, ksize); 

    	/* repeat for reference sequence */
        bzero(fcv, ipow(4, ksize) * sizeof(double));
        bzero(fcv_1, ipow(4, (ksize-1)) * sizeof(double));
        bzero(fcv_2, ipow(4, (ksize-2)) * sizeof(double)); 
	/* compute frequencies */
	if (ambig)
    	    compute_kstring_frequencies(r, rlen, fcv, fcv_1, fcv_2, ksize);
	else
    	    compute_freqs012(r, rlen, fcv, fcv_1, fcv_2, ksize);

    	/* compute corrected vector */
	b = calloc(sizeof(double), ipow(4, ksize));
	compute_corrected_CV(b, fcv, fcv_1, fcv_2, ksize); 

    	/* we can release memory now */
	free(fcv); free(fcv_1); free(fcv_2);
	
    	d = compute_k_distance(a, b, ksize);
	printf("Distance among sequences %s and %s is %1.9f\n", 
	        infile, reference, d);
    }
    return 0;
}


/**
 * help the user with some usage information
 */
void usage() {
    fprintf(stderr, "Usage: compdist -i infile -o outfile -k k-size -a -c\n");
    fprintf(stderr, "       -a toggles consideration of ambiguous codes\n");
    fprintf(stderr, "       -c toggles correction for random mutations\n");
    exit(1);
}

/**
 * print an error message to stdout and die
 */
void error(char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(1);
}

/**
 * read input sequence from the specified file.
 *
 *	This function opens the specified input file and reads in
 * the sequence in FASTA format, allocating memory as needed to
 * hold it.
 *
 *	On exit, seq will contain a newly allocated array with the
 * sequence, and slen the sequence length.
 */
void read_sequence(char *infile, char **seq, int *slen)
{
    char *s;
    int sl, i, memsize;
    char line[BUFSIZ+1];
    FILE *in;
        
    if ((infile[0] == '-') && (infile[1] == '\0'))
      in = stdin;
    else
      if ((in = fopen(infile, "r")) == NULL)
          error("Could not open input file");

    do {
      s = fgets(line, BUFSIZ, in);
      if (s == NULL) error("No  FASTA sequences in input file");
    } while ((line[0] != '>') && (line[0] != ';'));

    /* we'll ignore sequence name for the time being */
    
    /* we are ready to read in the sequence */
    memsize = BUFSIZ;
    s = malloc(memsize+1 * sizeof(char));
    sl = 0; s[sl] = '\0';
    if (s == NULL) error("Not enough memory");
    
    while (fgets(line, BUFSIZ, in) != NULL) {
    	if (line[0] == ';')	/* comment line */
	    continue;
	if (line[0] == '>')	/* new sequence */
	    break;
	/* add line to sequence */
	for (i = 0; line[i] != '\0'; i++) {
	    if (index(nt, line[i]) != NULL) {
	      if (sl >= memsize) {
	          memsize += BUFSIZ;
		  if ((s = realloc(s, memsize+1 * sizeof(char))) == NULL)
		      error("Not enough memory");
	      }
	      s[sl++] = line[i];
	    }
	}
    }
    s[sl] = '\0';
    *seq = s;
    *slen = sl;
    if (in != stdin) fclose(in);
}

#ifdef PRETTYPRINT3X3
/*
 * a function to print a 3-char long words CV
 * more neatly
 *
 * this is a legacy of early debugging
 */
void print_cv(char *outfile, double *cv, int ksize)
{
    int i, j, k;
    FILE *out;
    
    if ((outfile[0] == '-') && (outfile[1] == '\0'))
    	out = stdout;
    else
        if ((out = fopen(outfile, "w+")) == NULL) 
	    error("Cannot open output file");
    
    /* this is for debugging only, that's why we cheat */
    for (i = 0; i < 4; i++) {
        fprintf(out, "%c\n", nt[i]);
        for (j = 0; j < 4; j++) {
	    fprintf(out, "   %c", nt[j]);
            for (k = 0; k < 4; k++) {
                fprintf(out, "\t%.3f", cv[(i*4*4)+(j*4)+k]);
            }
            fprintf(out, "\n");
        }
	fprintf(out, "\t A\t T\t G\t C\n");
    }
}
#endif

/**
 * Compute power of an integer number
 */
#ifdef STRAIGHT_IPOW
long int ipow(int b, unsigned int exp)
{
    long int p = 0;

    if ((exp != 0) && (b != 0)) {
        p = 1;
        for (; exp > 0; exp--)
	    p *= b;
    }
    return p;
}
#else
long int ipow(int x, unsigned int n)
{
    int aux = 1;
    while (n > 0) {
	if (n & 1) {		//odd ? 
	    aux *= x;
	    if (n == 1)
		return aux;
	}
	x *= x;
	n /= 2;
    }
    return aux;
}
#endif


/**
 * Compute composition vectors (no ambiguity)
 *
 *	This subroutine computes the compositions vectors of 
 * k size, (k-1) size and (k-2) size words for a given sequence
 *
 */
void
compute_freqs012(char *s, int slen,
		 double *f, double *f1, double *f2, int ksize)
{
    int i, pow4, offset, j, nt;
    int *coords;
    double norm;
    
    if ((coords = calloc(sizeof(int), ksize)) == NULL)
    	error("Not enough memory");

    /* slen does not include the trailing '\0' */
    for (i = 0; i <= (slen - ksize); i++) {
	/* run over all the sequence a sliding window of ksize */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize; j++) {
	    switch (s[i+j]) {
	    case 'A':
		offset = 0;
		break;
	    case 'T':
		offset = 1;
		break;
	    case 'G':
		offset = 2;
		break;
	    case 'C':
		offset = 3;
		break;
	    default:
		error("Sequence contains ambiguity codes");
		break;		/* reject ambiguity codes */
	    }
	    coords[j] = offset;
	}
#ifdef READABLE
	/* this k-string */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize; j++) {
	    offset += coords[(ksize - 1) - j] * pow4;
	    pow4 <<= 2;		/* pow *= 4 */
	}
	f[offset]++;
	/* k-1 substring */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize - 1; j++) {
	    offset += coords[(ksize - 2) - j] * pow4;
	    pow4 <<= 2;		/* pow *= 4 */
	}
	f1[offset]++;
	/* k-2 substring */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize - 2; j++) {
	    offset += coords[(ksize - 3) - j] * pow4;
	    pow4 <<= 2;		/* pow *= 4 */
	}
	f2[offset]++;
#else
	/* this is an alternate, more efficient implementation
	   where I collate all computations in one loop: the
	   trick is to add first the tails of w and w-1, then
	   proceed with the substring common to w, w-1 and w-2: */
	{
	    int offset1, offset2;
	    offset = offset1 = offset2 = 0;
	    offset += coords[ksize - 1];
	    offset += coords[ksize - 2] * 4;
	    offset1 += coords[ksize - 2];
	    pow4 = 1;
	    for (j = 0; j < ksize - 2; j++) {
	        offset2 += coords[(ksize - 3) - j] * pow4;
	        pow4 <<= 2;		/* pow *= 4 */
		offset1 += coords[(ksize - 3) - j] * pow4;
		offset += coords[(ksize - 3) - j] * (pow4 << 2);
	    }
	    f[offset]++;
	    f1[offset1]++;
	    f2[offset2]++;
	}
#endif

    }
    /* At this point we still have three substrings to count: */
    /* last coords[1..n] [1..n-1] and [2..n] (remember coords's 0-offset) */
    /* k-1 substring */
    pow4 = 1;
    offset = 0;
    for (j = 1; j < ksize; j++) {
	offset += coords[(ksize-1) - (j-1)] * pow4;
	pow4 *= 4;
    }
    f1[offset]++;
    /* k-2 substrings */
    pow4 = 1;
    offset = 0;
    for (j = 1; j < ksize - 1; j++) {
	offset += coords[(ksize - 2) - (j-1)] * pow4;
	pow4 *= 4;
    }
    f2[offset]++;
    pow4 = 1;
    offset = 0;
    for (j = 2; j < ksize; j++) {
	offset += coords[(ksize-1) - (j-2)] * pow4;
	pow4 *= 4;
    }
    f2[offset]++;

    free(coords);

    /*
     * convert counts to freqs: 
     *	f(u) = n(u) / (N - k + 1)
     * Chan, Chan & Wang.
     */
    pow4 = ipow(4, ksize-2);
    for (i = 0; i < pow4; i++) {
    	f[i] /= slen - ksize + 1;
	f1[i] /= slen - (ksize - 1) + 1;
	f2[i] /= slen - (ksize - 2) + 1;
    }
    pow4 <<= 2;
    for ( ; i < pow4; i++) {
    	f[i] /=  slen - ksize + 1;
	f1[i] /= slen - (ksize - 1) + 1;
    }
    pow4 <<= 2;
    for ( ; i < pow4; i++) {
    	f[i] /= slen - ksize + 1;
    }
    /**/
}




/*
 * **************************
 * New attempt. Let's further subdivide the problem
 * and accept ambiguity.
 * **************************
 */

/**
 * Compute composition vectors (allowing ambiguity)
 *
 *	This subroutine computes the compositions vectors of 
 * k size, (k-1) size and (k-2) size words for a given sequence
 *
 */
void
compute_kstring_frequencies(char *s,	/* sequence */
			    int ssize,	/* sequence length */
			    double *fk,	/* frequencies of strings of length k */
			    double *fk_1,	/* freqs. of strings of length k-1 */
			    double *fk_2,	/* freqs. of strings of length k-2 */
			    int k)		/* substring length */
{
    int i, pow4;
  
    if (k < 3)
	return;
    /* slen does not count the trailing '\0' */
    for (i = 0; i <= (ssize - k); i++) {
	add_kstring(&s[i], k, fk);
	add_kstring(&s[i], k - 1, fk_1);
	add_kstring(&s[i], k - 2, fk_2);
    }
    /* at this point we still need to count 3 k-strings */
    add_kstring(&s[i], k - 1, fk_1);
    add_kstring(&s[i], k - 2, fk_2);
    add_kstring(&s[++i], k - 2, fk_2);
    /* and there we are */

    /*
     * convert counts to freqs: 
     *	f(u) = n(u) / (N - k + 1)
     * Chan, Chan & Wang.
     */
    pow4 = ipow(4, k-2);
    for (i = 0; i < pow4; i++) {
    	fk[i] /=   ssize - k + 1;
	fk_1[i] /= ssize - (k - 1) + 1;
	fk_2[i] /= ssize - (k - 2) + 1;
    }
    pow4 <<= 2;
    for ( ; i < pow4; i++) {
    	fk[i] /=   ssize - k + 1;
	fk_1[i] /= ssize - (k -1) + 1;
    }
    pow4 <<= 2;
    for ( ; i < pow4; i++) {
    	fk[i] /= ssize - k + 1;
    }
    /* */
    return;

}

/*
 * increase the frequency count of k-string 's' (of size k) 
 * and store it in the specified table t.
 */
add_kstring(s, k, t)
char *s;			/*
				 * s is a pointer to the k-string to count, it needs not be
				 * null-terminated and might happen to be shorter than k by
				 * error, so we should check
				 */
int k;				/* length of the k-string */
double *t;			/*
				 * t is a table to store the frequencies into, it must be of
				 * a hypercube of size 4**k (sizeof(alphabet)**k).
				 */
{
    double *hyperkstr;
    int i, sym;
    double p;
    
    /*
     * we want to account for ambiguity
     * one way to achieve this is by allocating an array 
     * sizeof(alphabet) * k
     * then for each position in the string we store the
     * corresponding frequencies of each nucleotide in the
     * corresponding positions in the table
     *
     * Then we process the hyper-k-string to update the
     * frequency of each combination.
     */
    /* 1. First fill in the table of symbol frequencies per position */
    if ((hyperkstr = calloc(sizeof(double), 4*k)) == NULL) {
	/* memory error */
	error("Not enough memory");
    }

    for (i = 0; i < k; i++) {
	switch (s[i]) {
	case '\0':
	    /* the string is shorter than k, nothing to count */
	    return;
	case 'A':
	    sym = A;
	    hyperkstr[(i * 4) + sym] = relfreq[sym];
 	    break;
	case 'T':
	    sym = T;
	    hyperkstr[(i * 4) + sym] = relfreq[sym];
 	    break;
	case 'G':
	    sym = G;
	    hyperkstr[(i * 4) + sym] = relfreq[sym];
	    break;
	case 'C':
	    sym = C;
	    hyperkstr[(i * 4) + sym] = relfreq[sym];
	    break;
	case 'R':
	    sym = R;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    break;
	case 'Y':
	    sym = Y;
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    break;
	case 'M':
	    sym = M;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    break;
	case 'K':
	    sym = K;
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	case 'S':
	    sym = S;
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    break;
	case 'W':
	    sym = W;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	case 'B':
	    sym = B;
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	case 'D':
	    sym = D;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	case 'H':
	    sym = H;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	case 'V':
	    sym = V;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    break;
	case 'N':
	    sym = N;
	    hyperkstr[(i * 4) + A] = relfreq[sym];
	    hyperkstr[(i * 4) + C] = relfreq[sym];
	    hyperkstr[(i * 4) + G] = relfreq[sym];
	    hyperkstr[(i * 4) + T] = relfreq[sym];
	    break;
	default:		/* not a valid k-string to count */
	    return;
	}
    }

    /* 2. At this point we have a hyper-k-string with the relative
       significance of each position, e.g. for ATRYG
          +---+---+---+---+---+
       A  | 1 |   |0.5|   |   |
          +---+---+---+---+---+
       T  |   | 1 |   |0.5|   |
          +---+---+---+---+---+
       C  |   |   |   |0.5|   |
          +---+---+---+---+---+
       G  |   |   |0.5|   | 1 |
          +---+---+---+---+---+
            A   T   R   Y   G

       we now need to find the positions in the table for the
       k-strings ATATG, ATACG, AGTG and ATGCG and store the
       corresponding value (0.25) on each of them four.

       To do this we need to compute the position in t for each
       which would be
       (A*4^4)+(T*4^3)+(R*4^2)+(Y*4^1)+(G*4^0)
     */

    add_hyperkstr(hyperkstr, k, t, 1., 0);
    free(hyperkstr);
}


/*
 * add the frequencies of a hyper k-string to the
 * appropriate cells in the composition vector,
 * accounting for uncertainty associated to ambiguity
 * codes.
 */
void add_hyperkstr(h, k, t, prob, pos)
double *h;
int k;
double *t;
double prob;
int pos;
{
    int h_offset, t_offset, i;
    double p;

    /*
     * Now we process the hyper-k-str using a recursive algorithm
     * to cover all combinations and update the corresponding
     * frequencies in t.
     *
     * For each combination, the complete frequency will be
     * p = p1 * p2 * .. * pk
     *
     * We can prune the combinations faster if we do not pursue
     * any that contains a 0. frequency.
     *
     * Since we have a combinatorial problem and k is not fixed
     * we can only use a recursive algorithm (or unrecurse it).
     *
     * Since we will be exploding a two dimensional array to find
     * coordinates in a hypercube, some complicated math will be 
     * needed here. Watch out.
     *
     * Caveat emptor!
     */
    k--;
    if (k < 0)
	return;			/* we are done */
    t_offset = ipow(4, k);
    /* check hyperkstr[pos][1..4] */
    for (i = 0; i < 4; i++) {
	p = h[i];
	if (p > 0.) {	
	    if (k > 0)
	        /* move one position up and repeat */
	        add_hyperkstr(h+4, k, t, prob * p, pos + (t_offset * i));
	    else {
	        /* we are done, update frequency table and return */
	        t[pos + (t_offset * i)] += prob * p;
	    }
	}
    }
}

/*
   We can rewrite the above procedure to make it more efficient
   by noting that we must do k-2 before we do k-1 and k.
   If we start from low offset and go higher we could swipe the
   three kstrings at once, at the cost of losing generality.
   
   This requires some care for the boundary cases and to avoid
   duplicate counts.
   
   To be done later on.
*/


void print_kfreqs(char *outfile, double *fcv, int k)
{
    int i, j, index, nelems;
    FILE *out;
#if DEBUG
#define FOURBYLINE 1
#else
#define FOURBYLINE 0
#endif

    if ((outfile[0] == '-') && (outfile[1] == '\0'))
    	out = stdout;
    else
        if ((out = fopen(outfile, "a+")) == NULL) 
	    error("Cannot open output file");
    
    fputs("\n", out);
    nelems = ipow(4, k);
    for (i = 0 ; i < nelems; i++) {
    	/* compute k-address of this element */
	for (j = 1; j <= k; j++) {
	    index = (i / ipow(4, k-j)) % 4;
	    fputc(nt[index], out);
	}
	fprintf(out, ": % .5f ", fcv[i]);
	if (FOURBYLINE == 1) {
	    if ((i % 4) == 3) fprintf(out, "\n");
	}
	else
	    fputc('\n', out);
    }
    if (out != stdout)
    	fclose(out);
}


/*
 * Those were raw counts.
 *
 *	For serious work we want to use corrected frequencies that account
 * for evolutionary processes removing random background.
 *
 * Chan, Chan, Wang (2009):
 *
 *	denoising formulas:
 *
 *	q(LwR) = (f(Lw) · f(wR)) / f(w)
 *
 *	f = observed frequency
 *	q = predicted frequency (noise)
 *	w = word
 *	L = prefix residue
 *	R = suffix residue
 *
 *	JR Note: if f(w) == 0, then f(Lw) and f(wR) must be 0 and the
 *	predicted frequency must be 0 too.
 *
 *	signal to noise ratio (significance):
 *
 *	c(u) = (f(u) - q(u)) / q(u)   <<if q(u) != 0, else c(u) = 0>>
 *
 *	JR Note: here, if q(u) == 0 then either f(w), f(Lw) or f(wR)
 *	must be zero, and therefore f(LwR) must be zero too.
 *
 *	distance for two taxa C1 and C2:
 *
 *	d(C1, C2) = (1 - cos<c1, c2>) / 2
 *
 *	cos<c1,c2> = (c1 · c2) / (||c1|| · ||c2||)
 *
 *	||c|| = L2-norm c = (sum_i=1..n x_i^2)^1/2
 *
 *	Chan, R. H., Chan, T. H., and Wang, R. W. (2009) Composition
 *	Vector Method based on Maximum Entropy Principle for Sequence
 *	Comparison. Chinese University of Hong Kong, Department of 
 *	Mathematics, Research report 2009-05(367).
 *	ftp://ftp.math.cuhk.edu.hk/report/2009-05.ps.Z
 *
 *
 *
 * Gao, Qi, Hao (2004):
 *
 *	f(a1..ak) = observed frequency
 *	f0(a1..ak) = expected frequency
 *	a(a1..ak) = normalized frequency
 *
 *
 *                    f(a1..ak-1) · f(a2..ak)      (L-K+1) · (L-K+3)
 *	f0(a1..ak) = ------------------------- · --------------------
 *	                     f(a2..ak-1)              (L-k+2)^2
 *
 *	where 
 *		f(a1..ak) = f(LwR) above (the complete word)
 *		f0(a1..ak) = predicted q(LwR) above (the complete word)
 *		f(a1..ak-1) = f(Lw) above (word but last letter)
 *		f(a2..ak) = f(wR) above (word but first letter)
 *		f(a2..ak-1) = f(w) above (word but extreme letters)
 *		L = Sequence length
 *		K = word length
 *
 *		When L >> K the second normalization factor can be ignored
 *		and the formula reverts to Chan-Chan-Wang.
 *
 *	"It is the difference between the actual counting result f and
 *	the predicted value f0 that really reflects the shaping role of
 *	selective evolution. Therefore, we collect
 *
 *	a(a1_ak) = (f(a1_ak) - f0(a1_ak)) / f0(a1_ak)
 *
 *	Then, let a be the CV for species A and b be the CV for species B
 *	The correlation (measured by the cosine, as in document similarities
 *	theory) between both species is
 *
 *		           sum_i=1..N a_i · b_i
 *	C(A,B) = ------------------------------------------
 *		  (sum_i=1..N a_i^2 · sum_i=1..N b_i^2)^1/2
 *
 *
 *	And then the distance between both of them would be
 *
 *	D(A,B) = (1 - C(A,B)) / 2
 *
 *	Gao, L., Qi, J. and Hao, B. (2006) Simple Markov Subtraction
 *	Essentially Improves Prokaryote Phylogeny. AAPPS Bulletin, June
 *	2006, pp. 3-7.
 *
 *	Hao, B. and Qi, J. (2004) Prokaryote phylogeny without sequence
 *	alignment: from avoidance signature to composition distance. 
 *	Journal of Bioinformatics and Computational Biology, Vol. 2,
 *	No. 1, pp. 1-19.
 *
 *
 * Hu and Wang (2001)
 *
 *	P_wk = Frequency already known
 *
 *	P0_wk = Predicted frequency (from k-1 words)
 *
 *	               P_c1..ck x P_c2..ck+1
 *	P0_c1..ck+1 = -----------------------		[Eq. 6]
 *	                     P_c2..ck
 *
 *	Which is same as above.
 *
 *	"With the frequencies of longer words, one can always obtain the 
 *	frequencies of shorter ones. On the other hand, the expected 
 *	frequencies of longer words, Eq. (6), is predicted from the 
 *	frequencies of shorter words, with no more information added.
 *	Therefore, the deviation of the measured frequencies from the 
 *	expected ones gives new information emerges only in the frequencies 
 *	of the longer words. In order to use this part of information, we 
 *	refer to the following signicance index
 *
 *	        P_wk - P0_wk
 *	I_wk = ---------------
 *		(P0_wk)^1/2
 *
 *	Which divides by the square root of the expected frequency instead.
 *	Without this square root, the normalized frequency will be 
 *	overstandarized.
 *
 *	Rui Hu and Bin Wang, (2001) Statistically significant strings are 
 *	related to regulatory elements in the promoter region of 
 *	Saccharomyces cerevisiae, Physica A290, 464-474 (2001).
 *
 *
 * Hao, Qi, Wang
 *
 *	p(a1..ak) = observed frequency
 *	p0(a1..ak) = expected frequency
 *
 *	             p(a1..ak-1) p(a2..ak)
 *	p0(a1..ak) = ---------------------
 *	                  p(a2..ak-1)
 *
 *	"This is nothing but a k-2th order Markov model...
 *
 *	"It is the difference between the actual counting result f and 
 *	the predicted value f0 that really reflects the shaping role of
 * 	selective evolution. Therefore, we collect"
 *
 *	a(a1..ak) = (f(a1..ak) - f0(a1..ak)) / max(f0(a1..ak), 1)
 *
 *	JR Note: if f0 == 1, the denominator has no meaning. Otherwise,
 *	we are ignoring the denominator, which is ODD as they do not
 *	mention it elsewhere!
 *
 *	This denominator is different from the other ones!
 *
 *	The correlation C is calculated as the cosine function, and the
 *	distance is defined as
 *
 *		D(A,B) = (1 - C(A,B)) / 2
 *
 *	Since C(A, B) may vary between -1 and 1, the distance is 
 *	normalized to the interval (0, 1). The collection of distances 
 *	for all species pairs comprises a distance matrix.
 *
 *	Hao, B., Qi, J. and Wang, B. (2003) Prokaryotic phylogeny based 
 *	on complete genomes without sequence alignment. Modern Physics 
 *	Letters B, Vol. 17, No. 2 (2003) 1-4
 *
 * Qi-Wang-Hao
 *
 *	a(a1..ak) = (p(a1..ak) - p0(a1..ak)) / p0(a1..ak) 
 *		if p0(a1..ak) != 0, else a(a1..ak) = 0
 *
 *	Same as Chan-Chan-Wang
 *
 * CCV method:
 *
 *	Similar, but use all CV from k_1 to k_2 (k_1 < k_2) sizes.
 *
 *	Wu, X., Wan, X. F., Wu, G., Xu, D., Lin, G. (2006) Phylogenetic 
 *	analysis using complete signature information of whole genomes 
 *	and clustered Neighbour-Joining method. Int J Bioinform Res Appl.
 *	2006;2(3):219-248.
 *
 *	X.Wu, X.Wan, G.Wu, D. Xu, and G.-H. Lin. Whole genome phylogeny
 * 	construction via complete composition vectors. Technical Report 
 *	TR05-06, Department of Computing Science, University of Alberta, 
 *	January 2005.
 *
 * Improved CV/CCV method
 *
 *	The true expected frequency for f(a1..ak) is
 *
 *	E[f(a1..ak)] = (N - k + 1) / 4^k
 *
 *	The variance for f(a1..ak) is
 *
 *	var[f(a1..ak)] = (N - k + 1) / 4^k ·
 *			 (1 - 1/4^k)  -
 *			 2 / 4^2k ·
 *			 (k - 1) ·
 *			 (N - 3/2 k + 1) +
 *			 2/4^k ·
 *			 sum_t=1..k-1 (N - k + 1 - t) · Jt / 4t
 *
 *	where Jt = 1 if (a1..ak-t) == (at+1..ak) else Jt = 0 for t=1..k-1
 *
 *	Normalization function is now
 *
 *	(f(a1..ak) - E(f(a1..ak)) / sqrt(var[f(a1..ak)]) 
 *
 *
 *
 *	Lu, G, Zhang, S, Fang, X (2008) An improved string composition 
 *	method for sequence comparison. BMC Bioinformatics 2008, 
 *	9(Suppl 6):S15
 *
 */

/**
 *	See Hu and Wang 2001.
 */

void compute_corrected_CV(double *a, 
			  double *f, double *f_1, double *f_2, 
			  int ksize)
{
    /* we need to compute f0 first, but the access pattern allows
       us to use a temporarily to store it. for the sake of
       readability, we will use an auxiliary pointer but one must
       keep in mind all along the algorithm this fact (actullay
       when finally computing a[i] = (f[i] - f0[i]) / f0[i] ) */
    double *f0;
   
    /* we might compute f_1 and f_2 on the fly off f
      and probably would save a lot in computation time... TOBEDONE */
    int i, n, offset;
   
    f0 = a;		/* to save space */
    n = ipow(4, ksize);
    offset = ipow(4, ksize-1);
    for (i = 0; i < n; i++) {
    	/* locate f_1(a1..ak-1), f_1(a2..ak) and f_2(a2..ak)
	 *
	 * We can locate these strings very easily and efficiently: 
	 * to see how, let's change from base 4 to base 10 with k=4:
	 *
	 *	Suppose the string f we have is at coordinates [9][7][3][6]
	 * since we are in base 10, it is located in position 9736 (well,
	 * actually 9-1, 7-1, 3-1, 6-1).
	 *
	 *	To find f(a1..ak-1) we want coordinates [9][7][3] in f_1,
	 * which are at position 973. Hence to get if from f(9736) we simply
	 * *integer* divide by the base: 9736 / 10 = 973
	 *
	 *	To find f(a2..ak) we want coordinates [7][3][6] in f_1, which
	 * are at position 736. Hence to get it from f(9736) we simply get
	 * the rest modulo 1000 (base^k-1): 9736 % (10^3) = 736.
	 *
	 *	Next we need to find f(a2..ak-1) in f_2. In this case that
	 * would be [7][3]. We get it by applying the two previous transforms
	 * i. e. (9736 % (10^3)) / 10 = 73.
	 *
	 *	In our case, the base is 4, so we use 4^(k-1) as an aid
	 */
	if (f_2[(i % offset) / 4] != 0.)
	    f0[i] = f_1[i % offset] * f_1[i / 4] / f_2[(i % offset) / 4];
	else
	    f0[i] = 0.;
	
	/* now we can just compute a */
	/* Note that we could as well do it in place over f if wanted 
	   by computing f0 as a temporary double (TODO: normalize_CV() ) */
	if (f0[i] != 0.)
	    a[i] = (f[i] - f0[i]) / f0[i];
	else
	    a[i] = 0.;
   }
   /* easy, wasn't it? */
   /* Whoa! if this works, I will be amazed at my own skill! */
}

/** 
 * compute distance
 */
double compute_k_distance(double *a, double *b, int k)
{
    int i, n;
    double sab, sa2, sb2;
    double c, d;
    
    n = ipow(4, k);
    sab = sa2 = sb2 = 0;
    for (i = 0; i < n; i++) {
    	sab += a[i] * b[i];
	sa2 += a[i] * a[i];
	sb2 += b[i] * b[i];
    }
    c = sab / sqrt(sa2 * sb2);
    d = (1 - c) / 2.0;
    return d;
    
}
