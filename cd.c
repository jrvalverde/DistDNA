/*
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


/*
 * Alphabet definition
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

char *nt = "ATGCRYMKSWBDHVN";

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




void usage();
void error(char *msg);
void read_sequence(char *infile, char **seq, int *slen);
void print_cv(char *outfile, double *cv, int ksize);
long int ipow(int b, unsigned int exp);
void compute_freqs012(char *s, int slen,
		double *f, double *f1, double *f2, int ksize);
void compute_kstring_frequencies(char *s, int ssize, double *fk, 
		double *fk_1, double *fk_2, int k);
void
add_hyperkstr(double *h, int k, double *t, double prob, int pos);




int main(int argc, char *argv[])
{
    int opt, i;
    extern char *optarg;
    extern int optind, optopt;
    int ksize, slen;
    char *infile, *outfile, *s;
    double *fcv, *fcv_1, *fcv_2;

    /* initialize values */
    opt = 0;
    ksize = slen = 0;
    infile = outfile = s = NULL;
    fcv = fcv_1 = fcv_2 = NULL;

    /* parse arguments */
    while ((opt = getopt(argc, argv, "i:o:k:")) != -1) {
    	switch (opt ) {
	    case 'i': infile = optarg; break;
	    case 'o': outfile = optarg; break;
	    case 'k': ksize = atoi(optarg); break;
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

    /* allocate frequency vectors */
    fcv = calloc(sizeof(double), ipow(4, ksize));
    fcv_1 = calloc(sizeof(double), ipow(4, (ksize - 1)));
    fcv_2 = calloc(sizeof(double), ipow(4, (ksize - 2)));
    if ((fcv == NULL) || (fcv_1 == NULL) || (fcv_2 == NULL)) {
	error("Not enough memory");
    }

    /* compute frequencies */
    compute_kstring_frequencies(s, slen, fcv, fcv_1, fcv_2, ksize);

    /* print out results for checking */
    print_cv("-", fcv, ksize);

    /* reinitialize f* */
    for (i = 0; i < ipow(4, ksize); i++) fcv[i] = 0.;
    for (i = 0; i < ipow(4, ksize-1); i++) fcv_1[i] = 0.;
    for (i = 0; i < ipow(4, ksize-2); i++) fcv_2[i] = 0.;

    compute_freqs012(s, slen, fcv, fcv_1, fcv_2, ksize);

    /* print out results for checking */
    print_cv(outfile, fcv, ksize);

    return 0;
}

void usage() {
    fprintf(stderr, "Usage: compdist -i infile -o outfile -k k-size\n");
    exit(1);
}

void error(char *msg) {
    fprintf(stderr, "Error: %s\n", msg);
    exit(1);
}

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
}


void print_cv(char *outfile, double *cv, int ksize)
{
    int i, j, k;
    FILE *out;
    
    if ((outfile[0] == '-') && (outfile[1] == '\0'))
    	out = stdout;
    else
        if ((out = fopen(outfile, "w+")) == NULL) 
	    error("Cannot open output file");
    
    /* this is for debugging only right now so we'll cheat */
    for (i = 0; i < 4; i++) {
        fprintf(out, "%c\n", nt[i]);
        for (j = 0; j < 4; j++) {
	    fprintf(out, "   %c", nt[j]);
            for (k = 0; k < 4; k++) {
                fprintf(out, "\t%.2f", cv[(i*4*4)+(j*4)+k]);
            }
            fprintf(out, "\n");
        }
	fprintf(out, "\t A\t T\t G\t C\n");
    }



}

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


void
compute_freqs012(char *s, int slen,
		 double *f, double *f1, double *f2, int ksize)
{
    int i, pow4, offset, j, nt;
    int *coords;
    
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
	/* this k-string */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize; j++) {
	    offset += coords[(ksize - 1) - j] * pow4;
	    pow4 *= 4;
	}
	f[offset]++;
	/* k-1 substring */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize - 1; j++) {
	    offset += coords[(ksize - 2) - j] * pow4;
	    pow4 *= 4;
	}
	f1[offset]++;
	/* k-2 substring */
	pow4 = 1;
	offset = 0;
	for (j = 0; j < ksize - 2; j++) {
	    offset += coords[(ksize - 3) - j] * pow4;
	    pow4 *= 4;
	}
	f2[offset]++;
    }

    /* At this point we still have three substrings to count: */
    /* last coords[1..n] [1..-n-1] and [2..n] */
    /* -1 substring */
    pow4 = 1;
    offset = 0;
    for (j = 1; j < ksize; j++) {
	offset += coords[(ksize - 1) - j] * pow4;
	pow4 *= 4;
    }
    f1[offset]++;
    /* -2 substrings */
    pow4 = 1;
    offset = 0;
    for (j = 1; j < ksize - 1; j++) {
	offset += coords[(ksize - 2) - j] * pow4;
	pow4 *= 4;
    }
    f2[offset]++;
    pow4 = 1;
    offset = 0;
    for (j = 2; j < ksize; j++) {
	offset += coords[(ksize - 1) - j] * pow4;
	pow4 *= 4;
    }
    f2[offset]++;

    free(coords);
}




/*
 * **************************
 * New atempt. Let's further subdivide the problem
 * **************************
 */


void
compute_kstring_frequencies(char *s,	/* sequence */
			    int ssize,	/* sequence length */
			    double *fk,	/* frequencies of strings of length k */
			    double *fk_1,	/* freqs. of strings of length k-1 */
			    double *fk_2,	/* freqs. of strings of length k-2 */
			    int k)		/* substring length */
{
    int i;
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

void
add_hyperkstr(h, k, t, prob, pos)
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
   we can rewrite the above procedure to make it more efficient
   by noting that we must do k-2 before we do k-1 and k.
   If we start from low offset and go higher we could swipe the
   three kstrings at once, at the cost of losing generality.
   
   To be done after debugging the above codes.
*/
