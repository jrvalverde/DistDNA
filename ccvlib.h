#ifndef _CCVLIB_

void ccv_error(char *msg);
void compute_freqs012(char *s, int slen,
		double *f, double *f1, double *f2, int ksize);
void compute_kstring_frequencies(char *s, int ssize, double *fk, 
		double *fk_1, double *fk_2, int k);
void add_hyperkstr(double *h, int k, double *t, double prob, int pos);
void print_kfreqs(char *outfile, double *fcv, int k);
void compute_corrected_CV(double *a, 
			  double *f, double *f_1, double *f_2, 
			  int ksize);
double compute_k_distance(double *a, double *b, int k);
void make_ccv(char *infile, char *outfile,int ksize, int ambig, int correct);
double dist_ccv(char *sequence, char *reference, int ksize, int ambig, int correct);


#endif
