/* dsp.c --- 
 * 
 * Filename: dsp.c
 * Description: 
 * Author: 
 * Maintainer: 
 * Created: Wed Jun  4 22:09:30 2014 (-0700)
 * Version: 
 * URL: 
 * Keywords: 
 * Compatibility: 
 * 
 */

/* Commentary: 
 * 
 * 
 * 
 */

/* Change Log:
 * 
 * 
 */

/* *This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth
 * Floor, Boston, MA 02110-1301, USA.
 */

/* Code: */

#include <stdio.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <erl_nif.h>

/* libdsp */
#include <fract_math.h>
#include <filter.h>
#include <stats.h>
#include <vector.h>
#include <window.h>
#include <limits.h>

/* XXX TODO
 *
 * 1. Use dirty nifs once they are supported on target.
 */

/* the Blackfin libdsp functions */

/*
 * Finite impulse response filer
 * 
 * The fir_fr16 function implements a finite impulse response (FIR)
 * filter.  The function generates the filtered response of the input
 * data input and stores the result in the output vector output. The
 * number of input sam- ples and the length of the output vector are
 * specified by the argument length.  The function maintains the
 * filter state in the structured variable filter_state, which must be
 * declared and initialized before calling the function. The macro
 * fir_init, defined in the filter.h header file, is available to
 * initialize the structure.
 */
extern void fir_fr16(const fract16 input[],
		     fract16 output[],
		     int length,
		     fir_state_fr16 *filter_state);

/*
 * Direct form I impulse response filter
 *
 * The iirdf1_fr16 function implements a direct form I infinite
 * impulse response (IIR) filter. It generates the filtered response
 * of the input data input and stores the result in the output vector
 * output. The number of input samples and the length of the output
 * vector is specified by the argu- ment length.
 */
extern void iirdf1_fr16(const fract16 input[], 
			fract16 output[],
			int length, 
			iirdf1_state_fr16 *filter_state);

/*
 * Convert coefficients for DF1 IIR filter
 *
 * This function transforms a set of A-coefficients and a set of
 * B-coefficients into a set of coefficients for the iirdf1_fr16
 * function (see on page 4-173), which implements an optimized, direct
 * form 1 infinite impulse response (IIR) filter.
 */
extern void coeff_iirdf1_fr16 (const float acoeff[],
			       const float bcoeff[],
			       fract16 coeff[ ], int nstages);

/* The autocoherance functions compute the autocoherence of the input
 * vec- tor samples[], which contain sample_length values. The
 * autocoherence of an input signal is its autocorrelation minus its
 * mean squared. The func- tions return the result in the output array
 * coherence[] of length lags .
 */
extern void autocoh_fr16(const fract16 samples[],
			 int sample_length,
			 int lags,
			 fract16 coherence[]);

/* The cross-coherance functions compute the cross-coherence of two
 * input vectors samples_x[] and samples_y[]. The cross-coherence is
 * the cross-correlation minus the product of the mean of samples_x
 * and the mean of samples_y. The length of the input vectors is given
 * by sample_length . The functions return the result in the array
 * coherence with lags elements.
*/ 
extern void crosscoh_fr16 (const fract16 samples_x[ ],
			   const fract16 samples_y[ ],
			   int sample_length,
			   int lags,
			   fract16 coherence[ ]);

/* The autocorrelation functions perform an autocorrelation of a
 * signal.  Autocorrelation is the cross-correlation of a signal with
 * a copy of itself.  It provides information about the time variation
 * of the signal. The signal to be autocorrelated is given by the
 * samples[] input array. The number of samples of the autocorrelation
 * sequence to be produced is given by lags .  The length of the input
 * sequence is given by sample_length .  Autocorrelation is used in
 * digital signal processing applications such as speech analysis.
 */
extern void autocorr_fr16 (const fract16 samples[],
			   int sample_length,
			   int lags,
			   fract16 correlation[]);

/* The cross-correlation functions perform a cross-correlation between
 * two signals. The cross-correlation is the sum of the scalar
 * products of the sig- nals in which the signals are displaced in
 * time with respect to one another.  The signals to be correlated are
 * given by the input vectors samples_x[] and samples_y[]. The length
 * of the input vectors is given by sample_length. The functions
 * return the result in the array correlation with lags elements.
 * Cross-correlation is used in signal processing applications such as
 * speech analysis.
 */
extern void crosscorr_fr16 (const fract16 samples_x[],
			    const fract16 samples_y[],
			    int sample_length,
			    int lags,
			    fract16 correlation[]);

/* The histogram functions compute a histogram of the input vector
 * samples[ ] that contains nsamples samples, and store the result in
 * the output vector histogram .  The minimum and maximum value of any
 * input sample is specified by min_sample and max_sample ,
 * respectively. These values are used by the function to calculate
 * the size of each bin as (max_sample – min_sample) / bin_count,
 * where bin_count is the size of the output vector histogram.  Any
 * input value that is outside the range [ min_sample, max_sample )
 * exceeds the boundaries of the output vector and is discarded.  L
 * out-of-bounds To preserve maximum checking, performance the
 * histogram_fr16 while performing function allocates a temporary work
 * area on the stack. The work area is allocated with (bin_count + 2)
 * elements and the stack may therefore overflow if the number of bins
 * is sufficiently large. The size of the stack may be adjusted by
 * making appropriate changes to the . ldf file.
*/
extern void histogram_fr16 (const fract16 samples[],
			    int histogram[],
			    fract16 max_sample,
			    fract16 min_sample,
			    int sample_length,
			    int bin_count);

/* mean
 *
 * The mean functions return the mean of the input array samples[ ].
 * The number of elements in the array is sample_length. */
extern fract16 mean_fr16(const fract16 samples[],
			 int sample_length);

/* var
 *
 * The variance functions return the variance of the elements within
 * the input vector samples[ ]. The number of elements in the vector
 * is sample_length .
 *
 * The var_fr16 function can be used to compute the variance of up to
 * 65535 input data with a value of 0x8000 before the sum a i
 * saturates.
 */
extern fract16 var_fr16(const fract16 samples[],
			int sample_length);

/*
 * Minimum value of two
 */
extern fract16 min_fr16(fract16 a, fract16 b);

/*
 * Maximum value of two
 */
extern fract16 max_fr16(fract16 a, fract16 b);

/* rms */
/* zero_cross*/

/*
 * Real Vector x Vector multiplication
 * 
 */
extern void vecvmlt_fr16(const fract16 vec_a[], 
			 const fract16 vec_b[],
			 fract16 out[], int length);
	
/*
 * Complex absolute value
 *
 * The cabs functions compute the complex absolute value of a complex
 * input and return the result.
 */
extern fract16 cabs_fr16(complex_fract16 a);

/*
 * Generate Hanning window
 *
 * The gen_hanning function generates a vector
 * containing the Hanning window. The length of the window required is
 * specified by the parameter window_size, and the parameter
 * window_stride is used to space the win- dow values within the
 * output vector hanning_window. The length of the output vector
 * should therefore be window_size*window_stride . This win- dow is
 * also known as the cosine window.
 */
extern void gen_hanning_fr16(fract16 hanning_window[],
			     int window_stride,
			     int window_size);


/* N-point radix-2 real input FFT
 *
 * This function transforms the time domain real input signal sequence
 * to the frequency domain by using the radix-2 FFT. The function
 * takes advantage of the fact that the imaginary part of the input
 * equals zero, which in turn eliminates half of the multiplications
 * in the butterfly.  The size of the input array input and the output
 * array output is fft_size , where fft_size represents the number of
 * points in the FFT. If the input data can be overwritten, the
 * optimum memory usage can be achieved by also specifying the input
 * array as the output array, provided that the mem- ory size of the
 * input array is at least 2*fft_size.  The twiddle table is passed in
 * the argument twiddle_table, which must contain at least fft_size/2
 * twiddle factors. The table is composed of +cosine and -sine
 * coefficients and may be initialized by using the func- tion
 * twidfftrad2_fr16 . For optimum performance, the twiddle table
 * should be allocated in a different memory section than the output
 * array.
 *
 * The argument twiddle_stride should be set to 1 if the twiddle table
 * was originally created for an FFT of size fft_size . If the twiddle
 * table was cre- ated for a larger FFT of size N*fft_size (where N is
 * a power of 2), then twiddle_stride should be set to N . This
 * argument therefore provides a way of using a single twiddle table
 * to calculate FFTs of different sizes.  The argument scale_method
 * controls how the function will apply scaling while computing a
 * Fourier Transform. The available options are static scaling
 * (dividing the input at any stage by 2), dynamic scaling (dividing
 * the input at any stage by 2 if the largest absolute input value is
 * greater or equal than 0.5), or no scaling.
 */
extern void rfft_fr16(const fract16 input[],
		      complex_fract16 output[],
		      const complex_fract16 twiddle_table[],
		      int twiddle_stride,
		      int fft_size,
		      int *block_exponent,
		      int scale_method);

/*
 * Generate FFT twiddle factors for radix-2 FFT
 *
 * The twidfftrad2_fr16 function calculates complex twiddle
 * coefficients for an FFT of size fft_size and returns the
 * coefficients in the vector twiddle_table . The size of the vector,
 * which is known as a twiddle table, must be at least fft_size/2. It
 * contains pairs of sine and cosine values that are used by an FFT
 * function to calculate a Fast Fourier Transform. The table generated
 * by this function may be used by any of the FFT functions cfft_fr16,
 * ifft_fr16, and rfft_fr16.
 */
extern void twidfftrad2_fr16(complex_fract16 twiddle_table[],
			     int fft_size);


/* NIF library resource structures */
struct fir_fr16_nif_res {
	/* cast pointer to this struct onto a resource allocated when
	 * the filter initialized. */
	/* must setup the data pointers to the appropriate blocks of
	 * data following this struct */
	fir_state_fr16 state;
	unsigned char mem[0];
};

struct rfft_fr16_nif_res {
	/* cast pointer to this struct onto a resource allocated when
	 * the filter initialized. */
	/* must setup the data pointers to the appropriate blocks of
	 * data following this struct */
	int twiddle_stride;
	int fft_size;
	int scale_method;
	complex_fract16 *twiddle_table;
	unsigned char mem[0];
};

struct iirdf1_fr16_nif_res {
	/* cast pointer to this struct onto a resource allocated when
	 * the filter initialized. */
	/* must setup the data pointers to the appropriate blocks of
	 * data following this struct */
	iirdf1_state_fr16 state;
	unsigned char mem[0];
};

/* resource types owned by this NIF library */
struct libdsp_priv_data {
	ErlNifResourceType *fir_fr16_res_type;
	ErlNifResourceType *iirdf1_fr16_res_type;
	ErlNifResourceType *rfft_fr16_res_type;
};

ERL_NIF_TERM
mk_atom(ErlNifEnv* env, const char* atom)
{
  ERL_NIF_TERM ret;

  if(!enif_make_existing_atom(env, atom, &ret, ERL_NIF_LATIN1))
    {
      return enif_make_atom(env, atom);
    }

  return ret;
}

/* the fir_fr16_init_nif interface function */
static ERL_NIF_TERM fir_fr16_init_nif(ErlNifEnv *env, 
				      int argc, 
				      const ERL_NIF_TERM argv[])
{
	const ERL_NIF_TERM *state_term;
	int state_tuple_len;
	char atom[80];
	const ERL_NIF_TERM *substate_term;
	int substate_tuple_len;
	ERL_NIF_TERM res_term;

	/* filter state terms */
	ErlNifBinary h_term;
	int l;

	struct fir_fr16_nif_res *res;

	struct libdsp_priv_data *priv = 
		(struct libdsp_priv_data *)enif_priv_data(env);

	/*
	 * State tuple is of the form:
	 *
	 * {{ h, binary }, ( binary of coeffs) 
	 *  { l, int }  (interp/decim index )
	 * }
	 *
	 */
	if (!enif_get_tuple(env, argv[0], &state_tuple_len, &state_term)) {
		printf("dsp: expected state tuple\n");
		return enif_make_badarg(env);
	}
	if (state_tuple_len > 2) {
		printf("dsp: invalid fir_fr16 state tuple\n");
		return enif_make_badarg(env);
	}		

	/* now go through the state tuple, and validate and extract
	 * the state */
	
	/* /\* grab and validate the substate atom *\/ */
	/* if (!enif_get_atom(env,state_term[0],atom,80,ERL_NIF_LATIN1)) { */
	/* 	printf("dsp: failed to get state atom\n"); */
	/* 	return enif_make_badarg(env); */
	/* } */
	/* if (strcmp(atom, "fir_state_fr16") != 0) { */
	/* 	printf("dsp: invalid state atom %s\n", atom); */
	/* 	return enif_make_badarg(env); */
	/* }		 */

	/*
	 * now try to get the coeffs tuple
	 */
	if (!enif_get_tuple(env, state_term[0], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get state atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "h") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		
	
	/* so far so good, grab the state value */	
	if (!enif_inspect_binary(env, substate_term[1], &h_term)) {
		printf("dsp: expected coeffs binary\n");
		return enif_make_badarg(env);
	}

	/* 
	 * now try to get the interp index tuple 
	 */
	if (!enif_get_tuple(env, state_term[1], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get state atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "l") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		

	if (!enif_get_int(env, substate_term[1], &l)) {
		printf("dsp: expected int l term\n");
		return enif_make_badarg(env);
	}

	/* set up our resource */
	/* note: we allocate a big enough block of memory to hold the
	 * state structure, coefficients, and delay line */
	res = enif_alloc_resource(priv->fir_fr16_res_type, 
				  sizeof(struct fir_fr16_nif_res) + 
				  h_term.size +
				  h_term.size/*dly line same length as coeffs*/);
	
	/* initialize the state struct ... */
	fir_init(res->state, 
		 (fract16 *)&res->mem[0], 
		 (fract16 *)&res->mem[h_term.size], 
		 h_term.size >> 1, 
		 l);

	/* copy in the coefficients and init the delay line */
	memcpy((unsigned char *)res->state.h, h_term.data, h_term.size);
	memset((unsigned char *)res->state.d, 0, h_term.size/*same length as coeffs*/);

	/* create the resource for this NIF function */
	res_term = enif_make_resource(env, res);

	/* this tells Erlang to garbage collect the resource when it
	 * is no longer referenced */
	enif_release_resource(res);
	
	/* make tuple { ok, Handle } and return it */
	return enif_make_tuple2(env, 
				enif_make_atom(env, "ok"),
				res_term);
}

/* the fir_fr16_nif interface function */
static ERL_NIF_TERM fir_fr16_nif(ErlNifEnv *env, 
				 int argc, 
				 const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	struct libdsp_priv_data *priv;
	struct fir_fr16_nif_res *res;

	/* grab the private data */
	priv = (struct libdsp_priv_data *)enif_priv_data(env);

	/* first arg is the resource handle */
	if (!enif_get_resource(env, 
			       argv[0], 
			       priv->fir_fr16_res_type, (void **)&res)) {
		printf("dsp: expected fir_fr16 resource handle as first arg\n");
		return enif_make_badarg(env);
	}

	/* second arg is the input data */
	if (!enif_inspect_binary(env, argv[1], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data */
	/* for non-interp/non-decimating FIR filter, output count == input count */
	if (!enif_alloc_binary(x_i_term.size, &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}

	/* do the computation */
	fir_fr16((fract16 *)x_i_term.data, /* input data */
		 (fract16 *)y_o_term.data, /* output data */
		 x_i_term.size >> 1,	   /* # of samples */
		 &res->state);		   /* the filter state */

	/* all done! */
	return enif_make_binary(env,&y_o_term);
}

/* convert floating point coefficients to fr16 into vector for use
 * with Direct-form I IIR filter function */
static ERL_NIF_TERM coeff_iirdf1_fr16_nif(ErlNifEnv *env, 
					  int argc, 
					  const ERL_NIF_TERM argv[])
{
	ErlNifBinary bin;
	unsigned int a_count;
	unsigned int b_count;
	float *a_coeffs;
	float *b_coeffs;
	double val;
	ERL_NIF_TERM head, tail;
	int  cidx;
	int nstages;
	int i;

	/* TODO do error checking on coeffs and lengths */
	/* free all memory */

	/* grab the A coefficients */

	if (!enif_get_list_length(env, argv[0], &a_count)) {
		fprintf(stderr, "dsp: expected list\n");
		return enif_make_badarg(env);
	}
	if (a_count == 0) {
		fprintf(stderr, "dsp: expected non-empty list\n");
		return enif_make_badarg(env);
	}

	/* allocate memory for the A coefficients */
	a_coeffs = (float *)enif_alloc(sizeof(float) * a_count);
	if (a_coeffs == NULL) {
		fprintf(stderr, "dsp: failed to allocate memory\n");
		return mk_atom(env, "alloc_failed");
	}

	/* grab all A coefficients */

	/* get first element */
	if (!enif_get_list_cell(env, argv[0], &head, &tail)) {
		fprintf(stderr, "dsp: expected valid list\n");
		enif_free(a_coeffs);
		return enif_make_badarg(env);
	}
	cidx = 0;
	do {
		/* grab the coefficient */
		if (!enif_get_double(env, head, &val)) {
			fprintf(stderr, "dsp: expected double\n");
			enif_free(a_coeffs);
			return enif_make_badarg(env);
		}

		/* copy in the value */
		a_coeffs[cidx] = (float)val;

		/* increment our index */
		cidx++;

		/* get next element */
	} while(enif_get_list_cell(env, tail, &head, &tail));

	/* grab the B coefficients */

	if (!enif_get_list_length(env, argv[1], &b_count)) {
		fprintf(stderr, "dsp: expected list\n");
		enif_free(a_coeffs);
		return enif_make_badarg(env);
	}
	if (b_count == 0) {
		fprintf(stderr, "dsp: expected non-empty list\n");
		enif_free(a_coeffs);
		return enif_make_badarg(env);
	}

	/* allocate memory for the B coefficients */
	b_coeffs = (float *)enif_alloc(sizeof(float) * b_count);
	if (b_coeffs == NULL) {
		fprintf(stderr, "dsp: failed to allocate memory\n");
		enif_free(a_coeffs);
		return mk_atom(env, "alloc_failed");
	}

	/* grab all B coefficients */

	/* get first element */
	if (!enif_get_list_cell(env, argv[1], &head, &tail)) {
		fprintf(stderr, "dsp: expected valid list\n");
		enif_free(a_coeffs);
		enif_free(b_coeffs);
		return enif_make_badarg(env);
	}
	cidx = 0;
	do {
		/* grab the coefficient */
		if (!enif_get_double(env, head, &val)) {
			fprintf(stderr, "dsp: expected double\n");
			enif_free(a_coeffs);
			enif_free(b_coeffs);
			return enif_make_badarg(env);
		}

		/* copy in the value */
		b_coeffs[cidx] = (float)val;

		/* increment our index */
		cidx++;

		/* get next element */
	} while(enif_get_list_cell(env, tail, &head, &tail));

	/* determine the number of stages */
	nstages = a_count >> 1;

	/* do some sanity checking */
	if (b_count != ((2*nstages)+1)) {
		fprintf(stderr, "dsp: B coefficients length must be 2*NSTAGES + 1\n");
		enif_free(a_coeffs);
		enif_free(b_coeffs);
		return enif_make_badarg(env);
	}

	for (i = 0; i < a_count; i++) {
		if ((long)a_coeffs[i] > LONG_MAX) {
			printf("dsp: error, A coefficients cannot be larger than %ld\n", LONG_MAX);
			enif_free(a_coeffs);
			enif_free(b_coeffs);
			return enif_make_badarg(env);
		}

		if ((long)a_coeffs[i] < LONG_MIN) {
			printf("dsp: error, A coefficients cannot be less than %ld\n", LONG_MIN);
			enif_free(a_coeffs);
			enif_free(b_coeffs);
			return enif_make_badarg(env);
		}
	}

	for (i = 0; i < b_count; i++) {
		if ((long)b_coeffs[i] > LONG_MAX) {
			printf("dsp: error, B coefficients cannot be larger than %ld\n", LONG_MAX);
			enif_free(a_coeffs);
			enif_free(b_coeffs);
			return enif_make_badarg(env);
		}

		if (b_coeffs[i] < LONG_MIN) {
			printf("dsp: error, B coefficients cannot be less than %ld\n", LONG_MIN);
			enif_free(a_coeffs);
			enif_free(b_coeffs);
			return enif_make_badarg(env);
		}
	}

	/* allocate a binary for the converted coefficients */
	if (!enif_alloc_binary((4*nstages + 2) * sizeof(fract16), &bin)) {
		fprintf(stderr, "dsp: failed to allocate binary\n");
		enif_free(a_coeffs);
		enif_free(b_coeffs);
		return mk_atom(env, "alloc_failed");
	}

	/* init the coeffs vector */
	coeff_iirdf1_fr16(a_coeffs, b_coeffs, (fract16 *)bin.data, nstages);

	/* free all alloc'ed memory */
	enif_free(a_coeffs);
	enif_free(b_coeffs);

	return enif_make_tuple2(env, mk_atom(env, "ok"), enif_make_binary(env, &bin));
}

/* the iirdf1_fr16_init_nif interface function */
static ERL_NIF_TERM iirdf1_fr16_init_nif(ErlNifEnv *env, 
					 int argc, 
					 const ERL_NIF_TERM argv[])
{
	const ERL_NIF_TERM *state_term;
	int state_tuple_len;
	char atom[80];
	const ERL_NIF_TERM *substate_term;
	int substate_tuple_len;
	ERL_NIF_TERM res_term;

	/* filter state terms */
	ErlNifBinary c_term;
	int n;

	struct iirdf1_fr16_nif_res *res;

	struct libdsp_priv_data *priv = 
		(struct libdsp_priv_data *)enif_priv_data(env);

	/*
	 * State tuple is of the form:
	 *
	 * {{ c, binary }, (binary of coeffs) 
	 *  { n, int }     (number of stages)
	 * }
	 *
	 */
	if (!enif_get_tuple(env, argv[0], &state_tuple_len, &state_term)) {
		printf("dsp: expected state tuple\n");
		return enif_make_badarg(env);
	}
	if (state_tuple_len > 2) {
		printf("dsp: invalid fir_fr16 state tuple\n");
		return enif_make_badarg(env);
	}		

	/* now go through the state tuple, and validate and extract
	 * the state */
	
	/* /\* grab and validate the substate atom *\/ */
	/* if (!enif_get_atom(env,state_term[0],atom,80,ERL_NIF_LATIN1)) { */
	/* 	printf("dsp: failed to get state atom\n"); */
	/* 	return enif_make_badarg(env); */
	/* } */
	/* if (strcmp(atom, "iirdf1_state_fr16") != 0) { */
	/* 	printf("dsp: invalid state atom %s\n", atom); */
	/* 	return enif_make_badarg(env); */
	/* }		 */

	/*
	 * now try to get the coeffs tuple
	 */
	if (!enif_get_tuple(env, state_term[0], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get state atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "c") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		
	
	/* so far so good, grab the coeffs binary */	
	if (!enif_inspect_binary(env, substate_term[1], &c_term)) {
		printf("dsp: expected coeffs binary\n");
		return enif_make_badarg(env);
	}

	/*
	 * now try to get the number of stages
	 */
	if (!enif_get_tuple(env, state_term[1],
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get state atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "n") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}

	if (!enif_get_int(env, substate_term[1], &n)) {
		printf("dsp: expected int n term\n");
		return enif_make_badarg(env);
	}

	/* set up our resource */
	/* note: we allocate a big enough block of memory to hold the
	 * state structure, coefficients, and delay line */
	res = enif_alloc_resource(priv->iirdf1_fr16_res_type, 
				  sizeof(struct iirdf1_fr16_nif_res) + 
				  c_term.size + /* length of coeffs vector */
				  c_term.size/*dly line same length as coeffs*/);
	
	/* initialize the state struct ... */
	iirdf1_init(res->state, 
		    (fract16 *)&res->mem[0], /* memory for coeffs vector */
		    (fract16 *)&res->mem[c_term.size], /* memory for delay line */
		    n);				       /* number of stages */

	/* copy in the coefficients and init the delay line */
	memcpy((unsigned char *)res->state.c, c_term.data, c_term.size);
	memset((unsigned char *)res->state.d, 0, c_term.size/*same length as coeffs*/);

	/* create the resource for this NIF function */
	res_term = enif_make_resource(env, res);

	/* this tells Erlang to garbage collect the resource when it
	 * is no longer referenced */
	enif_release_resource(res);
	
	/* make tuple { ok, Handle } and return it */
	return enif_make_tuple2(env, 
				enif_make_atom(env, "ok"),
				res_term);
}

/* the iirdf1_fr16_nif interface function */
static ERL_NIF_TERM iirdf1_fr16_nif(ErlNifEnv *env, 
				    int argc, 
				    const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	struct libdsp_priv_data *priv;
	struct iirdf1_fr16_nif_res *res;

	/* grab the private data */
	priv = (struct libdsp_priv_data *)enif_priv_data(env);

	/* first arg is the resource handle */
	if (!enif_get_resource(env, 
			       argv[0], 
			       priv->iirdf1_fr16_res_type, (void **)&res)) {
		printf("dsp: expected iirdf1_fr16 resource handle as first arg\n");
		return enif_make_badarg(env);
	}

	/* second arg is the input data */
	if (!enif_inspect_binary(env, argv[1], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data */
	if (!enif_alloc_binary(x_i_term.size, &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}

	/* do the computation */
	iirdf1_fr16((fract16 *)x_i_term.data, /* input data */
		    (fract16 *)y_o_term.data, /* output data */
		    x_i_term.size >> 1,	   /* # of samples */
		    &res->state);		   /* the filter state */

	/* all done! */
	return enif_make_binary(env,&y_o_term);
}

/* little helper function to validate power of two fft_size */
int is_power_of_two (unsigned int x)
{
	return ((x != 0) && !(x & (x - 1)));
}

/* the rfft_fr16_init_nif interface function */
static ERL_NIF_TERM rfft_fr16_init_nif(ErlNifEnv *env, 
				       int argc, 
				       const ERL_NIF_TERM argv[])
{
	const ERL_NIF_TERM *state_term;
	int state_tuple_len;
	char atom[80];
	const ERL_NIF_TERM *substate_term;
	int substate_tuple_len;
	ERL_NIF_TERM res_term;
	ErlNifBinary twiddle_table_term;
	int fft_size, twiddle_stride, scale_method;
	struct rfft_fr16_nif_res *res;
	struct libdsp_priv_data *priv = 
		(struct libdsp_priv_data *)enif_priv_data(env);

	/*
	 * State tuple is of the form:
	 *
	 * {{ fft_size, int },
	 *  { twiddle_stride, int },
	 *  { scale_method, int }
	 * }
	 *
	 */
	if (!enif_get_tuple(env, argv[0], &state_tuple_len, &state_term)) {
		printf("dsp: expected state tuple\n");
		return enif_make_badarg(env);
	}
	if (state_tuple_len > 4) {
		printf("dsp: invalid rfft_fr16 state tuple\n");
		return enif_make_badarg(env);
	}		

	/* now go through the state tuple, and validate and extract
	 * the state */
	
	/* /\* grab and validate the substate atom *\/ */
	/* if (!enif_get_atom(env,state_term[0],atom,80,ERL_NIF_LATIN1)) { */
	/* 	printf("dsp: failed to get state atom\n"); */
	/* 	return enif_make_badarg(env); */
	/* } */
	/* if (strcmp(atom, "rfft_state_fr16") != 0) { */
	/* 	printf("dsp: invalid state atom %s\n", atom); */
	/* 	return enif_make_badarg(env); */
	/* }		 */

	/*
	 * now try to get the fft_size tuple
	 */
	if (!enif_get_tuple(env, state_term[0], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected fft_size substate tuple\n");
		return enif_make_badarg(env);
	}
	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get fft_size atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "fft_size") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		
	/* so far so good, grab the value */	
	if (!enif_get_int(env, substate_term[1], &fft_size)) {
		printf("dsp: expected fft_size value (int)\n");
		return enif_make_badarg(env);
	}
	/* check that the fft_size is a sane power of two */
	if (!is_power_of_two(fft_size)) {
		printf("dsp: invalid fft_size (%d), must be a power of two!\n", fft_size);
		return enif_make_badarg(env);
	}

	/* 
	 * now try to get the twiddle_stride tuple 
	 */
	if (!enif_get_tuple(env, state_term[1], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected twiddle_stride substate tuple\n");
		return enif_make_badarg(env);
	}
	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get twiddle_stride atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "twiddle_stride") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		
	if (!enif_get_int(env, substate_term[1], &twiddle_stride)) {
		printf("dsp: expected twiddle_stride value (int)\n");
		return enif_make_badarg(env);
	}

	/* 
	 * now try to get the scale_method tuple 
	 */
	if (!enif_get_tuple(env, state_term[2], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected scale_method substate tuple\n");
		return enif_make_badarg(env);
	}
	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,80,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get scale_method atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "scale_method") != 0) {
		printf("dsp: invalid substate atom %s\n", atom);
		return enif_make_badarg(env);
	}		
	if (!enif_get_int(env, substate_term[1], &scale_method)) {
		printf("dsp: expected scale_method value (int)\n");
		return enif_make_badarg(env);
	}

	/* set up our resource */

	/* allocate a binary for the twiddle table */
	if (!enif_alloc_binary(sizeof(complex_fract16) * (fft_size/2), &twiddle_table_term)) {
		printf("dsp: failed to allocate twiddle table binary\n");
		return enif_make_atom(env, "alloc_failure");
	}

	/* note: we allocate a big enough block of memory to hold the
	 * state structure and twiddle table */
	res = enif_alloc_resource(priv->rfft_fr16_res_type, 
				  sizeof(struct rfft_fr16_nif_res) + 
				  twiddle_table_term.size);
	
	/* initialize the state struct ... */
	res->fft_size = fft_size;
	res->twiddle_stride = twiddle_stride;
	res->scale_method = scale_method;
	res->twiddle_table = (complex_fract16 *)&res->mem[0];
	twidfftrad2_fr16(res->twiddle_table, fft_size);

	/* create the resource for this NIF function */
	res_term = enif_make_resource(env, res);

	/* this tells Erlang to garbage collect the resource when it
	 * is no longer referenced */
	enif_release_resource(res);
	
	/* make tuple { ok, Handle } and return it */
	return enif_make_tuple2(env, 
				enif_make_atom(env, "ok"),
				res_term);
}

/* the rfft_fr16_nif interface function */
static ERL_NIF_TERM rfft_fr16_nif(ErlNifEnv *env, 
				  int argc, 
				  const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	struct libdsp_priv_data *priv;
	struct rfft_fr16_nif_res *res;
	int block_exponent;

	/* grab the private data */
	priv = (struct libdsp_priv_data *)enif_priv_data(env);

	/* first arg is the resource handle */
	if (!enif_get_resource(env, 
			       argv[0], 
			       priv->rfft_fr16_res_type, (void **)&res)) {
		printf("dsp: expected rfft_fr16 resource handle as first arg\n");
		return enif_make_badarg(env);
	}

	/* second arg is the input data */
	if (!enif_inspect_binary(env, argv[1], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}
	/* check that the input is fft_size */
	if (x_i_term.size != (res->fft_size * sizeof(fract16))) {
		printf("dsp: input length must equal fft_size (%d)\n", res->fft_size);
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data */
	if (!enif_alloc_binary(res->fft_size * sizeof(complex_fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* do the computation */
	rfft_fr16((fract16 *)x_i_term.data,	     /* input data */
		  (complex_fract16 *)y_o_term.data, /* output data */
		  res->twiddle_table,  /* the twiddle table */
		  res->twiddle_stride, /* twiddle stride */
		  res->fft_size,       /* FFT size */
		  &block_exponent,     /* block exponent */
		  res->scale_method); /* scale method */

	/* all done! */
	return enif_make_tuple2(env, 
				enif_make_binary(env,&y_o_term),
				enif_make_int(env, block_exponent));
}

/* Autocoherence NIF */
static ERL_NIF_TERM autocoh_fr16_nif(ErlNifEnv *env, 
				     int argc, 
				     const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	fract16 *x;
	fract16 *y;
	int lags;

	/* first arg is vector of fr16 data samples */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* secong arg is number of lags */
	if (!enif_get_int(env, argv[1], &lags)) {
		printf("dsp: expected unsigned integer\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data (output
	 * count is equal to number of lags) */
	if (!enif_alloc_binary(lags * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* set input and output pointers to binary data buffers */
	x = (fract16 *)x_i_term.data;
	y = (fract16 *)y_o_term.data;

	/* to the computation */
	autocoh_fr16(x, x_i_term.size >> 1, lags, y);

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* Cross-coherence NIF */
static ERL_NIF_TERM crosscoh_fr16_nif(ErlNifEnv *env, 
				      int argc, 
				      const ERL_NIF_TERM argv[])
{
	ErlNifBinary x1_i_term;
	ErlNifBinary x2_i_term;
	ErlNifBinary y_o_term;
	fract16 *x1;
	fract16 *x2;
	fract16 *y;
	int lags;

	printf("%s: WARNING: EXPERIMENTAL!\n", __PRETTY_FUNCTION__);

	/* first arg is vector of fr16 data samples (X) */
	if (!enif_inspect_binary(env, argv[0], &x1_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* second arg is vector of fr16 data samples (Y) */
	if (!enif_inspect_binary(env, argv[1], &x2_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	if (x1_i_term.size != x2_i_term.size) {
		printf("dsp: both input binaries must be the same length\n");
		return enif_make_badarg(env);
	}

	/* third arg is number of lags */
	if (!enif_get_int(env, argv[2], &lags)) {
		printf("dsp: expected unsigned integer\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data (output
	 * count is equal to number of lags) */
	if (!enif_alloc_binary(lags * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* set input and output pointers to binary data buffers */
	x1 = (fract16 *)x1_i_term.data;
	x2 = (fract16 *)x2_i_term.data;
	y = (fract16 *)y_o_term.data;

	/* to the computation */
	crosscoh_fr16(x1, x2, x1_i_term.size >> 1, lags, y);

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* Autocorrelation NIF */
static ERL_NIF_TERM autocorr_fr16_nif(ErlNifEnv *env, 
				      int argc, 
				      const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	fract16 *x;
	fract16 *y;
	int lags;

	printf("%s: WARNING: EXPERIMENTAL!\n", __PRETTY_FUNCTION__);

	/* first arg is vector of fr16 data samples */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* secong arg is number of lags */
	if (!enif_get_int(env, argv[1], &lags)) {
		printf("dsp: expected unsigned integer\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data (output
	 * count is equal to number of lags) */
	if (!enif_alloc_binary(lags * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* set input and output pointers to binary data buffers */
	x = (fract16 *)x_i_term.data;
	y = (fract16 *)y_o_term.data;

	/* to the computation */
	autocorr_fr16(x, x_i_term.size >> 1, lags, y);

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* Cross-correlation NIF */
static ERL_NIF_TERM crosscorr_fr16_nif(ErlNifEnv *env, 
				       int argc, 
				       const ERL_NIF_TERM argv[])
{
	ErlNifBinary x1_i_term;
	ErlNifBinary x2_i_term;
	ErlNifBinary y_o_term;
	fract16 *x1;
	fract16 *x2;
	fract16 *y;
	int lags;

	/* first arg is vector of fr16 data samples (X) */
	if (!enif_inspect_binary(env, argv[0], &x1_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* second arg is vector of fr16 data samples (Y) */
	if (!enif_inspect_binary(env, argv[1], &x2_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	if (x1_i_term.size != x2_i_term.size) {
		printf("dsp: both input binaries must be the same length\n");
		return enif_make_badarg(env);
	}

	/* third arg is number of lags */
	if (!enif_get_int(env, argv[2], &lags)) {
		printf("dsp: expected unsigned integer\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data (output
	 * count is equal to number of lags) */
	if (!enif_alloc_binary(lags * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* set input and output pointers to binary data buffers */
	x1 = (fract16 *)x1_i_term.data;
	x2 = (fract16 *)x2_i_term.data;
	y = (fract16 *)y_o_term.data;

	/* to the computation */
	crosscorr_fr16(x1, x2, x1_i_term.size >> 1, lags, y);

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* the hostogram_fr16_nif interface function */
/* NOTE: the output vector is a binary of uint32_t values, not fract16!! */
static ERL_NIF_TERM histogram_fr16_nif(ErlNifEnv *env, 
				       int argc, 
				       const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	fract16 min_sample, max_sample;
	int bin_count;

	/* first arg is vector of fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* second arg is max sample value */
	if (!enif_get_int(env, argv[1], (int *)&max_sample)) {
		printf("dsp: expected integer\n");
		return enif_make_badarg(env);
	}

	/* arg 3 is min sample value */
	if (!enif_get_int(env, argv[2], (int *)&min_sample)) {
		printf("dsp: expected integer\n");
		return enif_make_badarg(env);
	}

	/* arg 4 is bin count */
	if (!enif_get_int(env, argv[3], &bin_count)) {
		printf("dsp: expected integer\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data */
	if (!enif_alloc_binary(bin_count * sizeof(int), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}

	/* do the computation */
	histogram_fr16((fract16 *)x_i_term.data, /* input vec */
		       (int *)y_o_term.data, /* output vec  */
		       max_sample,
		       min_sample,
		       x_i_term.size >> 1, /* sample_length */
		       bin_count);
	
	/* all done! */
	return enif_make_binary(env, &y_o_term);
}


/* the vecvmlt_fr16_nif interface function */
static ERL_NIF_TERM vecvmlt_fr16_nif(ErlNifEnv *env, 
				     int argc, 
				     const ERL_NIF_TERM argv[])
{
	ErlNifBinary xa_i_term;
       	ErlNifBinary xb_i_term;
	ErlNifBinary y_o_term;

	/* first arg is vector A of fract16s */
	if (!enif_inspect_binary(env, argv[0], &xa_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* second arg is vector B of fract16s */
	if (!enif_inspect_binary(env, argv[1], &xb_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* sanity check on lengths */
	if (xa_i_term.size != xb_i_term.size) {
		printf("dsp: vector lengths must be equal!\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data (output
	 * length equals length of input vectors) */
	if (!enif_alloc_binary(xa_i_term.size, &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* do the computation */
	vecvmlt_fr16((fract16 *)xa_i_term.data, /* input vec A */
		     (fract16 *)xb_i_term.data, /* input vec B */
		     (fract16 *)y_o_term.data,	/* output vec  */
		     xa_i_term.size/sizeof(fract16));/* vector len  */

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* the cabs_fr16_nif interface function */
/* NOTE: this function supports single or vector arguments (i.e., the
 * binary can contain a single sample, or a be an entire vector of
 * samples */
static ERL_NIF_TERM cabs_fr16_nif(ErlNifEnv *env, 
				  int argc, 
				  const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	ErlNifBinary y_o_term;
	complex_fract16 *x;
	fract16 *y;
	int length;
	int i;

	/* struct libdsp_priv_data *priv; */

	/* /\* grab the private data *\/ */
	/* priv = (struct libdsp_priv_data *)enif_priv_data(env); */

	/* first arg is vector of complex_fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* compute complex_fract16 vector length */
	length = x_i_term.size / sizeof(complex_fract16);

	/* next, allocate an output binary for the result data (output
	 * count is equal to input count, but size of elements is
	 * fract16 */
	if (!enif_alloc_binary(length * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* set input and output pointers to binary data buffers */
	x = (complex_fract16 *)x_i_term.data;
	y = (fract16 *)y_o_term.data;
	
	/* do the computation */
	for (i = 0; i < length; i++) {
		y[i] = cabs_fr16(x[i]);
	}

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* the gen_hanning_fr16_nif interface function */
static ERL_NIF_TERM gen_hanning_fr16_nif(ErlNifEnv *env, 
					 int argc, 
					 const ERL_NIF_TERM argv[])
{
	ErlNifBinary y_o_term;
	int stride;
	int size;

	/* struct libdsp_priv_data *priv; */

	/* /\* grab the private data *\/ */
	/* priv = (struct libdsp_priv_data *)enif_priv_data(env); */

	if (!enif_get_int(env, argv[0], &stride)) {
		printf("dsp: expected window stride value (int)\n");
		return enif_make_badarg(env);
	}

	if (!enif_get_int(env, argv[1], &size)) {
		printf("dsp: expected window size value (int)\n");
		return enif_make_badarg(env);
	}

	/* next, allocate an output binary for the result data */
	if (!enif_alloc_binary(stride * size * sizeof(fract16), &y_o_term)) {
		printf("dsp: failed to allocate output binary\n");
		return enif_make_atom(env, "alloc_failure");
	}
	
	/* generate the window */
	gen_hanning_fr16((fract16 *)y_o_term.data, stride, size);

	/* all done! */
	return enif_make_binary(env, &y_o_term);
}

/* the mean_fr16_nif interface function */
static ERL_NIF_TERM mean_fr16_nif(ErlNifEnv *env, 
				  int argc, 
				  const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	fract16 *x;
	fract16 y;

	/* first arg is vector of fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* set input and output pointers to binary data buffers */
	x = (fract16 *)x_i_term.data;
	
	/* do the computation */
	y = mean_fr16(x, x_i_term.size >> 1);

	/* all done! */
	return enif_make_int(env, y);
}

/* the var_fr16_nif interface function */
static ERL_NIF_TERM var_fr16_nif(ErlNifEnv *env, 
				 int argc, 
				 const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	fract16 *x;
	fract16 y;

	/* first arg is vector of fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* set input and output pointers to binary data buffers */
	x = (fract16 *)x_i_term.data;
	
	/* do the computation */
	y = var_fr16(x, x_i_term.size >> 1);

	/* all done! */
	return enif_make_int(env, y);
}


/* the min_fr16_nif interface function */
/* NOTE: this function supports single or vector arguments (i.e., the
 * binary can contain a single sample, or a be a vector of samples */
static ERL_NIF_TERM min_fr16_nif(ErlNifEnv *env, 
				 int argc, 
				 const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	fract16 *x;
        fract16 min;
	int i;

	/* first arg is vector of fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* set input pointer to binary data buffer */
	x = (fract16 *)x_i_term.data;
	
	/* do the computation */
	min = 0x7fff;
	for (i = 0; i < (x_i_term.size >> 1); i++) {
		min = min_fr16(min, x[i]);
	}

	/* all done! */
	return enif_make_int(env, (int)min);
}

/* the max_fr16_nif interface function */
/* NOTE: this function supports single or vector arguments (i.e., the
 * binary can contain a single sample, or a be a vector of samples */
static ERL_NIF_TERM max_fr16_nif(ErlNifEnv *env, 
				 int argc, 
				 const ERL_NIF_TERM argv[])
{
	ErlNifBinary x_i_term;
	fract16 *x;
        fract16 max;
	int i;

	/* first arg is vector of fract16s */
	if (!enif_inspect_binary(env, argv[0], &x_i_term)) {
		printf("dsp: expected input binary\n");
		return enif_make_badarg(env);
	}

	/* set input pointer to binary data buffer */
	x = (fract16 *)x_i_term.data;
	
	/* do the computation */
	max = 0x8000;
	for (i = 0; i < (x_i_term.size >> 1); i++) {
		max = max_fr16(max, x[i]);
	}

	/* all done! */
	return enif_make_int(env, (int)max);
}

/* XXX TODO implement reload() and upgrade() */

/* Loads the NIF module and initializes private data */
static int load(ErlNifEnv *env, void **priv_data, ERL_NIF_TERM load_info)
{
	struct libdsp_priv_data *priv;

	if (!(priv = enif_alloc(sizeof(struct libdsp_priv_data)))) {
		printf("dsp: failed to allocate private data!\n");
		return -1;
	}

	/* create resources */
	priv->fir_fr16_res_type = 
		enif_open_resource_type(env, 
					NULL/*module name not used*/, 
					"fir_fr16_state", 
					NULL /*no DTOR needed*/, 
					ERL_NIF_RT_CREATE, 
					NULL/*no needed*/);
	if (!priv->fir_fr16_res_type) {
		printf("dsp: failed to open fir_fr16 resource type\n");
		return -1;
	}

	priv->iirdf1_fr16_res_type = 
		enif_open_resource_type(env, 
					NULL/*module name not used*/, 
					"iirdf1_fr16_state", 
					NULL /*no DTOR needed*/, 
					ERL_NIF_RT_CREATE, 
					NULL/*no needed*/);
	if (!priv->iirdf1_fr16_res_type) {
		printf("dsp: failed to open iirdf1_fr16 resource type\n");
		return -1;
	}

	priv->rfft_fr16_res_type = 
		enif_open_resource_type(env, 
					NULL/*module name not used*/, 
					"rfft_fr16_state", 
					NULL /*no DTOR needed*/, 
					ERL_NIF_RT_CREATE, 
					NULL/*no needed*/);
	if (!priv->rfft_fr16_res_type) {
		printf("dsp: failed to open rfft_fr16 resource type\n");
		return -1;
	}


	*priv_data = priv;

	return 0;
}

static void unload(ErlNifEnv *env, void *priv_data)
{
	struct libdsp_priv_data *priv = 
		(struct libdsp_priv_data *)priv_data;
	
	/* free the private data */
	enif_free(priv);
}

static ErlNifFunc nif_funcs[] = {
	{ "fir_fr16_init", 1, fir_fr16_init_nif },
	{ "fir_fr16", 2, fir_fr16_nif },
	{ "coeff_iirdf1_fr16", 2, coeff_iirdf1_fr16_nif },
	{ "iirdf1_fr16_init", 1, iirdf1_fr16_init_nif },
	{ "iirdf1_fr16", 2, iirdf1_fr16_nif },
	{ "rfft_fr16_init", 1, rfft_fr16_init_nif },
	{ "rfft_fr16", 2, rfft_fr16_nif },
	{ "vecvmlt_fr16", 2, vecvmlt_fr16_nif },
	{ "cabs_fr16", 1, cabs_fr16_nif },
	{ "gen_hanning_fr16", 2, gen_hanning_fr16_nif },
	{ "autocoh_fr16", 2, autocoh_fr16_nif },
	{ "crosscoh_fr16", 3, crosscoh_fr16_nif },
	{ "autocorr_fr16", 2, autocorr_fr16_nif },
	{ "crosscorr_fr16", 3,  crosscorr_fr16_nif },
	{ "histogram_fr16", 4, histogram_fr16_nif },
	{ "mean_fr16", 1, mean_fr16_nif },
	{ "var_fr16", 1, var_fr16_nif },
	{ "max_fr16", 1, max_fr16_nif },
	{ "min_fr16", 1, min_fr16_nif },
};

ERL_NIF_INIT(dsp, nif_funcs, load, NULL, NULL, unload)

/* dsp.c ends here */
