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
#include "filter.h"

/* XXX TODO
 *
 * integrate new capabilities from libserial:
 * 
 * 1. store PID of process which instantiated the dsp function object (filter, FFT, etc.)
 *
 * 2. can then spawn threads for each function and then a call to the
 * function enqueues data to be processed which gets consumed by the
 * thread. When the thread is finishes processing the data, it can
 * send a message to the caller that the data is ready to be read (or
 * could just send the data, either way).
 *
 * 3. for this model, it probably makes sense to have erlang wrapper
 * processes for each DSP function, which implements a simple blocking
 * call API to the user. So user does something like:
 *
 * ok = dsp:start_fir_fr16(taps),
 * { ok, Result } = libdsp:fir_fr16(Data),
 * %% do something with Result....
 *  
 */

/* the Blackfin libdsp functions */
extern void fir_fr16(const fract16 input[],
		     fract16 output[],
		     int length,
		     fir_state_fr16 *filter_state);

/* NIF library resource structures */
struct fir_fr16_nif_res {
	/* cast pointer to this struct onto a resource allocated when
	 * the filter initialized. */
	/* must setup the data pointers to the appropriate blocks of
	 * data following this struct */
	fir_state_fr16 state;
	unsigned char mem[0];
};

/* resource types owned by this NIF library */
struct libdsp_priv_data {
	ErlNifResourceType *fir_fr16_res_type;
};


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
	 * { fir_state_fr16, { h, binary }, ( binary of coeffs) 
	 *                   { l, int }  (interp/decim index )
	 * }
	 *
	 */
	if (!enif_get_tuple(env, argv[0], &state_tuple_len, &state_term)) {
		printf("dsp: expected state tuple\n");
		return enif_make_badarg(env);
	}

	/* now go through the state tuple, and validate and extract
	 * the state */
	if (state_tuple_len > 3) {
		printf("dsp: invalid fir_fr16 state tuple\n");
		return enif_make_badarg(env);
	}		
	
	/* grab and validate the substate atom */
	if (!enif_get_atom(env,state_term[0],atom,20,ERL_NIF_LATIN1)) {
		printf("dsp: failed to get state atom\n");
		return enif_make_badarg(env);
	}
	if (strcmp(atom, "fir_state_fr16") != 0) {
		printf("dsp: invalid state atom %s\n", atom);
		return enif_make_badarg(env);
	}		

	/*
	 * now try to get the coeffs tuple
	 */
	if (!enif_get_tuple(env, state_term[1], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,20,ERL_NIF_LATIN1)) {
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
	if (!enif_get_tuple(env, state_term[2], 
			    &substate_tuple_len, &substate_term)) {
		printf("dsp: expected substate tuple\n");
		return enif_make_badarg(env);
	}

	/* grab and validate the substate atom */
	if (!enif_get_atom(env,substate_term[0],atom,20,ERL_NIF_LATIN1)) {
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

/* XXX TODO implement reload() and upgrade() */

static ErlNifFunc nif_funcs[] = {
	{ "fir_fr16_init", 1, fir_fr16_init_nif },
	{ "fir_fr16", 2, fir_fr16_nif },
};

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
		printf("dsp: failed to open fir_fr16_state resource type\n");
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

ERL_NIF_INIT(dsp, nif_funcs, load, NULL, NULL, unload)

/* dsp.c ends here */
