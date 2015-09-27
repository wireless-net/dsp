%%% @author Lumenosys Robotics <dbutter@lumenosys.com>
%%% @copyright (C) 2015, Lumenosys Robotics
%%% @doc
%%%

-module(dsp).
-export([fir_fr16_init/1, 
	 fir_fr16/2, 
         coeff_iirdf1_fr16/2,
	 iirdf1_fr16_init/1, 
	 iirdf1_fr16/2, 
	 rfft_fr16_init/1, 
	 rfft_fr16/2, 
         vecvmlt_fr16/2, 
         vecdot_fr1x32/2, 
         vecdot_fr16_sr/3, 
	 cabs_fr16/1,
	 gen_hanning_fr16/2,
	 autocoh_fr16/2,
	 crosscoh_fr16/3,
	 autocorr_fr16/2,
	 crosscorr_fr16/3,
	 histogram_fr16/4,
	 mean_fr16/1,
	 var_fr16/1,
	 max_fr16/1,
	 min_fr16/1,
	 interleave/2,
	 deinterleave/3]).

-on_load(init/0).

init() ->
    Lib = filename:join(code:priv_dir("dsp"), "libdsp"),
    ok = erlang:load_nif(Lib, 0).

fir_fr16_init(_State) ->
    exit(nif_library_not_loaded).

fir_fr16(_Handle, _InputData) ->
    exit(nif_library_not_loaded).

coeff_iirdf1_fr16(_A, _B) ->
    exit(nif_library_not_loaded).    

iirdf1_fr16_init(_State) ->
    exit(nif_library_not_loaded).

iirdf1_fr16(_Handle, _InputData) ->
    exit(nif_library_not_loaded).

rfft_fr16_init(_State) ->
    exit(nif_library_not_loaded).

rfft_fr16(_Handle, _InputData) ->
    exit(nif_library_not_loaded).

vecvmlt_fr16(_InputVecA, _InputVecB) ->
    exit(nif_library_not_loaded).

vecdot_fr1x32(_InputVecA, _InputVecB) ->
    exit(nif_library_not_loaded).

vecdot_fr16_sr(_InputVecA, _InputVecB, _SRand) ->
    exit(nif_library_not_loaded).

cabs_fr16(_InputVec) ->
    exit(nif_library_not_loaded).

gen_hanning_fr16(_Stride, _Size) ->
    exit(nif_library_not_loaded).

autocoh_fr16(_Samples, _Lags) ->
    exit(nif_library_not_loaded).

crosscoh_fr16(_SamplesA, _SamplesB, _Lags) ->
    exit(nif_library_not_loaded).

autocorr_fr16(_Samples, _Lags) ->
    exit(nif_library_not_loaded).

crosscorr_fr16(_SamplesA, _SamplesB, _Lags) ->
    exit(nif_library_not_loaded).

histogram_fr16(_Samples, _MaxSample, _MinSample, _BinCount) ->
    exit(nif_library_not_loaded).

mean_fr16(_Samples) ->
    exit(nif_library_not_loaded).

var_fr16(_Samples) ->
    exit(nif_library_not_loaded).

max_fr16(_Samples) ->
    exit(nif_library_not_loaded).

min_fr16(_Samples) ->
    exit(nif_library_not_loaded).

interleave(_BinList, _ChunkSize) ->
    exit(nif_library_not_loaded).

deinterleave(_Bin, _Count, _ChunkSize) ->
    exit(nif_library_not_loaded).


%%% @end
