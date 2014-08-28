%%% @author Devin Butterfield <dbutter@db>
%%% @copyright (C) 2014, Devin Butterfield
%%% @doc
%%%

-module(dsp).
-export([fir_fr16_init/1, fir_fr16/2]).
-on_load(init/0).

init() ->
    Lib = filename:join(code:priv_dir("dsp"), "libdsp"),
    ok = erlang:load_nif(Lib, 0).

fir_fr16_init(_State) ->
    exit(nif_library_not_loaded).

fir_fr16(_Handle, _InputData) ->
    exit(nif_library_not_loaded).

%%% @end
