dsp
===

**dsp** is an Erlang NIF library that provides a collection of optimized Digital Signal Processing functions including filters (IIR, FIR), transforms (FFT, IFFT), correlation functions, statistics, and more. This library allows Erlang applications to efficiently perform signal processing on the Analog Devices Blackfin processors. It can be used for things like sensor signal processing, feature extaction, robot control, and data visualization on the Lumenosys Robotics [BMOD][1] board. 

Dependencies
------------

To build you will need a working installation of Erlang 17 (or
later). <br/>
Please refer to [Erlang/OTP](http://www.erlang.org) for information on building and installing Erlang/OTP.

This application is built using [rebar](https://github.com/rebar/rebar). Refer to [building rebar](https://github.com/rebar/rebar/wiki/Building-rebar) for information on building and using rebar.

You will also need libbfdsp, a fork of the Analog Devices DSP library for the Blackfin processor. This should be pulled in by rebar automatically.

You will also need the Analog Devices [Blackfin GNU toolchain][2] installed in your path. 

Downloading
-----------

```sh
$ git clone git://github.com/lumenosys/dsp.git
```
Building
--------

Compile:

```sh
$ cd dsp
$ make all
...
==> dsp (compile)
```

Usage example
-------------

```erlang
%%
%% Intialize Fourier Transform
%%
fft_fr16_init(FFTSize) ->
    {ok, Handle} = dsp:rfft_fr16_init({{fft_size, FFTSize}, 
				       {twiddle_stride, 1},
				       {scale_method, 1}}),
    %% generate hanning window
    Win = dsp:gen_hanning_fr16(1, FFTSize),
    %% Note: (FFTSize/2) * 2 since fract16 is 2 bytes
    {Handle, Win}.

%%% 
%%% Compute the Fourier Transform on a binary vector of Q15
%%% samples. The first argument is a tuple constructed by first
%%% calling stft_fr16_init(WindowLen, Overlap). The second argument is
%%% a binary containing the Q15 formatted samples.
%%%
fft_fr16({FFTHandle, Win}, NewData) ->
    %% window the signal
    WinSig = dsp:vecvmlt_fr16(Win, NewData),
    %% do the FFT
    {Cresult, _BlockExponent} = dsp:rfft_fr16(FFTHandle, WinSig),
    %% compute the magnitude (and keep only the real half of the FFT result)
    RealSize = ?FFT_SIZE, % Real-side size in bytes == (FFT_SIZE / 2) * 2
    <<NewMagFrame:RealSize/binary, _Img/binary>> = dsp:cabs_fr16(Cresult),
    {ok, {FFTHandle, Win}, NewMagFrame}.

```

Copyright and License
---------------------

> %CopyrightBegin%
>
> Copyright Lumenosys Robotics 2014-2015. All Rights Reserved.
>
> Licensed under the Apache License, Version 2.0 (the "License");
> you may not use this file except in compliance with the License.
> You may obtain a copy of the License at
>
>     http://www.apache.org/licenses/LICENSE-2.0
>
> Unless required by applicable law or agreed to in writing, software
> distributed under the License is distributed on an "AS IS" BASIS,
> WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
> See the License for the specific language governing permissions and
> limitations under the License.
>
> %CopyrightEnd%


[1]: https://lumenosys.com/products
[2]: http://sourceforge.net/projects/adi-toolchain/