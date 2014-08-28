

#The Lumenosys Obsidian DSP NIF library port application#



### Interface description ###

### Dependencies ###

To build you will need a working installation of Erlang 17 (or
later). <br/>
Please refer to [Erlang/OTP](http://www.erlang.org) for information on building and installing Erlang/OTP.

The crc16 NIF library application is built using [rebar](https://github.com/rebar/rebar). Refer to [building rebar](https://github.com/rebar/rebar/wiki/Building-rebar) for information on building and using rebar.

### Downloading

```sh
$ git clone git://github.com/dbutter/dsp.git
```
### Configuration
#### Concepts
...
#### Files
...
### Building

Compile:

```sh
$ cd dsp
$ rebar compile
...
==> dsp (compile)
```

### Usage example
...
```sh
$ export ERL_LIBS=/path/to/dsp
$ erl -sname obsidian@ion1 -boot start_sasl -eval "application:start(dsp)"
...
```


