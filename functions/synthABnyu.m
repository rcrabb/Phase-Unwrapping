function [ output_args ] = synthABnyu( input_args )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%[results0, state0, data0, params0] = go_MIT('params.EVAL_NAMES = {''frog1''}; params.USE_COLOR_IMAGE = 0; params.NATURAL_ILLUMINATION = 1;params.USE_INIT_Z = 0; params.INIT_Z_SIGMA = 30;');
b = input_args*1000;
output_args = localfun(b);

end

function [ localout ] = localfun(localin)

localout = localin*10/3;

end

[r0, state0, d0, params] = go_MIT('params.EVAL_NAMES = {''nyu1''}; params.USE_COLOR_IMAGE = 0; params.NATURAL_ILLUMINATION = 1;params.USE_INIT_Z = 0; params.INIT_Z_SIGMA = 1;');
[r1, state1, d1, params] = go_MIT('params.EVAL_NAMES = {''nyu1''}; params.USE_COLOR_IMAGE = 0; params.NATURAL_ILLUMINATION = 1;params.USE_INIT_Z = 1; params.INIT_Z_SIGMA = 1;');
