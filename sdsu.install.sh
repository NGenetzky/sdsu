#!/bin/sh
dest=../
src=./
sub0=sdsu
data_dest=${dest}_data/$sub0

mkdir -p $data_dest
cp -t $data_dest ${src}_data/*