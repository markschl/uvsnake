#!/bin/bash

set -e

# install seqtool
outdir="$CONDA_PREFIX/bin"

tag=v0.3.0
prefix=https://github.com/markschl/seqtool/releases/download/$tag/seqtool-$tag

u="$(uname -s)"
case "${u}" in
    Linux*)     url=$prefix-$(arch)-unknown-linux-gnu.tar.gz;;
    Darwin*)    url=$prefix-$(arch)-apple-darwin.tar.gz;;
    *)          echo "Linux or Mac expected" >&2 && exit 1
esac

tmp=_seqtool_tmp
f=$tmp/seqtool.tar.gz
mkdir -p $tmp

wget -O $f "$url"
tar -xzf $f -C $tmp
mv $tmp/st "$outdir"

rm -R $tmp
