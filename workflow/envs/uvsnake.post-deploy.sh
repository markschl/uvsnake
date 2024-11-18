#!/bin/bash

set -e

# install seqtool
outdir="$CONDA_PREFIX/bin"

tag=0.4.0-beta.2
prefix=https://github.com/markschl/seqtool/releases/download/$tag/seqtool

u="$(uname -s)"
case "${u}" in
    Linux*)     url=$prefix-$(arch)-unknown-linux-gnu.tar.xz;;
    Darwin*)    url=$prefix-$(arch)-apple-darwin.tar.xz;;
    *)          echo "Linux or Mac expected" >&2 && exit 1
esac

tmp=_seqtool_tmp
f=$tmp/seqtool.tar.xz
mkdir -p $tmp

wget -O $f "$url"
tar -xJf $f -C $tmp
mv $tmp/*/st "$outdir"

rm -R $tmp
