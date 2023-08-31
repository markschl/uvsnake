#!/usr/bin/env bash

./uvsnake test --rulegraph config quality unoise3 uparse sintax clean clean_all clean_qc |
    dot -Tpng > rulegraph.png
