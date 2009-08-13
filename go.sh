#!/bin/sh
set -x
rm -f ../dat/st*
rm -f ../img/img*
g++ nbody.cc && ./a.out && idl -e '.r movnbody'

