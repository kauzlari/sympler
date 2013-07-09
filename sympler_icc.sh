#!/bin/bash
export LD_LIBRARY_PATH=/simstorage/software/intel/cpp11/lib/intel64
export LANG=C
script_path=${BASH_SOURCE-$0}
abs_path="$(cd "${script_path%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
path_only="`dirname $abs_path`"
exec $path_only/sympler "$@"
