#!/bin/env bash

declare work_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
declare binaries_dir=${PWD}/binaries
declare bench_result=${PWD}/results/csv;
mkdir -p ${bench_result};

export OMP_NUM_THREADS=`nproc`

declare GROUP_NAME=(
    MEM_DP
    L2
    L3
    L2CACHE
    L3CACHE
)

declare maps=$(find -path './binaries/*' -prune -name "map_*");
declare stencils=$(find -path './binaries/*' -prune -name "stencil_*");

declare SIZES=(
    200
    400
    600
)

declare MAX_IT=420

cmd="likwid-perfctr -O -C S0"

set -x;

for group_name in ${GROUP_NAME[@]};
do
    for stencil in ${stencils[@]};
    do
        if ! [[ "$stencil" == *"CL"* ]] || [[ "$stencil" == *"VIENNACL"* ]];
        then
            for SIZE in ${SIZES[@]};
            do
                echo ${group_name} ${stencil##*/}" "${SIZE} ${MAX_IT}
                $cmd -m -f -G 0 -g ${group_name} ${stencil} ${SIZE} ${MAX_IT} >> ${bench_result}/${stencil##*/}.csv
            done
        else if [[ "$stencil" == *"VEXCL_"* ]];
            then
                for SIZE in ${SIZES[@]};
                do
                    stencil_name=${stencil##*/}
                    echo ${group_name} ${stencil_name}' '${SIZE} ${MAX_IT}
                    $cmd -f -G 0 -g ${group_name} ${stencil} ${SIZE} ${MAX_IT} |
                    sed 's/TABLE,/TABLE,Region VEXCL_'${stencil_name:14}'-'${SIZE}'-'${MAX_IT}',/g' \
                    >> ${bench_result}/${stencil_name}.csv
                done
            fi
        fi
    done
done

