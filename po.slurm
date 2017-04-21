#!/bin/bash -l

#SBATCH -p commons
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J minisweep
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --contiguous

ulimit -c unlimited

export nproc_x=2
export nproc_y=2
export OMP_NUM_THREADS=2

# export PES_PER_NODE=2
# export NPES=$(($PES_PER_NODE * $SLURM_NNODES))
# export CPUS_PER_PE=12
# export OMP_NUM_THREADS=$CPUS_PER_PE

for DIR in $HOME/minisweep/openmp_build $HOME/minisweep/tasks_build; do
#     for nthread_e in 1 2 4 6 8 12 24; do
    for nthread_e in 8 12 24 48; do
        pushd $DIR

        export OMP_NUM_THREADS=$nthread_e

        echo "===== $DIR $nthread_e ====="
        srun --ntasks=1 --ntasks-per-node=1 --cpus-per-task=$OMP_NUM_THREADS \
            ./sweep --ncell_x  8 --ncell_y 16 --ncell_z 32 \
            --nthread_e $nthread_e --ne 64 --na 32 --nblock_z 32

        popd
        echo
    done

    echo
done
