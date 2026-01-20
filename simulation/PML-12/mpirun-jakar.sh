# Default values
nprocs=40
fname="traditional_6_pml-off_anelastic-c2"


# Parse named arguments
for arg in "$@"; do
    case $arg in
        np=*)
            nprocs="${arg#*=}"
            ;;
        fname=*)
            fname="${arg#*=}"
            ;;
        *)
            echo "Unknown argument: $arg"
            echo "Usage: $0 [np=<number>] [fname=<filename>]"
            echo "Example: $0 np=8 fname=WholeSpace_external"
            exit 1
            ;;
    esac
done


module purge
module load compiler-rt/2024.0.0 ifort/2024.0.0 mpi/2021.13 R

mpirun -np $nprocs ../../build/./waveqlab3d ./input/${fname}.in | tee ./output/${fname}.out


