#!/bin/bash

PROG=`basename $0`
if [ $# -ne 3 ]
then
    echo "Usage: $PROG comment ptDep weight"
    exit;
fi

function setenv(){ export $1=$2; }

setenv comment $1
setenv ptDep $2
setenv weight $3
setenv noFileToRun 100
setenv noEvents 1000

setenv DoWhat       toyFlow
setenv oname        $DoWhat_${comment}

setenv Main_DIR     `pwd`
setenv Out_DIR      $Main_DIR/outputs/${oname}/data
setenv LOG_DIR      $Main_DIR/outputs/${oname}/logs
setenv OUT_ERRORS   $Main_DIR/outputs/${oname}/errors

mkdir -p $Out_DIR
mkdir -p $OUT_ERRORS
mkdir -p $LOG_DIR

# simplify this !!!
cat << EOF > exec_toyFlow_$comment
#!/bin/bash -f
cd $Main_DIR
source setup.sh
export what=$DoWhat
export sedN=1000
export iseg=\$SLURM_ARRAY_TASK_ID
sedN=\`expr \$sedN + \${iseg}\`
export outfile=${comment}_\$sedN
export Log=$LOG_DIR/$DoWhat-\$sedN.log
./\${what} $noEvents 1000 $ptDep $weight \$outfile $Out_DIR/ \$sedN >& \$Log
cd $Main_DIR
EOF
#\${what} \$outfile \$sedN 5000000 >& \$Log
chmod +x exec_dijet_$comment
    sbatch -v --array=1-$noFileToRun exec_toyFlow_$comment -J $comment  -e $OUT_ERRORS -o $OUT_ERRORS
