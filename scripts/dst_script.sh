#! /bin/sh

# Run like this:
# >dst_script.sh run29468_cooked.A00
# to make single dst from run29468_cooked.A00.00, run29468_cooked.A00.01, .. etc

 
# get listfile of the cooked files (ie all *cooked* with no .root appended)
ls -1 $1* | grep -v root > flist

# run dstmaker_no_cuts to make <file>.root
# then run dstmaker_2pos_cut on <file>.root to make <file>.2pos.root
    
date
echo "Running dstmaker_no_cuts -Lflist $1.root"
echo "flist contains:"
more flist
dstmaker_no_cuts -Lflist $1.root

echo "Running dstmaker_2pos_cut $1.root $1.2pos.root"
dstmaker_2pos_cut $1.root $1.2pos.root

rm flist
exit
