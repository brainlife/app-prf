#!/bin/bash
module unload matlab && module load matlab/2017a

log=compiled/commit_ids.txt
true > $log
echo "/N/u/brlife/git/jsonlab" >> $log
(cd /N/u/brlife/git/jsonlab && git log -1) >> $log
echo "/N/u/brlife/git/mrTools" >> $log
(cd /N/u/brlife/git/mrTools && git log -1) >> $log
echo "/N/u/brlife/git/vistasoft" >> $log
(cd /N/u/brlife/git/vistasoft && git log -1) >> $log

#maybe dangerous
#echo "compiling with following env" >> $log
#env >> $log

cat > build.m <<END
addpath(genpath('/N/u/brlife/git/jsonlab'))
addpath(genpath('/N/u/brlife/git/mrTools'))
addpath(genpath('/N/u/brlife/git/vistasoft'))
mcc -m -R -nodisplay -d compiled main
exit
END
matlab -nodisplay -nosplash -r build
