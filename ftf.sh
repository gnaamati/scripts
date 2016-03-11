aaa=$1
bbb=$(echo $aaa | sed 's|::|/|g')
echo $bbb
cmd="use $aaa; print \$INC{$bbb}"
echo $cmd

cmd='use Bio::EnsEMBL::Variation::Study; print $INC{"Bio/EnsEMBL/Variation/Study.pm"'

#perl -e 'use Bio::EnsEMBL::Variation::Study; print $INC{"Bio/EnsEMBL/Variation/Study.pm"}'
perl -e '$ENV{cmd}'

#echo 'file1=' $file1
#echo 'f2=' $file2
#echo Bio::EnsEMBL::EGPipeline::Common::RunnableDB::EGSpeciesFactory | sed 's|::|/|g'


