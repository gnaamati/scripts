file=$1
echo "fixing file:"
echo $file

perl -p -i -e 's/Triticum_aestivum_CS42_//' $file
perl -p -i -e 's/$/;/' $file

cat $file | perl -an -E 'if ($F[2] ne 'CDS'){ print }' > temp1
cat temp1 | perl -an -E 'if ($F[0] =~ /_U_/) {$F[0] = $`.'_U'}; say join("\t", @F)' > temp2
cat temp2 | perl -an -E 'if ($F[6] eq "."){$F[6] = "+"}; say join("\t", @F)' > temp3

mv temp3 $file
