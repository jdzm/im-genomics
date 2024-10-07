conda activate tester

PASS_FOLDER=/media/bioinfoserver/Gdrive/BKPs/data_bkp/metagenomics/funghi_indoor_alpeggi_1_20210611_20210618/funghi_indoor_alpeggi_1_20210611/20210611_1115_MN24748_FAO28578_24238651/fastq_pass
#`ls -1 $PASS_FOLDER`
for barcode in `ls -1 $PASS_FOLDER`; do 
    echo "processing $barcode"
    for fastq in `ls -1 ${PASS_FOLDER}/${barcode}/*fastq`; do
        bioawk -c fastx '{print $name, length ($seq), meanqual($qual)}' $fastq >> ${barcode}.summary.lq.tsv
        # bioawk -v bc="$barcode" -c fastx '{print $name, length ($seq), meanqual($qual), bc}' $fastq >> ${barcode}.summary.lq.tsv
    done
done

# bioawk -v bc="$barcode" -c fastx '{print $name, length ($seq), meanqual($qual), bc}' $fastq | head
# bioawk -c fastx '{print $name, length ($seq), meanqual($qual)}' ${PASS_FOLDER}/${barcode} >> barcode01.summary.lq.tsv

# add this to also get gc % gc($seq)
# bioawk -c fastx 'length($seq) > 100{ print ">"$name; print $seq }'  input.fasta 