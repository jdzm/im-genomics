## GC content calculator
conda activate qcmeta
sample=AXEK_8
ref=${cada}/assemblies/clean/flye_circ_${sample}.fasta

# uses `infoseq`
cat $ref | infoseq -auto -only -name -length -pgc stdin

# Loop for our genomes
for sample in AXEK_6 AXEK_7 AXEK_8 AXEK_9; do
    echo "Calculating %GC content for $sample"
    ref=${cada}/assemblies/clean/flye_circ_${sample}.fasta
    cat $ref | infoseq -auto -only -name -length -pgc stdin
    echo ""
done
