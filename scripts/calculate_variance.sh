
THRESHOLD=200

for i in "$@"
do
case $i in
    -b=*|--bamfile=*)
    BAMFILE="${i#*=}"
    ;;
    -l=*|--listfile=*)
    COVFILE="${i#*=}"
    ;;
    -f=*|--reference=*)
    REFERENCE="${i#*=}"
    ;;
    -t=*|--threshold=*)
    THRESHOLD="${i#*=}"
    ;;
    -o=*|--output_dir=*)
    OUT="${i#*=}"
    ;;
    *)
            # unknown option
    ;;
esac
done


BAMNAME=$(echo $(basename $BAMFILE) |  cut -d '.' -f 1)
rm -f $OUT
touch $OUT
echo "rname,snps,deep_sites" >> $OUT
cat $COVFILE | tail -n +2|  awk "\$6 > $THRESHOLD" | cut -f1 | while read line;
do
    echo $line
    n_var=$(freebayes -f $REFERENCE $BAMFILE -p 1 -r $line \
    | bcftools filter -e "QUAL <= 20" | grep -v '#' | wc -l)
    depth=$(samtools depth  $BAMFILE -r $line | awk '$3 > 5' | wc -l)
    echo "$line, $n_var, $depth" >> ${OUT}
done
