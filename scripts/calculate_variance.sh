
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
cat $COVFILE | tail -n +2|  awk "\$6 > $THRESHOLD" | cut -f1 | while read line;
do
    echo $line
    n_var=$(bcftools mpileup  -Ov -f $REFERENCE \
    $BAMFILE -r $line |  \
    bcftools call -mv -Ov --ploidy 1 | grep -v '#' | awk '$6 > 200' | wc -l)
    echo "$line, $n_var" >> ${OUT}
done
