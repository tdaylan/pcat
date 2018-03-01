
thisdate=`date "+%Y%m%d%H%M%S"`
echo `python $1 $2 > thisdate.$3.log &2>1 &`


