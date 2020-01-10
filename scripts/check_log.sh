path=$1
notdone=0
for name in `ls ${path} | grep log`
do
    info=`tail -n 1 ${path}/$name`
    if [[ $info != *"all done"* ]]
    then
        echo $name" not done"
        notdone=1
        break
    fi
done

data_num=`ls ${path}/splited | wc -l`
result_num=`ls ${path} | grep "oo_"|wc -l`

if [ $data_num -ne $result_num ]
then
    notdone=1
    echo data_num=$data_num
    echo result_num=$result_num
fi

if [[ $notdone -eq 1 ]]
then
    echo "not done"
else
    echo "done"
fi 
