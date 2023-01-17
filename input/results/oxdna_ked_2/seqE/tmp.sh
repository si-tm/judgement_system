for i in 1 2 3 4 5 6 7 8 9 10
do
echo "E"$i
# mkdir "E"$i
# cp Dockerfile "E"$i"/Dockerfile"
# mkdir "E"$i"/test_e"$i"_200000_1"
done

DATA=`cat files.txt`
i=0
# echo $i
while read line
do
        i=$((i + 1))
        echo "sed -i -e 's|TARGET|"$line"|' "E"$i"/Dockerfile""
done << FILE
$DATA
FILE


