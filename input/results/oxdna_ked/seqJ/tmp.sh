for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
echo "J"$i
mkdir "J"$i
cp Dockerfile "J"$i"/Dockerfile"
mkdir "J"$i"/test_j"$i"_200000_1"
done

DATA=`cat files.txt`
i=0
# echo $i
while read line
do
        i=$((i + 1))
        echo "sed -i -e 's|TARGET|"$line"|' "J"$i"/Dockerfile""
done << FILE
$DATA
FILE


