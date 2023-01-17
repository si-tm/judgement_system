for i in 1 2 3 4 5 6 7 8 9 10
do
echo "B"$i
mkdir "B"$i
cp Dockerfile "B"$i"/Dockerfile"
mkdir "B"$i"/test_b"$i"_200000_1"
done

DATA=`cat files.txt`
i=0
# echo $i
while read line
do
        i=$((i + 1))
        echo "sed -i -e 's|TARGET|"$line"|' "B"$i"/Dockerfile""
done << FILE
$DATA
FILE


