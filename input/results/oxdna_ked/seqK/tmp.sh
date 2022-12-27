for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
echo "K"$i
mkdir "K"$i
cp Dockerfile "K"$i"/Dockerfile"
mkdir "K"$i"/test_k"$i"_200000_1"
done

DATA=`cat files.txt`
i=0
# echo $i
while read line
do
        i=$((i + 1))
        echo "sed -i -e 's|TARGET|"$line"|' "K"$i"/Dockerfile""
done << FILE
$DATA
FILE


