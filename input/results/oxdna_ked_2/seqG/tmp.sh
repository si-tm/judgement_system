for i in 1 2 3 4 5 6 7 8 9 10
do
echo "G"$i
mkdir "G"$i
cp Dockerfile "G"$i"/Dockerfile"
mkdir "G"$i"/test_g"$i"_200000_1"
done

DATA=`cat files.txt`
i=0
# echo $i
while read line
do
        i=$((i + 1))
        echo "sed -i -e 's|TARGET|"$line"|' "G"$i"/Dockerfile""
done << FILE
$DATA
FILE


