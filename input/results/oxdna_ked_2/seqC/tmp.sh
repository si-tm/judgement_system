DATA=`cat files.txt`
i=0
# echo $i
while read line
do
	i=$((i + 1))
	echo "sed -i -e 's|TARGET|"$line"|' "C"$i"/Dockerfile""
        # echo "sed -i -e"
	# echo "'s|TARGET|"$line"|'" 
	# echo "C"$i"/Dockerfile"
	# sed -i -e "'s|TARGET|"$line"|' " "C"$i"/Dockerfile"
	# sed: -e expression #1, char 1: unknown command: `''
done << FILE
$DATA
FILE

