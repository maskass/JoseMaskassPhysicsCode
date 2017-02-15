#! /bin/sh

if [ -s old.list ]  
then  
    echo "old.list file exists."
else  
    echo "Creating dummy old.list file."
    touch old.list
fi  
    

while true
do
    ls dst_ssh/*.dat > new.list
    diff --unchanged-line-format="" --old-line-format="" old.list new.list > ascii.list
    
    mv new.list old.list

    if [ -s ascii.list ]  
    then  
	echo "---------------------------------------- "
	echo "CREATING ROOT FILES FROM NEW DST FILES"
	./do_dst2root.sh ascii.list
	echo " "
	echo "NEW ROOT FILES CREATED"
	echo "---------------------------------------- "
    else  
	echo "---------------------------------------- "
	echo " "
	echo "NO NEW DST FILES PRODUCED IN THE LAST 5 MINUTES"
	echo " "
	echo "---------------------------------------- "
    fi  
    
    sleep 300
     
done

exit 0
