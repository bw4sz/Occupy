#!/bin/bash 

dns=$(aws ec2 describe-instances --query 'Reservations[*].Instances[*].PublicDnsName' --output text | grep a)

#convert to array
d=($dns)

echo $d

#list of dispersions
disp=(0.5 0.75)

for i in ${!disp[@]}; do
	echo DNS: ${d[$i]} $i Dispersion: ${disp[$i]}
	
	#upload the init file
	scp -i "C:/Users/Ben/Dropbox/Amazon/ec2.pem" init.bash ubuntu@${d[$i]}:~

	#run the file, password needs to be changed
	ssh -i "C:/Users/Ben/Dropbox/Amazon/ec2.pem" ubuntu@${d[$i]} "bash init.bash ${disp[$i]}" 2> out.txt &
	#ssh -i "C:/Users/Ben/Dropbox/Amazon/ec2.pem" ubuntu@${d[$i]} rm -rf Occupy

	# Need to confirm host?
	echo "Yes"
done

