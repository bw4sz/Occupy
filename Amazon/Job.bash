#!/bin/bash
#cd into occupy, if directory doesn't exist, kill the run.

git clone git@github.com:bw4sz/Occupy.git --depth 1

cd Occupy||sudo halt

#make new branch
#name it the instance ID
iid=$(ec2metadata --instance-id)

git checkout -b $iid

#render script
Rscript -e "rmarkdown::render('Abundance.Rmd')" &> run.txt

#push results
git add --all
git commit -m "ec2 run complete"
git push -u origin $iid

#kill instance
sudo halt
