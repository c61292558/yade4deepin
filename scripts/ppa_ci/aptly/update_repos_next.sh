#!/bin/bash

# All the next uploads of the packages

PATHDEB="/home/anton/deb/"

set -e

# check, if $PATHDEB/OLD_ID exists and if not, create it
if [ ! -f ${PATHDEB}/OLD_ID ]
then
	touch ${PATHDEB}/OLD_ID
	echo "0" > ${PATHDEB}/OLD_ID
fi

OLD_ID=$(cat ${PATHDEB}/OLD_ID)

NEW_ID=$(curl -L "https://gitlab.com/api/v4/projects/10133144/repository/commits/master" | jq --raw-output '.short_id')
if [ "$OLD_ID" != "$NEW_ID" ]
then
	for i in bionic bullseye focal jammy bookworm trixie noble
	do
		cd ${PATHDEB}
		wget https://gitlab.com/api/v4/projects/10133144/jobs/artifacts/master/download?job=deb_$i -O yade.zip
		unzip yade.zip
		aptly repo remove yadedaily-$i 'yadedaily'
		aptly repo remove yadedaily-$i 'libyadedaily'
		aptly repo remove yadedaily-$i 'yadedaily-doc'
		aptly repo remove yadedaily-$i 'python3-yadedaily'
		aptly repo add yadedaily-$i deb/*.deb
		aptly repo add yadedaily-$i deb/*.dsc
		rm -rf ${PATHDEB}/*
		aptly publish update $i filesystem:yadedaily:
		aptly db cleanup
	done
	echo $NEW_ID > ${PATHDEB}/OLD_ID
else
	echo "No new Yade version is available, exiting"
	exit
fi
