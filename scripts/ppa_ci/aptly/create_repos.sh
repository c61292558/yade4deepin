#!/bin/bash

set -e
for i in bionic bullseye focal jammy bookworm trixie noble
do
    aptly repo create -distribution=$i -component=main yadedaily-$i
done
