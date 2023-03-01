#!/bin/bash

while getopts ":d:a:b:x:y:m:n:" opt; do
    case $opt in
	d) dir="$OPTARG"
	    ;;
	a) irods_1="$OPTARG"
	    ;;
	b) irods_2="$OPTARG"
	    ;;
	x) local_1="$OPTARG"
	    ;;
	y) local_2="$OPTARG"
	    ;;
	m) name_1="$OPTARG"
	    ;;
	n) name_2="$OPTARG"
	    ;;
	\?) echo "Invalid option -$OPTARG" >&2
	    exit 1
	    ;;
    esac

    case $OPTARG in
	-*) echo "Option $opt needs a valid argument"
	    exit 1
	    ;;
    esac
done

raw_1=${dir}/${name_1}
raw_2=${dir}/${name_2}

echo "the raw_1 file is: ${raw_1}"
echo "the raw_2 file is: ${raw_2}"
echo "the local_1 file is: ${local_1}"
echo "the local_2 file is: ${local_2}"

iget -vK ${irods_1} ${dir}
iget -vK ${irods_2} ${dir}

if [[ "${raw_1}" != "${local_1}" ]]; then
    echo "moving ${raw_1} to ${local_1}!"
    mv ${raw_1} ${local_1}
fi

if [[ "${raw_2}" != "${local_2}" ]]; then
    echo "moving ${raw_2} to ${local_2}!"
    mv ${raw_2} ${local_2}
fi
