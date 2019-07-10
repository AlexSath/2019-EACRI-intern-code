#!/bin/bash
#The goal of this script is to organize all of the .gz files created by the
#fingerprint script based on the file extension
pathOfFingerprints=$1

outnFiles=()
outFiles=()
outCloseFiles=()

createDirIfNecessary () {
  if [ ! -x "$pathOfFingerprints/$1" ]
  then echo "Folder $1 not found, creating $1..."
    mkdir $pathOfFingerprints/$1
  else echo "Folder $1 exists; pass..."
  fi;
}

searchAndMoveFiles () {
  for file in "$pathOfFingerprints/../*.$1.gz";
  do
    mv $file "$pathOfFingerprints/$1"
  done
}


createDirIfNecessary outn
echo "Moving .outn.gz files..."
searchAndMoveFiles outn
createDirIfNecessary out
echo "Moving .out.gz files..."
searchAndMoveFiles out
createDirIfNecessary out.close
echo "Moving .out.close.gz files..."
searchAndMoveFiles out.close

exit
