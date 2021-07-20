#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
run "g++ -o selectionM selectionM.cpp `root-config --cflags --glibs` -I. "
