#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
run "g++ -o selection selection.cpp `root-config --cflags --glibs` -I. "
