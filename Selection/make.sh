#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
	
run "g++ -O3 -c selection.cpp `root-config --cflags` -I. -o selection.o"
#run "g++ -O3 -c DecayTree.cpp `root-config --cflags` -I. -o DecayTree.o"
#run "g++ -o MiniMaker MiniMaker.cpp `root-config --cflags --glibs` -I. DecayTree.o"
