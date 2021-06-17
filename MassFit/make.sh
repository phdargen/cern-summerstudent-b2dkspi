#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
	
run "g++ `root-config --cflags --libs`  -lRooFit -lRooFitCore -lHtml -lMinuit -c fit.cpp -o fit.o"
