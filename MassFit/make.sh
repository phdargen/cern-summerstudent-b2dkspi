#!/bin/bash
function run(){
	        echo $1
	        $1
	        echo
}
run "g++ -o fit fit.cpp `root-config --cflags --libs`  -lRooFit -lRooFitCore -lHtml -lMinuit"
