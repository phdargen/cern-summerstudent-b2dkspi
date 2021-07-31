This file contains steps in issue6.

B_2Gauss_CB+Cheby:
	fitting with double Gaussian as signal & CB + 4-ordered Chebychev Polynomial as bkg.
	chi2 = 0.828441.
	fitB_CB_2Gauss_Cheby.cpp: code
	fitB_CB_2Gauss_Cheby.txt: output
	fitB_CB_2Gauss_Cheby.eps: plot

B_Gauss_CB+3-Cheby:
	fitting with a Gaussian as signal & CB + 3-ordered Chebychev Polynomial as bkg.
	chi2 = 1.04889
	fitB_CB_Gauss_Cheby.cpp: code
	fitB_CB_Gauss_Cheby.txt: output
	fitB_CB_Gauss_Cheby.eps: plot

B_Gauss_CB+4-Cheby:
	fitting with a Gaussian as signal & CB + 4-ordered Chebychev Polynomial as bkg.
	chi2 = 1.02132
	fitB_CB_Gauss_Cheby.cpp: code
	fitB_CB_Gauss_Cheby.txt: output
	fitB_CB_Gauss_Cheby.eps: plot

B_Gauss_Exp:
	fitting with a Gaussian as signal & Exponential as bkg.
	chi2 = 7.0232
	fitB_Gauss_Exp.cpp: code
	fitB_Gauss_Exp.txt: output
	fitB_Gauss_Exp.eps: plot

B_Gauss_Gauss+4-Cheby:
	fitting with a Gaussian as signal & Gaussian + 4-ordered Chebychev Polynomial as bkg.
	chi2 = 1.23247
	fitB_2Gauss_Cheby.cpp: code
	fitB_2Gauss_Cheby.txt: output
	fitB_2Gauss_Cheby.eps: plot

B_Gauss_Gauss+exp:
	fitting with a Gaussian as signal & Gaussian + Exponential as bkg.
	chi2 = 1.52509
	fitB_2Gauss_Exp.cpp: code
	fitB_2Gauss_Exp.txt: output
	fitB_2Gauss_Exp.eps: plot

B_MCFit:
	Fitting signal with triple Gaussian.
	chi2 = 0.646647
	fitB_MC.cpp: code
	fitB_MC.txt: output
	fitB_MC.eps: plot
