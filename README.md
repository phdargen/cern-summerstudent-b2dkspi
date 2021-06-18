# Cern-SummerStudent-B2DKspi

Feasibility study for a measurement of CP violation in B -> D Ks pi decays at LHCb.

# How to get started

- login

ssh username@lxplus.cern.ch -X -Y -o ServerAliveInterval=60

- go to your work directory

cd /afs/cern.ch/work/u/username

- clone repo

git clone ssh://git@gitlab.cern.ch:7999/phdargen/cern-summerstudent-b2dkspi.git

cd cern-summerstudent-b2dkspi

- setup ROOT

/bin/bash

source setup_root.sh

- start interactive ROOT session and look at data

root -l

TFile *_file0 = TFile::Open("/eos/lhcb/user/p/phdargen/summerStudents21/Stripped/Data_B2DKspi_DD_11.root")

new TBrowser

t = (TTree*) _file0->Get("DecayTree")

t->Draw("B_DTF_MM")

- implement selection

cd Selection

source make.sh

./selection

- setup TMVA

cd myTMVA

source setup.sh ../TMVA-v4.2.0/

- train classifier

root -l

.L TMVAClassification.cpp

TMVAClassification( "BDTG", "MC" , "B2DKspi", "all",  "LL" )

- apply classifier

.L TMVAClassificationApplication.cpp+

applyToAll()

# Useful stuff

Git tutorial

https://hsf-training.github.io/analysis-essentials/git/README.html

Root tutorial

https://root.cern/primer/#histograms-in-root

About the B -> D Ks pi decay

https://arxiv.org/pdf/hep-ph/0605129v2.pdf

https://arxiv.org/pdf/0712.3469v1.pdf

TMVA

https://root.cern.ch/download/doc/tmva/TMVAUsersGuide.pdf

RooFit

https://root.cern.ch/download/doc/RooFit_Users_Manual_2.91-33.pdf

Random LHCb resources 

http://cds.cern.ch/record/2702424/files/CERN-THESIS-2019-216.pdf

https://www.sciencedirect.com/science/article/pii/S2405601415006811?via%3Dihub

# Meetings

Sign up to B2OC mailing list:  lhcb-phys-cp-measurements-trees@cern.ch

https://indico.cern.ch/category/6734/overview?period=week

B2OC-TD Tue 14.00 (bi-weekly)

B2OC Amplitude Analyses Thu 10.00 (bi-weekly)

B2OC Thu 14.00

# How to get help

- Mattermost channels:

https://mattermost.web.cern.ch/

- Mailing lists: 

https://groups.cern.ch/directory/Pages/categoryresults.aspx?Column=Experiment&ColumnDisplayName=Experiment&Value=LHCb
