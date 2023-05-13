Based on the codes here: http://opendata.cern.ch/record/466

To use it please, do (example for accessing charged particles):

1) After Docker setup, copy the folder "HiForest" inside "~/CMSSW_3_9_2_patch5/src/"

2) Compile the code doing the command: scram b

3) The code that get tracks information is here: "HiForest/HiForestProducer/src/AnalyzerFlow.cc"
   
   The configuration file is here: "HiForest/HiForestProducer/hiforestanalyzer_Flow_cfg.py". 
  
   It executes the C++ code above running the command "cmsRun". Then, do the following (before decide in the configuration file the number of events that will process):
   
   cd HiForest/HiForestProducer/
   
   cmsRun hiforestanalyzer_Flow_cfg.py
   
   The output will be a ROOT file: "HiForestAOD_DATA_FlowAnalysis.root" with 3 ROOT trees. 
   
   The tree in "demo/Flow" has charged particle information and also event information.



