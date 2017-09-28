# LostLepton_avgTF
Estimation of LostLepton BG using an average TransferFactor approach

- Optimized for CMSSW_7_6_0 or later
## How to set up the project
1. Check out code:

  ```
  git clone https://github.com/skurz/LostLepton_avgTF.git
  ```

## Run Package 

- Step-by-step walkthrough (details on how to recalculate the systematics and other adjustments like the definitions of SearchBins, see step 3):

### 1) Create the SFs for Signal- and Control-Region. The code has to be run twice since for the SR, the isotrack veto is applied, but not for the CR:

a) Controlregion

	1. Check SFMaker.h if
	  
	    ```
	    const int useDeltaPhiCut = 1			// -1 for prediction in low delta Phi region
	    const bool includeIsotrkVeto = false;	// true: needed for SR, false: needed for CR
	    ```

	2. Check path to Skims in SFMaker.h (The Skims are necessary since the weights for b-tag reweighting are included.):

	    ```
	    const string path_toSkims("/nfs/dust/cms/user/kurzsimo/LostLepton/skims_v12/SLe/tree_");
	    ```

	3. Check source files (full nTuples, no Skims) in

	    ```
	    MakeSFs.C  // has to be run on full Trees
	    ```

	4. Since you want to produce the SFs for the SR, set the name of the output files in MakeSFs.C to

		```
	    Effchain[i]->Process("SFMaker", TString::Format("SFCR_%d.root",i));
	    ```

	5. Run code:
	  
	    ```
	    root -l MakeSFs.C+
	    ```

	This should produce 4 files with SFs (independent SFs for ttbar, wjets, singlet, other)

b) Signalregion

	1. Check SFMaker.h if
	  
	    ```
	    const int useDeltaPhiCut = 1 // -1 for prediction in low delta Phi region
	    const bool includeIsotrkVeto = true;  // true: needed for SR, false: needed for CR
	    ```

	2., 3. same as in a)

	4. Since you want to produce the SFs for the CR, set the name of the output files in MakeSFs.C to

		```
	    Effchain[i]->Process("SFMaker", TString::Format("SFSR_%d.root",i));
	    ```

	5. same as in a)

### 2) Create the TFs (after completing step 1):
   
1. Check TFMaker.h if
  
    ```
    const int useDeltaPhiCut = 1  // -1 for prediction in low delta Phi region
    ```
    
2. It might come in handy to set the luminosity if you want to have a look at some absolute numbers; the TFs will obviously stay the same:
    
    ```
    const double scaleFactorWeight = 35862.351;
    ```
    
3. Right now the filename of the trees has to follow a convention to make sure that the correct SF is picked (for ttbar, wjets,...). See TFMaker.C:

    ```
    if((std::string(currFileName.Data()).find(std::string("TTJets"))) != std::string::npos){
        SFCR_histFile_path = "SFCR_0.root";
        SFSR_histFile_path = "SFSR_0.root";
    }else if((std::string(currFileName.Data()).find(std::string("WJetsToLNu"))) != std::string::npos){
        SFCR_histFile_path = "SFCR_1.root";
        SFSR_histFile_path = "SFSR_1.root";
    }else if((std::string(currFileName.Data()).find(std::string("ST_"))) != std::string::npos){
        SFCR_histFile_path = "SFCR_2.root";
        SFSR_histFile_path = "SFSR_2.root";
    }else{
        SFCR_histFile_path = "SFCR_3.root";
        SFSR_histFile_path = "SFSR_3.root";
    }
    ```

4. Check source files (full nTuples, no Skims) in

    ```
    MakeTF.C  // has to be run on full Trees
    ```
    
5. Run code:
  
    ```
    root -l MakeTF.C+
    ```
  
### 3) Create the Prediction (after completing step 1). Can be done based on MC or data:
   
a) MC

	1. Check Prediction.h if
	  
	    ```
	    const int useDeltaPhiCut = 1  // -1 for prediction in low delta Phi region
	    ```
	    
	2. Check Prediction.h if
	    
	    ```
	    const bool runOnData = false;
	    const bool runOnStandardModelMC = true;
	    const bool runOnSignalMC = false; 
	    ```
	    
	3. Check source files (Trees/Skims) in

	    ```
	    MakePrediction_Data.C  // can be run on either full Trees or Skims
	    ```
	    
	4. Run code:
	  
	    ```
	    root -l MakePrediction_separate.C+					// Do the prediction based on MC
	    hadd Prediction.root Prediction_separate/*.root		// hadd the histograms
	    ```

b) Data

	1. same as in a)
	    
	2. Check Prediction.h if
	    
	    ```
	    const bool runOnData = true;
	    const bool runOnStandardModelMC = false;
	    const bool runOnSignalMC = false; 
	    ```
	    
	3. Check source files (Trees/Skims) in

	    ```
	    MakePrediction_Data.C  // can be run on either full Trees or Skims
	    ```
	    
	4. Run code:
	  
	    ```
	    root -l MakePrediction_Data.C+  // Do the prediction based on Data
	    ```
	  
### 3) Produce final output 'LLPrediction.root' for e.g. integration etc.:

For the avg. TF method this histogram is already produced in the prediction step (Prediction.root). However uncertainties, plot labels, etc missing. Can be copied from event-by-event approach.
  
### 4) Produce data/MC comparison for paper ('LLPrediction.root' has to exist):

1. Check 'Plot_searchBin_comparison.C' if:
   
    ```
    void Plot_searchBin_comparison(string option="", int pull=0){...  // use option="QCD" to show in QCD binning -> needs LLPrediction.root in QCD binning!    
    bool doDataVsMC = false;                                    // true if you want to show data/MC instead
    ```

    Some additional changes might be necessary if the SearchBins were changed (separation lines in plot etc)
  
2. Run code:
  
    ```
    root -l Plot_searchBin_comparison.C
    ```

### 5) Do signal contamination studies ('Efficiencies.root' has to exist). 

Should already be implemented, however never tested it. In case it doesn't work, or if instructions are needed, have a look at the event by event approach.

### 6) Modifications: SearchBins, recalculation of uncertainties etc

SearchBins can be modified in SearchBins.h. However, none of the uncertainties are implemented yet. The basic idea is the following:

1. The necessary histograms can be generated as in the event-by-event approach with the existing scrips, e.g. uncertainty on MTW cut efficiency, or acceptance.

2. The uncertainty on the final prediction coming from every single one of those effects can be done in a similar way than the application of the SF.

	Example: Acceptance
	Vary the acceptance efficiency up/down by the values that are stored in the histograms that you generated in the step before. This is done like in the very first step where the SFs for the Signal- and Controlregion are produced. (If the number of leptons that are within the acceptance increases, the number of leptons that are out-of-acceptance has to go down by the exact amount so the overall number of single lepton events in the MC sample stays constant).
	In the end, you get two different factors for signal- and controlregion, that can then be applied, and a new TF is calculated. Comparing this to the nomial TF directly leads to the uncertainty on the central prediction coming from the variation of the acceptance.

	Generally speaking, this means that the three steps (MakeSFs, MakeTF, MakePrediction) have to be repeated for every single uncertainty.