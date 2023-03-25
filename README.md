# Reproducible synthetic faces for research in forensic anthropology

Source code to generate all the analyzes in the article "The use of reproducible synthetic faces for research in forensic anthropology".


This paper presents reproducible methods for generating synthetic face databases drawing on forensic anthropology expertise both for the design of these methods and for proposing criteria that allow the scrutiny of these technologies in a transparent and accessible way.    

From this forensic anthropology perspective, we propose measures of global and specific utility of the synthetic databases generated to mitigate diversity. The global utility of synthetic databases depends on whether they preserve similar morphological diversity to that observed in the original database, while their specific utility depends on whether they preserve a similar type of morphological shape variation to that in the original database.    

Applying the following workflows to the set of facial images of East Asian females found in **./Data/** folder will deliver the results presented in the article. These images are part of the Multi-Racial Mega-Resolution database (Strohminger et al. 2016) licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License, and can be found at: https://osf.io/skbq2/    


Both the global (Global_utility_safety) and the specific (Specific_utility) workflows contain three main files:
- **./R/Functions.R** contains all the custom functions used for the evaluations.    
- **./_targets.R** integrates all the steps to follow for each of the evaluations (see workflow.png).      
- **./tar_exe.R** runs all the steps to follow for each of the evaluations.    

To run each workflow execute the **tar_exe.R** file.    

#### References    
Strohminger N, Gray K, Chituc V, Heffner J, Schein C, Heagins TB. The MR2: A multi-racial, mega-resolution database of facial stimuli. Behav Res. 2016 Sep;48(3):1197–204.

