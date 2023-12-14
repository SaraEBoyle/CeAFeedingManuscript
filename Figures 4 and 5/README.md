# Code used in the creation of figures 4 and 5. 
See https://github.com/SaraEBoyle/MetabolicCageAnalysis

# Metabolic Cage Analysis
- Takes a .csv file from the CLAMS metabolic cage monitoring software (CLAX) and cleans, averages, and plots group data
- Allows users to specify whether cages contain experimental or control mice.
- Allows users to screen for and exclude system artifacts
- Consolidates and graphs key data by group, including respiratory exchange ratio (RER), oxygen, movement, and food and water intake.
- Plots experimental vs. control groups and allows customizable smoothing and time window selection.
- Saves .mat matrix of all individual ("full_graphing_variables.mat") and averaged ("full_metabolic_matrix.mat") data
- Also saves parameters used in session ("RERParameters.mat")

1. Export a .csv of the metabolic data from CLAX software (click data, check "all", click export and save as "all.csv")
2. Open "full_metabolic_analysis.m"
3. Follow the intructions in the comments to specify your desired parameters. You can select certain cages to exclude, the time frame to examine, and downsampling and smoothing options.
4. Run "full_metabolic_analysis.m" from the same folder as your "all.csv" experiment information.
5. On the first run, you will have the option to scan your data for artifacts from the system. The CLAMS metabolic cages often are influenced by electrical noise that will be visible across all cage sensors. This program will display all cage RER data in discrete windows. If you see artifacts that are present across all cages, click to the left, then the right boundary of where you see the artifacts. Selected artifacts will be replaced with average values of the surrounding data. If you see no artifacts, press enter to continue.
![Select artifacts fig](https://github.com/SaraEBoyle/MetabolicCageAnalysis/assets/83416542/03f7f437-3632-438e-a042-daf816f7e662)

6. The program will then graph individual mouse data. Periods during the dark cycle (when mice are most active) are represented as well. Press enter to continue.
![Individual RERs](https://github.com/SaraEBoyle/MetabolicCageAnalysis/assets/83416542/84a52afb-f3b0-4e6c-9332-0b04f55874c3)
 ![food intake](https://github.com/SaraEBoyle/MetabolicCageAnalysis/assets/83416542/a1ae0dc4-5712-4eca-87e2-d29e1ee2bcc6)
  
7. Then it will graph group average data, saving all figures and extracted data. The program calculates group and individual averages from the data after removing artifacts, and also splits the data by light/dark cycle. All data is stored in full_metabolic_matrix.mat, and the data used for graphing is stored in full_graphing_variables.mat, for convenient import into data visualization software.
 ![group rer figure](https://github.com/SaraEBoyle/MetabolicCageAnalysis/assets/83416542/b0ab40ec-b43e-48e0-8d94-ff904d486f36)
  
<img width="1022" alt="Screen Shot 2023-06-29 at 2 52 55 PM" src="https://github.com/SaraEBoyle/MetabolicCageAnalysis/assets/83416542/96521ee5-b366-4685-afb8-3065e1ac8df4">
