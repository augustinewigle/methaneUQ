helper_functions.R Helper functions for the paper

uq_model_fitting_code.R Example code to analyse controlled release data for two technologies: QOGI C and airborne NIR HSI. 
                        Code used to implement Application section 4.2.

trial1_anon.csv and trial2_anon.csv - Anonymized controlled release data from the first and second field trials, respectively.

In the csv files, a value of estimate_kgh equal to NA indicates that the technology was unable to quantify emissions although they may have been detected A value of estimate_kgh equal to zero (0) indicates that the technology did not detect emissions in that release

Both of the above situations are removed by the data_prep function in helper_functions.R
