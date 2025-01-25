# Tensor Analysis with Matricization and Feature Selection
This is the Github page for the code used to generate the data and visualizations showed in paper "Tensor analysis of animal behavior by matricization and feature selection", by Beichen Wang, Jiazhang Cai, Luyang Fang, Ping Ma, Yuk Fai Leung.

## Reproducibility Statement
Random seed were given specific numbers so that results and visualizations can be reproduced. All intermediate output and plots can be found in their correponding folders. Please refer to the "Usage" section for the detailed descriptions on the steps and folder content.

## Table of Contents
- Get Started
  - Packages and libraries
  - Data availability
- Usage
- Contact
- Acknowledgements

## Get Started
### Packages and libraries
Both R and Python are used to generate the results. Below are the packages and libraries used for running the code:
- R: "e1071" (SVM), "naivebayes" (naive bayes), "caret" (k nearest neighbor), "rpart" (decision tree), "randomForest" (random forest), xgboost (extrem gradient boost), "ggplot2" (plotting), "ggrepl" (plotting), "dplyr" (data process), "pROC" (calculate AUROC), "doParallel" (parallel computing), "foreach" (parallel computing), "rstatix" (statistics)
- Python: "pandas" (data frame), "io" (write output), "os" (operating system), "glob" (path names), "sklearn" (machine learning), "numpy" (array), "tensorflow" (neural network), "keras" (neural network), "multiprocessing" (parallel computing).


Please install the packages/libraries following the their instructions prior to running the code.

### Data Availability
The raw data can be accessed by the following Harvard Dataverse link: https://doi.org/10.7910/DVN/B8HBU9

## Usage
Please save the data and the code under the same parent directory. Prior to execution, please double-check the working directory for each code to ensure that the correct files will be loaded. The sequence of exececution follows the order below:

1. Clean and normalize the Data.
  - "Preprocessing"
    - "Preprocessing_normalization_optimized_all_Lighton.R" and "Preprocessing_normalization_optimized_all_Lightoff.R"
    - "combine_control.R"


2. Transform the 3-D tensor data to 2-D tensor data by matricization.
  - "transform_data"
    - "transform_data.R"


3. Plot the distance curve corresponding to Fig 3.
  -  "ave_plot"
    - "curve.R"


4. Conduct the CP decomposition and generate plots corresponding Supp Fig 3.
  - "tensor_decomposition"
    -  "tensor_decomp_on.R" and "tensor_decomp_off.R"


5. Apply the filter methods for feature selection on transformed data.
  - "filter_wilc"
    - "filter_wilc.R"


6. Embedded methods: apply the Embedded methods for feature selection.
  - "embedded_wilc"
    - "tuneRF.1_wilc.R"


7. Tally the number of features for each feature set.
  - "survey"
    - "feature_number_wilc.R"


8. Conduct the 10-fold cross-validation (CV) using different feature sets.
  - "CV_wilc"
    - "\*.R" and "\*. py"


9. Conduct the leave-out validation using the best hyperparameters in CV.
  - "final_test_all_wilc"
    -  "\*.R" and "\*. py"


## Contact
Please contact Beichen "BC" Wang via the following email for any question:
"wang4537 at purdue.edu"

## Acknowledgements
We thank Mengrui Zhang for thoughtful discussions during the early development of this work. Please refer to the publication for additional acknowledgements on funding sources.
