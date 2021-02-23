# OmicsPipeline
Analysis Tools for Omics Data

### omics_pipeline
The omics pipeline performs single marker regression for omics data. It supports continuous, binary and time to event outcomes. 
The pipeline also supports interaction analyses and GEE for continuous outcomes.
Analyses are can be parallelized across multiple threads.
The pipeline supports some common imputation and transformation methods.

### The main functions is runPipeline()
```R
runPipeline(outcome, covariates, omic_fn, phenofile, sample_id = 'sample_id', subject_id = 'subject_id', prefix ='out', interaction=FALSE,model='continuous', ttevent = '', transform=FALSE, winsor=FALSE,  num_cores=10,  impute=FALSE)
```

* outcome  string column name
* covariates  should be a c() list e.g.  "c('sex','age','bmi')"
* omic_fn File path for omics data.  Should be matrix with omic as row, sample as col
* phenofile File path for phenotype file.  Can be .csv or .dta (Stata)
* sample_id should match the labels on the omics matrix
* subject_id used for multiple observation GEE for grouping
* prefix output file name prefix.  Will output prefix_results.csv and prefix_QQ.png
* winsor FALSE, integer  # integer is the number of sd for winsorizing omic measurements
* impute FALSE , 'zero', 'halfmin' # imputation stategy for omic.  Missing data replaced with zeros, replaced with half the minimum value or not imputed
* transform FALSE, 'invnt', 'std01', 'log2', 'log2_std01'   # Transformations for omic data.  Inverser normal transform ( invnt ), standardizatio mean-0 var-1 ( std01 ), log2 or log2 then standardize ( log2_std01 )
* model 'continuous','binary','survival', 'continuous_gee' 
* ttevent  string column name - only used for survival, otherwise ignored


