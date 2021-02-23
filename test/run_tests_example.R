source('omics_pipeline_v01.R')
pheno = read.csv('test_pheno_singleobs.csv',as.is=T)
omic = t(read.csv('test_omic_data.csv',as.is=T,row.names=1))
df = merge(pheno,omic[,c('omic1','omic3'),drop=F],by.x='samp_id',by.y='row.names')


## Single Observation Tests
## linear regression
runPipeline(outcome='outcome', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction=FALSE,prefix ='linear_test',
            model='continuous', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)
  
summary(lm(outcome~omic1+x1+m30,data=df))
system('grep -w omic1 linear_test_results.csv')



## logistic regression
runPipeline(outcome='event', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction=FALSE, prefix ='binary_test',
            model='binary', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(glm(event~omic1+x1+m30,data=df,family='binomial'))
system('grep -w omic1 binary_test_results.csv')


## time-to-event 
runPipeline(outcome='event', ttevent = 'tt_event', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction=FALSE, prefix ='ttevent_test',
            model='survival', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(coxph(Surv(tt_event, event) ~ omic3 + x1 + m30,data=df))
system('grep -w omic3 ttevent_test_results.csv')


## linear regression interaction
runPipeline(outcome='outcome', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction='x1',prefix ='linear_test_interaction',
            model='continuous', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(lm(outcome~omic1*x1+m30,data=df))
system('grep -w omic1 linear_test_interaction_results.csv')


## logistic regression interaction
runPipeline(outcome='event', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction='x1',prefix ='binary_test_interaction',
            model='binary', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(glm(event~omic3*x1+m30,data=df,family='binomial'))
system('grep -w omic3 binary_test_interaction_results.csv')




## time-to-event interaction
runPipeline(outcome='event', ttevent = 'tt_event', covariates='m30+x1', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno_singleobs.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction="x1", prefix ='ttevent_test_interaction',
            model='survival', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(coxph(Surv(tt_event, event) ~ omic3*x1 + m30,data=df))
system('grep -w omic3 ttevent_test_interaction_results.csv')




## Multiple Observation Tests

pheno = read.csv('test_pheno.csv',as.is=T)
omic = t(read.csv('test_omic_data.csv',as.is=T,row.names=1))
df = merge(pheno,omic[,c('omic1','omic3'),drop=F],by.x='samp_id',by.y='row.names')


## linear regression

runPipeline(outcome='outcome', covariates='x1+m30', omic_fn = 'test_omic_data.csv', phenofile= 'test_pheno.csv', 
            sample_id = 'samp_id', subject_id = 'subj_id',  interaction=FALSE, prefix ='continuous_gee',
            model='continuous_gee', transform=FALSE, winsor=FALSE,  num_cores=10, impute=FALSE)

summary(gee(outcome~omic3+x1+m30, id=as.factor(df$subj_id), data=df, silent = FALSE))
system('grep -w omic3 continuous_gee_results.csv')
