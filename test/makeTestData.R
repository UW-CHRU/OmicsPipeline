set.seed(7)
N = 100 # N samples
P = 30 # N omics

subjects = rep(paste0('subj',1:(N/2)),2)  ## 2 obs per sample
samples = paste0('samp',1:N) 
event =  sample(c(0,1),N,prob=c(.8,.2), replace=TRUE)
tt_event = sample(1:100,N, replace=TRUE)
outcome = abs(rnorm(N,100,15))
x1 = log(abs(rnorm(N,100,100)))
x2 = log(abs(rnorm(N,100,100)))

missing = function(n,prob){
  x = round(abs(rnorm(n,100,15)),2)
  print(x)
  drop = sample(1:n,floor(n*prob))
  print(drop)
  x[drop] = NA
  x
  
}

pheno = data.frame('subj_id' = subjects, 'samp_id' = samples, 'event' = event, 'tt_event' = tt_event , 
           'outcome' = outcome, 'x1' = x1, 'x2' = x2, 
           'm10' = missing(N,.10),
           'm30' = missing(N,.30),
           'm60' = missing(N,.60),
           'm90' = missing(N,.90)
           )
write.csv(pheno,file='test_pheno.csv',row.names=F,quote=F)
pheno = pheno[!duplicated(pheno$subj_id),]

write.csv(pheno,file='test_pheno_singleobs.csv',row.names=F,quote=F)



res = list()
for(i in 1:P){
 res[[i]] =  missing(N,runif(1)) # each omic can have 0-100% missing
}

data = as.matrix(do.call(rbind,res))
colnames(data) = pheno$samp_id
row.names(data) = paste0('omic',1:P)
write.csv(data,file="test_omic_data.csv",quote=F)

