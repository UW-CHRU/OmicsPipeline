

## outcome  string column name
## covariates  should be a c() list e.g.  "c('sex','age','bmi')"
## omic_fn File path for omics data.  Should be matrix with omic as row, sample as col
## phenofile File path for phenotype file.  Can be .csv or .dta (Stata)
## sample_id should match the labels on the omics matrix
## subject_id used for multiple observation GEE for grouping
## prefix output file name prefix.  Will output prefix_results.csv and prefix_QQ.png
## winsor FALSE, integer  # integer is the number of sd for winsorizing omic measurements
## impute FALSE , 'zero', 'halfmin' # imputation stategy for omic.  Missing data replaced with zeros, replaced with half the minimum value or not imputed
## transform FALSE, 'invnt', 'std01', 'log2', 'log2_std01'   # Transformations for omic data.  Inverser normal transform ( invnt ), standardizatio mean-0 var-1 ( std01 ), log2 or log2 then standardize ( log2_std01 )
## model 'continuous','binary','survival', 'continuous_gee' 
## ttevent  string column name - only used for survival, otherwise ignored


cat("Omics Pipeline Version v01\n") 
print(sessionInfo())

doOneCox = function(i,mb,pheno,compformula,outcome,estLab='y'){

    pheno$y <- as.numeric(mb[i,]);
    cpheno = pheno[!is.na(pheno$y),]
    NM = NROW(cpheno)
    
    if(NM > 4 &  NROW(cpheno[cpheno[[outcome]]==1,]) >= 5 & var(cpheno$y) > 0){
        
        tryCatch({lres = coef(summary(coxph(compformula,data=cpheno)))
        return( c(row.names(mb)[i], lres[estLab,1], lres[estLab,3], lres[estLab,5],NM))},
        error=function(e) return( c(row.names(mb)[i], NA, NA, NA,-9)))
    }else{
        return( c(row.names(mb)[i], NA, NA, NA,0))
        
    }
}

doOne = function(i,mb,pheno,compformula,model='continuous',estLab='y'){
    pheno$y <- as.numeric(mb[i,]);
    cpheno = pheno[!is.na(pheno$y),]
    
    NM = NROW(cpheno)
    
    if(NM > 4 & var(cpheno$y) > 0){
        if(model == 'binary'){
            tryCatch({lres = coef(summary(glm(compformula,data=cpheno,family='binomial')))
                return( c(row.names(mb)[i], lres[estLab,1], lres[estLab,2], lres[estLab,4],NM))},
                error=function(e) return( c(row.names(mb)[i], NA, NA, NA,-9)))
            
        }else{
            tryCatch({lres = coef(summary(lm(compformula,data=cpheno)))
                return( c(row.names(mb)[i], lres[estLab,1], lres[estLab,2], lres[estLab,4],NM))},
                error=function(e) return( c(row.names(mb)[i], NA, NA, NA,-9)))
        }
    
    }else{
        return( c(row.names(mb)[i], NA, NA, NA,0))
    }   
}


doOneGEE = function(i,mb,pheno,compformula,idtype_cluster=NA,estLab='y'){
    pheno$y <- as.numeric(mb[i,]);
    cpheno = pheno[!is.na(pheno$y),]
    NM = NROW(cpheno)
    environment(compformula) <- environment()
    if(NM > 4 & var(cpheno$y) > 0.001){
      
        tryCatch({result <- gee(compformula, id=as.factor(cpheno[[idtype_cluster]]), data=cpheno, silent = FALSE);
        tbl <- summary(result)$coef;
        
        return( c(row.names(mb)[i], tbl[estLab,1], tbl[estLab,2], 2*pnorm(-abs(tbl[estLab,1])/tbl[estLab,2]),NM))},
               error=function(e) return( c(names(mb)[i], NA, NA, NA,-9)))
    
    }else{
        return( c(names(mb)[i], NA, NA, NA,0))
    }   
}

ln.std <- function(x) {
    x <- log(x, 2)
    x <- (x - mean(x))/sd(x)
    return(x)
}

qqpval <- function(x, main=""){
    x<-sort(-log(x[x>0],10))
    n<-length(x)
    plot(x=qexp(ppoints(n))/log(10), y=x, xlab="Expected", ylab="Observed", main=main )
    abline(0,1,lty=3, col="purple")
}


runPipeline = function(outcome, covariates, omic_fn, phenofile, sample_id = 'sample_id', subject_id = 'subject_id', prefix ='out', interaction=FALSE,model='continuous', ttevent = '', transform=FALSE, winsor=FALSE,  num_cores=10,  impute=FALSE){


    # Read in phenotype 
    pheno_suff = strsplit(phenofile,"\\.")[[1]][-1]
    cat('reading ',phenofile,'  file type: ',pheno_suff,'\n')
    if(pheno_suff == 'dta'){
        library(foreign)
        pheno = read.dta(phenofile,convert.factors=F)
    }else if(pheno_suff == 'csv'){
        pheno = read.csv(phenofile,as.is=T)
    }
    cat('Npheno =',NROW(pheno),'\n')
    
    
    
    mmat = lm(as.formula(paste(outcome, covariates, sep=" ~ ")),data=pheno)
    cat('Design model\n')
    print(summary(mmat))
    if (!is.null(mmat$na.action)) pheno <- pheno[-mmat$na.action,]
    cat('Npheno complete cases=',NROW(pheno),'\n')
    
    
    library(parallel)
    library(data.table)
    
    (num_cores <- detectCores(logical=TRUE))
    options(mc.cores = num_cores)

    cat('reading omic data...\n')
    #samps = read.csv(paste0(datadir,omic_root,'_samples.csv'),as.is=T)
    #info = fread(paste0(datadir,omic_root,'_info.csv'))
    
    omic_suff = strsplit(omic_fn,"\\.")[[1]][-1] 
    
    if(pheno_suff == 'Rda'){
        print('Loading omic data from Rda file')
        data = get(load(omic_fn))
    }else{
        print('Reading omic data from csv file')
        data = read.csv(omic_fn,as.is=T,row.names=1,check.names=F)
    }
    
    print('omic as rows, samples as columns:')
    print(data[1:3,1:3])
    

    
    keep = intersect(pheno[[sample_id]],colnames(data))
    cat('N analysis after omic merge=',length(keep),'\n')
    pheno = pheno[!is.na(pheno[[sample_id]]),]
    pheno = pheno[match(keep,pheno[[sample_id]]),]
    if(model == 'continuous_gee' | model == 'cox_cluster'){
       
        pheno = pheno[order(pheno[[subject_id]]),]
    }
    data = data[,match(pheno[[sample_id]],colnames(data))]

    ## WINSORIZATION OPTIONS -- v0.6 had this after IMPUTE OPTIONS
    winsor = as.numeric(winsor)
    if(winsor > 0){
    	cat('\nWinsorization at', winsor, ' SD\n')
        winz = function(x,nsd=3){
            mm = mean(x,na.rm=T)
            msd = sd(x,na.rm=T)
            top = mm+msd*nsd
            bot = mm-msd*nsd
            
            x[x<bot] = bot
            x[x>top] = top
            x
        }
        data = t(apply(data,1,winz,nsd=winsor))
    }else{
    	cat('\nNO Winsorization \n')
    }
    
    ## IMPUTE OPTIONS
    if(impute == 'zero'){
        cat('\nIMPUTING MISSING VALUES to ZERO\n')
        data[is.na(data)] = 0
    }else if(impute == 'halfmin'){
        cat('\nIMPUTING MISSING VALUES to half the minimum\n')
        data = t(apply(data,1,FUN=function(x){ if(length( x[!is.na(x)]) > 1){ x[is.na(x)] = min(x,na.rm=T)/2}; return(x)}))
    }else{
        cat('\nNo imputation of omics data\n')
    }

    
    ## STANDARDIZATION OPTIONS
    if(transform == 'invnt'){
        cat('\nInverse Normal Transformation\n')
        invnt = function(x){
            qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
        }
        data = t(apply(data,1,invnt))
    }else if(transform == 'std01'){
        cat('\nStandardize 0-1\n')
        data = t(scale(t(data)))
    }else if(transform == 'log2'){
        cat('\nTransform log2\n')
        data <- log(data, 2)
    }else if(transform == 'log2_std01'){
        cat('\nTransform log2 then standardize 0-1\n')
        data = t(apply(as.matrix(data),1,FUN=ln.std))
    }else{
        cat('\nNo Standardization\n')
    }

    
    if(! identical(as.character(pheno[[sample_id]]),as.character(colnames(data)))){
        stop("Something wrong wth sample ordering!!!")
    }
    
    if(interaction != FALSE & model %in% c('binary','continuous')){
    	  print("RUNNING INTERACTION")
        compformula = as.formula(paste(outcome, paste(c(paste0('y*',interaction),covariates), collapse=" + "), sep=" ~ "))
        system.time(result_all <- do.call(rbind, mclapply(setNames(seq_len(nrow(data)), dimnames(data)[[1]] ), doOne,  mb=data, pheno=pheno,compformula, model=model,estLab=paste0('y:',interaction)))) 
    }else if(interaction != FALSE & model == 'survival'){
        print('Running survival interaction...\n')
        library(survival)
        survformula = paste('Surv(',ttevent,',',outcome,')')
        compformula = as.formula(paste(survformula, paste(c(paste0('y*',interaction),covariates), collapse=" + "), sep=" ~ "))
        system.time(result_all <- do.call(rbind, mclapply(setNames(seq_len(nrow(data)), dimnames(data)[[1]] ), doOneCox,  mb=data, pheno=pheno,compformula,outcome=outcome,estLab=paste0('y:',interaction))))
    }else if(model == 'survival'){
        print('Running survival...\n')
        library(survival)
        survformula = paste('Surv(',ttevent,',',outcome,')')
        compformula = as.formula(paste(survformula, paste(c('y',covariates), collapse=" + "), sep=" ~ "))
        system.time(result_all <- do.call(rbind, mclapply(setNames(seq_len(nrow(data)), dimnames(data)[[1]] ), doOneCox,  mb=data, pheno=pheno,compformula,outcome=outcome))) 
    }else if(model == 'continuous_gee'){
        library(gee)
        compformula = as.formula(paste(outcome, paste(c('y',covariates), collapse=" + "), sep=" ~ "))
        doOneGEE(3,mb=data, pheno=pheno,compformula,idtype_cluster=subject_id) 
        system.time(result_all <- do.call(rbind, mclapply(setNames(seq_len(nrow(data)), dimnames(data)[[1]] ), doOneGEE,  mb=data, pheno=pheno,compformula,idtype_cluster=subject_id))) 
        
    }else if(model %in% c('binary','continuous')){
        print(paste0('Running regression - model= ',model,'...'))
        compformula = as.formula(paste(outcome, paste(c('y',covariates), collapse=" + "), sep=" ~ "))
        system.time(result_all <- do.call(rbind, mclapply(setNames(seq_len(nrow(data)), dimnames(data)[[1]] ), doOne,  mb=data, pheno=pheno,compformula, model=model))) 
    }else{
      print(paste0('Model not recognized: ',model))
    }
    
    colnames(result_all) <- c('Name',"beta", "SE", "P", 'N')
    out = as.data.frame(result_all,stringsAsFactors=F)
    out$P = as.numeric(out$P)
    
    write.csv(out,file=paste0(prefix,'_results.csv'),quote=T,row.names=F) 
    out$beta = as.numeric(out$beta)
    out$SE = as.numeric(out$SE)
    lambda50 <- with(out,median(beta^2/SE^2,na.rm=TRUE)/qchisq(0.5, df=1))
    
    png(paste0(prefix,'_QQ.png'))
    if(any(out$P[!is.na(out$P)] ==0)) out[out$P == 0 & !is.na(out$P),]$P = 5e-324 
    qqpval(out$P,main=prefix)
    legend('topleft',paste('lambda ',round(lambda50,2)))
    dev.off()
}
