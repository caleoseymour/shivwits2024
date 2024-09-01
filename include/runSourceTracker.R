## Methodology source:
#  https://lakarstens.github.io/ControllingContaminants16S/Analyses/ControllingContaminants16S_SourceTrackerPrep.html
source('include/sourcetracker-1.0.1/src/SourceTracker.r')

runSourceTracker = function(st_otus,metadata,outdir,filebase,rarefaction)
{
  
  # extract the source environments and source/sink indices
  train.ix <- which(metadata$ControlType=='negative')
  test.ix <- which(metadata$ControlType=='experimental')
  envs <- metadata$Env
  if(is.element('Description',colnames(metadata))) desc <- metadata$Description

  #skip tuning (takes a long time), can determine if it is worth it later
  alpha1 <- alpha2 <- 0.001

  ## Run SourceTracker
  # train SourceTracker on training data ('source' samples)
  st <- sourcetracker(st_otus[train.ix,], envs[train.ix], rarefaction_depth = rarefaction)
  # predict / estimate source proportions on the test data ('sink' samples)
  results <- predict(st, st_otus[test.ix,], alpha1=alpha1, alpha2=alpha2,full.results = TRUE, rarefaction_depth = rarefaction)

  ## Export results
  # get average of full results across all runs of sourcetracker
  res.mean <- apply(results$full.results,c(2,3,4),mean)

  # Get depth of each sample for relative abundance calculation
  sample.depths <- apply(results$full.results[1,,,,drop=F],4,sum)

  # create directory to store the results
  subdir <- paste(outdir,'full_results',sep='/')
  dir.create(subdir,showWarnings=FALSE, recursive=TRUE)

  # write each environment as a separate file
  for(i in 1:length(results$train.envs)){
    env.name <- results$train.envs[i]
    res.mean.i <- res.mean[i,,]
    # handle the case where there is only one sink sample
    if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)

    # make rows be samples, columns be features
    res.mean.i <- t(res.mean.i)

    # ensure proper names are retained
    colnames(res.mean.i) <- colnames(st_otus)
    rownames(res.mean.i) <- results$samplenames

    filename.fractions <- sprintf('%s/%s_%s_contributions-rabund.txt', subdir, filebase, env.name)
    # calculate and save relative abundance
    res.mean.i.ra <- sweep(res.mean.i.ra,1,sample.depths,'/')
    sink(filename.fractions)
    cat('SampleID\t')
    write.table(res.mean.i.ra,quote=F,sep='\t')
    sink(NULL)

    filename.fractions <- sprintf('%s/%s_%s_contributions-absolute.txt', subdir, filebase, env.name)
    # calculate and save relative abundance
    sink(filename.fractions)
    cat('SampleID\t')
    write.table(res.mean.i,quote=F,sep='\t')
    sink(NULL)
  }

  #generate summary plots
  if(dim(results$draws)[2] > 1) {
    plot.types <- c('pie')
  } else plot.types <- c('pie', 'bar')
  envs<-metadata[rownames(results$proportions),'Env']
  envs<-unlist(envs)
  envs <- as.factor(envs)
  labels = sprintf('%s_%s',envs, rownames(results$proportions))
  plotixs <- sort(as.numeric(envs),index=TRUE)$ix
  for(plot.type in plot.types){
    # plot each environment separately
    for(env in unique(envs)){
      plotixs <- which(envs == env)
      pdf(sprintf('%s/%s_%s_%s.pdf',outdir,filebase,plot.type,env),width=5,height=5)
      plot(results, type=plot.type, labels=labels, include.legend=TRUE, indices=plotixs)
      dev.off()
    }
  }
  return(results)
}