library(CMF)
library(pROC)
library(data.table)

pca_gCMF <- function(fileName, input) {
	conn <- file(fileName,open="r")
	linn <-readLines(conn)
	test_set=strsplit(linn[1], '\t') [[1]][1]
	train_set=strsplit(linn[2], '\t') [[1]][1]

	K=length(linn)-1
	size=length(linn) - 2
	countList=list()
	matrixList=list()
	likelihood=vector()
	likelihood[1]="gaussian"

	inds <- matrix(0,nrow=K,ncol=2)
	inds[1,]=c(as.numeric(strsplit(linn[1], '\t') [[1]][3]),as.numeric(strsplit(linn[1], '\t') [[1]][4]))

	triplets=list()
	X=list()
	F1_F2_F3_SL_binary_all=read.table(paste(input,train_set,sep="/"),sep="\t",header=0)
	F1_F2_F3_SL_binary_all=as.matrix(F1_F2_F3_SL_binary_all)
	countList[[1]]=ncol(F1_F2_F3_SL_binary_all)

	X[[1]]<-matrix(F1_F2_F3_SL_binary_all,nrow=countList[[1]],ncol=countList[[1]])
	triplets[[1]]=matrix_to_triplets(X[[1]])

	triplets_test=list()
	Y=list()

	F3_SL_binary_test=read.table(paste(input,test_set,sep="/"),sep="\t",header=0)
	F3_SL_binary_test=as.matrix(F3_SL_binary_test)

	Y[[1]]<-matrix(F3_SL_binary_test,nrow=countList[[1]],ncol=countList[[1]])
	triplets_test[[1]]=matrix_to_triplets(Y[[1]])

	train=list()
	test=list()
	m1=triplets[[1]]
	m2=triplets_test[[1]]

	m1 <- data.table(m1)
	setkey(m1, 'V1', 'V2')
	m1[,"index1" := .I]
	m2 <- data.table(m2)
	setkey(m2, 'V1', 'V2')
	m2[,"index2" := .I]

	# Join the tables by key #
	m3 <- m1[m2]

	overlap <- m3[is.na(index1)==FALSE & is.na(index2)==FALSE,]

	myIndex=as.matrix(overlap[,4])

	train[[1]]=triplets[[1]][-myIndex,]
	test[[1]]=triplets[[1]][myIndex,]

	myIndex=m2[,V1]

	for (i in 3:length(linn)){

	fields =strsplit(linn[i], '\t') [[1]]
	if(!is.na(fields[1]))
	{
		inds[i-1,]=c(as.numeric(strsplit(linn[i], '\t') [[1]][3]),as.numeric(strsplit(linn[i], '\t') [[1]][4]))
		likelihood[i-1]="gaussian"

		#ADDED ESSENTIALITY
		exprfeatures=read.table(paste(input,fields[1],sep="/"),sep="\t",header=0)
		#replace NA with 0
		exprfeatures[is.na(exprfeatures)] <- 0
		colcount2=ncol(exprfeatures)
		countList[[i-1]]=colcount2
		isPCA = fields[2]
			if (isPCA == 1)
			{
				train1=as.matrix(exprfeatures)
				pca=prcomp(train1)
				pcasummary=summary(pca)$importance[3,]
				pc_keep_len=length(pcasummary[pcasummary<0.96]) #proportion of the variance explained=0.96
				pca_retain=as.matrix(pca$x[,1:pc_keep_len],ncol=pc_keep_len)
				countList[[i-1]]=pc_keep_len

				X[[i-1]]=matrix(as.matrix(pca_retain),nrow=nrow(pca_retain),ncol=ncol(pca_retain))
				triplets[[i-1]]=matrix_to_triplets(X[[i-1]])
				m2<- data.table(triplets[[i-1]])
				test[[i-1]]=triplets[[i-1]][m2[,V1] %in% myIndex,]
				`%not_in%` <- purrr::negate(`%in%`)
				train[[i-1]]=triplets[[i-1]][m2[,V1] %not_in% myIndex,]

			}
			else {
				X[[i-1]]=matrix(as.matrix(exprfeatures),nrow=nrow(exprfeatures),ncol=ncol(exprfeatures))
				triplets[[i-1]]=matrix_to_triplets(X[[i-1]])

				m2<- data.table(triplets[[i-1]])
				test[[i-1]]=triplets[[i-1]][m2[,V1] %in% myIndex,]
				`%not_in%` <- purrr::negate(`%in%`)
				train[[i-1]]=triplets[[i-1]][m2[,V1] %not_in% myIndex,]
			}
		}
	}

	close(conn)
	K=2
	D=as.numeric(countList)
	opts <- getCMFopts()
	opts$iter.max <- 10 # Less iterations for faster computation
	model <- CMF(train,inds,K,likelihood,D,opts=opts)

	outm <- predictCMF(test, model)

	truth=triplets_test[[1]][,3]

	temp_testerr1=model$errors[1,1]
	temp_results=outm$out[[1]]

	num_iter <- 10
	pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
	models <- foreach(i = 1:num_iter, .packages='CMF') %dopar% {
		m <- CMF(train,inds,K,likelihood,D,opts=opts)
		setTxtProgressBar(pb, i)
		m
	}
	
	for (model in models) {
	  newmodelError=model$errors[1,1]
	  if (newmodelError < temp_testerr1)
	  {
	    temp_testerr1=newmodelError
	    out <- predictCMF(test, model)
	    temp_results=out$out[[1]]
	  }
	}
	temp_results
}


