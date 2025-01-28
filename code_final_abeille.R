#author = asyl_scy
#------------------------------------------------------------------------------------------------fonction

auc_roc<-function(the_model,test_data, name_model,suffixe){
  #' fonction qui va gérer la génération/exportation des fichiers liés aux AUC/ROC
  #' @param the_model le modèle
  #' @param test_data les données test
  #' @param name_model le nom du model= soit 'knn', 'svmL', 'svmP', 'svmR'
  #' @param suffixe fait référence au type de test fait (Attention! faut qu'il y ait une direction "/results/suffixe/")|| pour moi suffixe = pyme_300,pyme600 ou pyme_300_600
  #' 


  model_predict <- predict(the_model, newdata = test_data[,-length(test_data)], type = "prob")[, 2]
  
  roc_model <- roc(test_data$experience, model_predict)
  #[,-length(testKNN)]
  # Calcul de l'AUC pour le model
  auc_model <- auc(roc_model)
  print(paste("AUC pour le modèle",name_model,":", auc_model))
  #----export auc

  
  stock_auc<-as.data.frame(matrix(ncol=1,nrow=1))
  colnames(stock_auc)<-"AUC"
  stock_auc[1]<-auc_model
  write.csv2(stock_auc, paste("results/",suffixe,"/auc_",suffixe,"_",name_model,".csv",sep=""),row.names = FALSE)
  
  # Tracer la courbe ROC pour SVM
  png(filename=paste("results/",suffixe,"/ROC_", name_model,"_", suffixe,".png", sep="" ), width=700, height=600, units="px", pointsize=7, res=80)
  plot(roc_model, main = paste("Courbe ROC - ",name_model,sep=""), col = "red")
  dev.off()
  
}

CM_and_co_export<-function(the_model,test_data, name_model,suffixe){
  #' fonction qui va gérer la génération/exportation des fichiers liés aux tables/matrices de confusions
  #' @param the_model le modèle
  #' @param test_data les données test
  #' @param name_model le nom du model= soit 'knn', 'svmL', 'svmP', 'svmR'
  #' @param suffixe fait référence au type de test fait (Attention! faut qu'il y ait une direction "/results/suffixe/")|| pour moi suffixe = pyme_300,pyme600 ou pyme_300_600
  #' 
  
  predfit=predict(the_model,test_data[,-length(test_data)])
  table_model<-table(true=test_data$experience, pred=predfit)
  CM_model<-confusionMatrix(data = factor(predfit), reference = factor(test_data$experience))$byClass[c("Sensitivity","Specificity", "Recall")]
  CM_TOUT<-confusionMatrix(data = factor(predfit), reference = factor(test_data$experience))
  #print(CM_TOUT)
  
  write.csv2(as.data.frame(t(CM_TOUT$overall)), paste("results/",suffixe,"/Confusion_matrix_TOUT_Global", name_model,".csv",sep=""), row.names = TRUE)
  write.csv2(as.data.frame(t(CM_TOUT$byClass)), paste("results/",suffixe,"/Confusion_matrix_TOUT_cl_", name_model,".csv",sep=""),
             row.names = TRUE)
  write.csv2(as.data.frame(table_model), file = paste("results/",suffixe,"/table_", name_model,".csv",sep=""),
             row.names = FALSE)
}


train_probs_AUC<-function(name_model,train_data,test_data,suffixe){
  #' fonction qui va gérer l'entrainement, le test et l'exportation des modèles KNN/SVM
  #' @param the_model le modèle
  #' @param test_data les données test
  #' @param name_model le nom du model= soit 'knn', 'svmL', 'svmP', 'svmR'
  #' @param suffixe fait référence au type de test fait (Attention! faut qu'il y ait une direction "/results/suffixe/")|| pour moi suffixe = pyme_300,pyme600 ou pyme_300_600
  #' 
  the_model<-switch(name_model,
                    "knn"= train(experience~.,data=na.omit(train_data),method="knn",tuneGrid = data.frame(k = c(1,1.5,2.5,2.7)),trControl = trainControl(method = "repeatedcv", number = 3, repeats = 10, classProbs = TRUE)),
                    "svmL" = train(experience~.,data=train_data,method="svmLinear", trControl = trainControl(method = "repeatedcv", number = 3, repeats = 10,classProbs = TRUE ),tuneGrid = expand.grid( C = c(0.1,0.5,1:5))),
                    "svmR"=train(experience~.,data=train_data,method="svmRadial",tuneGrid = expand.grid(sigma = c(0.01,0.1,1), C = c(0.1,0.5,1:5)), trControl = trainControl(method = "repeatedcv", number = 3, repeats = 10,classProbs = TRUE)),
                    "svmP"=train(experience~.,data=train_data,method="svmPoly",tuneGrid = expand.grid(scale = c(0.01,0.1,1), C = c(0.1,0.5,1:5), degree=c(2,3,4)), trControl = trainControl(method = "repeatedcv", number = 3, repeats = 10,classProbs = TRUE)),
                    stop("model pas bon")
  )
  print("train fini!!")
  
  CM_and_co_export(the_model,test_data, name_model,suffixe)
  print("exportation CM fini")
  
  auc_roc(the_model,test_data, name_model,suffixe)
  print("exportation AUC fini")
  
}


#-------------------------------------------------------------------------------------------------code
#-------------import dataset
stock_tout=read.csv2("stock_tout.csv")

#-------------300
#creation dataset
data300=stock_tout[which(stock_tout$experience=="ctrl300"),]#recupère controle 300
data300exp=stock_tout[which(stock_tout$experience=="Pym300exp"),]#recupère exp 300
data300tt=rbind(data300,data300exp)#rassemblement
set.seed(100)
data300tt<-data300tt[-1]#j'enlève visite

#separation train test
idx <- sample(1:nrow(data300tt), as.integer(0.7 * nrow(data300tt)))
train_data_1=na.omit(data300tt[idx,])
test_data_1=na.omit(data300tt[-idx,])

#
set.seed(NULL)
for (i in c("knn", "svmL", "svmR", 'svmP')){
train_probs_AUC(i,train_data_1,test_data_1,"pyme_300")
}

#------------600
#creation dataset
ctrl600=stock_tout[which(stock_tout$experience=="ctrl600"),]
exp600=stock_tout[which(stock_tout$experience=="Pym600exp"),]
TOUT_600=rbind(ctrl600,exp600)
TOUT_600=TOUT_600[-c(7:9)]
TOUT_600=TOUT_600[-1]

set.seed(100)
#separation train test
idx <- sample(1:nrow(TOUT_600), as.integer(0.7 * nrow(TOUT_600)))
train_data_2=na.omit(TOUT_600[idx,])
test_data_2=na.omit(TOUT_600[-idx,])

set.seed(NULL)
for (i in c("knn", "svmL", "svmR", 'svmP')){
  train_probs_AUC(i,train_data_2,test_data_2,"pyme_600")
}

#---------300/600
#creation dataset
stock_tout_2<-stock_tout[,-c(7,8,9)]#"IRT.SD","Cumu.time","Cumu.load"
tout_300_600<-na.omit(stock_tout_2)
tout_300_600<-tout_300_600[-1]#visit

#separation train test
idx <- sample(1:nrow(tout_300_600), as.integer(0.7 * nrow(tout_300_600)))
train_data_3=na.omit(tout_300_600[idx,])
test_data_3=na.omit(tout_300_600[-idx,])

set.seed(NULL)
for (i in c("knn", "svmL", "svmR", 'svmP')){
  train_probs_AUC(i,train_data_3,test_data_3,"pyme_300_600")
}


