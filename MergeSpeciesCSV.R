############################################################################
# Setup environment and variables
############################################################################
# Variable designiation
# Workspace directory
workspace = "D:/ebird/datadownload/post/mergethese/"
# Input ArcGIS Model csv file
in1 = "cogo"
in2 = "bago"
#in3 = "susc"
bird1 = read.csv(paste(workspace, in1, ".csv", sep=""), header=TRUE)
bird2 = read.csv(paste(workspace, in2, ".csv", sep=""), header=TRUE)
#bird3 = read.csv(paste(workspace, in3, ".csv", sep=""), header=TRUE)
# Join arcData to coverData based on cover_type
mergebird = rbind(bird1, bird2)#, bird3)
write.csv(mergebird, file=paste(workspace, in1, "_", in2, ".csv", sep=""), quote=FALSE, row.names=FALSE)
