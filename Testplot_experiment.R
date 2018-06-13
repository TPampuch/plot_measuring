# Here I want to compare 2 different methods of measuring a plot, and their precision
# This is the first script that is still very rough
rm(list = ls())
setwd("D:/PhD stuff/Alaska/plot-measuring/Precision-test")

library(dplyr)
library(geometry)

# read the vertex data first
# this data is written by hand and contains azimut and distance to every point measured from the bases
# b1 for base1 and b2 for base2 measurements
# i use two seperated files
testplot = read.csv("testplot.csv", sep=";")

# 1. calculate plot from both bases seperately and compare the overlap
# 2. scale it with automatic functions and compare to prior results

# 1.:
testplotb1_dist = testplot %>% filter(base == "Base1")
testplotb2_dist = testplot %>% filter(base == "Base2")
testplotwo1_dist = testplot %>% filter(base != "Base1")
testplotwo2_dist = testplot %>% filter(base != "Base2")

testplot_mat = matrix(nrow = length(unique(testplot$tree))+1, ncol = length(unique(testplot$tree))+1)
rownames(testplot_mat) = c("Base1",as.character(unique(testplot$tree)))
colnames(testplot_mat) = c("Base1",as.character(unique(testplot$tree)))
#### now fill it up automatically
for(a in 1:length(testplot$base)){
  base = as.character(testplot$base[a])
  tree = as.character(testplot$tree[a])
  dist = as.numeric(testplot$distance[a])
  testplot_mat[which(colnames(testplot_mat) == tree),
                 which(rownames(testplot_mat) == base)] = dist
  testplot_mat[which(rownames(testplot_mat) == base),
                 which(colnames(testplot_mat) == tree)] = dist
}

# kicks out all values with <2 measurements
counter=1
list = as.numeric()
for(a in 1:length(testplot_mat[,1])){
  if(length(unique(testplot_mat[a,])) < 4){
    list[counter] = a
    counter=counter+1
  }
}

testplot_mat = testplot_mat[-list,-list]


# cleanig data before by setting distance to themself to zero
for(i in 1:length(testplot_mat[1,])){
  testplot_mat[i,i] = 0
}

###### i have to add here the other base stations!!!!!!
for(y in 1:length(testplot_mat[1,])){
  for(x in 1:length(testplot_mat[,1])){
    searcher = testplot_mat[y,x]
    if(identical(searcher,NA_real_)){
      # b = distance to base
      b1 = testplot_mat[which(colnames(testplot_mat) == "Base1"), y]
      b2 = testplot_mat[which(colnames(testplot_mat) == "Base1"), x]
      # a = distance between tree and a common point
      treeSame = names(which(!is.na(testplot_mat[y,]) & !is.na(testplot_mat[x,]) & colnames(testplot_mat) != "Base1"))[1]
      a1 = testplot_mat[which(colnames(testplot_mat) == treeSame), y]
      a2 = testplot_mat[which(colnames(testplot_mat) == treeSame), x]
      # c = distance from treesame to base
      c = testplot_mat[which(colnames(testplot_mat) == "Base1"), which(rownames(testplot_mat) == treeSame)]
      # now it comes to the calculation (remember to add *pi/180 to angles - or in this case - /(pi/180))
      alpha1 = acos((b1^2+c^2-a1^2)/(2*b1*c))/(pi/180)
      alpha2 = acos((b2^2+c^2-a2^2)/(2*b2*c))/(pi/180)
      alphaG = alpha1+alpha2
      aG = sqrt(b1^2+b2^2-2*b1*b2*cos(alphaG*pi/180))
      # now replace the NA by aG
      if(identical(aG,numeric())){
      }else{
        testplot_mat[y,x] = aG
        testplot_mat[x,y] = aG
      }
    }
  }
}

tp1_mat = as.dist((testplot_mat[-c(10,12),-c(10,12)]))

tp1_cmd = cmdscale(tp1_mat, k=2)





#smacof
library(smacof)

tp1_smacof_sp = smacofSym(tp1_mat, type = "mspline", spline.intKnots = 5, spline.degree = 10)
tp1_sp_jackknife = jackknife(tp1_smacof_sp)
sp.perm <- permtest(tp1_smacof_sp, nrep = 1000, verbose = FALSE)

tp1_smacof_int = smacofSym(tp1_mat, type = "interval")
tp1_int_jackknife = jackknife(tp1_smacof_int)
int.perm <- permtest(tp1_smacof_int, nrep = 1000, verbose = FALSE)

par(mfrow=c(2,2))
plot(tp1_smacof_sp, plot.dim = c(1,2))
plot(tp1_smacof_sp, plot.type = "Shepard",
     main = "Shepard Diagram", ylim = c(0.1, 2))
plot(tp1_sp_jackknife)
plot(sp.perm)

plot(tp1_smacof_int, plot.dim = c(1,2))
plot(tp1_smacof_int, plot.type = "Shepard",
     main = "Shepard Diagram", ylim = c(0.1, 2))
plot(tp1_int_jackknife)
plot(int.perm)



##














# Warning: R works with rad and not degrees, this causes problems
# add *pi/180 to degrees to get rad
# here i create a distance-matrix using mostly the "kosinus-satz" that is
# a? = b? + c? - 2 *b*c * cos(alpha)
# i want to try to make everything more easy and automatic next time
dists_b1 = data.frame(base1b1 = c(0, 
                                  b1$dist[3], 
                                  b1$dist[1], 
                                  b1$dist[2]),
                      base2b1 = c(NA, 
                                  0, 
                                  sqrt((b1$dist[1]^2+b1$dist[3]^2)-(2*b1$dist[1]*b1$dist[3]*cos(angles_b1$beta*pi/180))),
                                  sqrt((b1$dist[3]^2+b1$dist[2]^2)-(2*b1$dist[3]*b1$dist[2]*cos(angles_b1$alpha*pi/180)))),
                      p1b1 = c(NA, 
                               NA, 
                               0, 
                               sqrt((b1$dist[2]^2+b1$dist[1]^2)-(2*b1$dist[2]*b1$dist[1]*cos(angles_b1$gamma*pi/180)))),
                      p2b1 = c(NA, 
                               NA, 
                               NA, 
                               0))
row.names(dists_b1) = c("base1b1","base2b1","p1b1","p2b1")

dists_b1d = as.dist(dists_b1)
dists_b1sc = cmdscale(dists_b1d)

plot(dists_b1sc[,1],dists_b1sc[,2])


# do the same for the 2nd base
dists_b2 = data.frame(base2b2 = c(0, 
                                  b2$dist[3], 
                                  b2$dist[1], 
                                  b2$dist[2]),
                      base1b2 = c(NA, 
                                  0, 
                                  sqrt((b2$dist[1]^2+b2$dist[3]^2)-(2*b2$dist[1]*b2$dist[3]*cos(angles_b2$beta*pi/180))),
                                  sqrt((b2$dist[3]^2+b2$dist[2]^2)-(2*b2$dist[3]*b2$dist[2]*cos(angles_b2$alpha*pi/180)))),
                      p1b2 = c(NA, 
                               NA, 
                               0, 
                               sqrt((b2$dist[2]^2+b2$dist[1]^2)-(2*b2$dist[2]*b2$dist[1]*cos(angles_b2$gamma*pi/180)))),
                      p2b2 = c(NA, 
                               NA, 
                               NA, 
                               0))
row.names(dists_b2) = c("base2b2","base1b2","p1b2","p2b2")

dists_b2d = as.dist(dists_b2)
dists_b2sc = cmdscale(dists_b2d)
dists_b2sc[,1] = dists_b2sc[,1]*-1

plot(dists_b2sc[,1],dists_b2sc[,2], xlim = c(-4,4), ylim=c(-2.5,2))

# create a dataframe with coordinates
dists = as.data.frame(rbind(dists_b1sc,dists_b2sc))

dists = data.frame(Tree = c("base1","base2","p1","p2","base2","base1","p1","p2"),
                   x = dists$V1,
                   y = dists$V2,
                   measured_from = c("b1","b1","b1","b1","b2","b2","b2","b2"))

plot(dists$x,dists$y)  


# at last calculate the relative distances again, this time between the same points
differences_vert = data.frame(Tree=dists$Tree, difference_x = 0, difference_y = 0, difference_total = 0)
for(i in 1:length(dists[,1])){
  differences_vert$difference_x[i] = diff(dists$x[which(dists$Tree==differences_vert$Tree[i])])
  differences_vert$difference_y[i] = diff(dists$y[which(dists$Tree==differences_vert$Tree[i])])
}

for(i in 1:length(dists[,1])){
  differences_vert$difference_total[i] = sqrt(differences_vert$difference_x[i]^2+differences_vert$difference_y[i]^2)
}

#get the mean
mean_vert = mean(differences_vert$difference_total)




# now the laser stuff, should be easier
# but they are 3 dimensional and we made a mistake here, therefore the vertical distance will be 1 at least
# we should ignore the 3rd dimension for the moment
laser = read.delim("prezision.txt", skip=6, dec = ",")

plot(laser$X,laser$Y)

cmdLas = cmdscale(dist(laser[,2:4]))
cmdLas[,2] = cmdLas[,2]*-1

laser = as.data.frame(cmdLas)

laser = data.frame(Tree = c("base1","p1","p2","base2","base1","p1","p2"),
                   x = laser$V1,
                   y = laser$V2,
                   measured_from = c("b1","b1","b1","b1","b2","b2","b2"))


# at last calculate the relative distances again, this time between the same points
differences_las = data.frame(Tree=laser$Tree, difference_x = 0, difference_y = 0, difference_total = 0)
for(i in 1:length(laser[,1])){
  
  diffx = diff(laser$x[which(laser$Tree== differences_las$Tree[i])])
  
  if(identical(diffx,numeric())){
    differences_las$difference_x[i] = 0
  }else{differences_las$difference_x[i] = diffx}
  
  diffy = diff(dists$y[which(laser$Tree== differences_las$Tree[i])])
  
  if(identical(diffy,numeric())){
    differences_las$difference_y[i] = 0
  }else{differences_las$difference_y[i] = diffy}
}

for(i in 1:length(laser[,1])){
  differences_las$difference_total[i] = sqrt(differences_las$difference_x[i]^2+differences_las$difference_y[i]^2)
}


#get the mean
mean_las = mean(differences_las$difference_total[c(1:3,5:7)])

#one last comparison
par(mfrow=c(1,2))
plot(laser$x, laser$y, col="red", axes = F, xlab = paste("mean diff. =",round(mean_las, digits = 4)), ylab = "", main = "TruPulse Laser")
abline(v=0,h=0)
plot(dists$x, dists$y, col="blue", axes = F, xlab = paste("mean diff. =",round(mean_vert, digits = 4)), ylab = "", main = "Hagl?f Vertex")
abline(v=0,h=0)

