rm(list = ls())
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)


# now to run the simulation
n = 1000 # number of fish
reps=100

x_bins=matrix(data=NA,nrow=reps,ncol=10)
y_bins=matrix(data=NA,nrow=reps,ncol=10)
z_bins=matrix(data=NA,nrow=reps,ncol=10)



for (i in 1:reps){
  # now we add some random fish
  xf=runif(n)*10
  yf=runif(n)*10
  zf=runif(n)*10
  counter=1
  # for(j in 1:10){
  #   for (k in 1:10){
  #     for (l in 1:10){
  #
  #       bins[counter]=length(which(xf>j-1 & xf<=j & yf>k-1 & yf<=k & zf>l-1 & zf<=l))
  #       counter=counter+1
  #     }
  #   }
  # }
  for(j in 1:10){
    x_bins[i,j]=length(which(xf>j-1 & xf<=j))
    y_bins[i,j]=length(which(yf>j-1 & yf<=j))
    z_bins[i,j]=length(which(zf>j-1 & zf<=j))

  }

}
dfx <- as.data.frame(x_bins)
colnames(dfx) <- c("01","02","03","04","05","06","07","08","09","10")
dfx_long <- dfx %>%
  pivot_longer(cols = everything(), names_to = "bins", values_to = "reps")

dfy <- as.data.frame(y_bins)
colnames(dfy) <- c("01","02","03","04","05","06","07","08","09","10")
dfy_long <- dfy %>%
  pivot_longer(cols = everything(), names_to = "bins", values_to = "reps")

dfz <- as.data.frame(z_bins)
colnames(dfz) <- c("01","02","03","04","05","06","07","08","09","10")
dfz_long <- dfz %>%
  pivot_longer(cols = everything(), names_to = "bins", values_to = "reps")

# 3. (Optional) Assign meaningful column names



p1 <- ggplot(data=dfx_long, aes(x=bins, y=reps)) +
  geom_boxplot()+labs(title = "X dimension", x = "bin", y = "target counts") +
  theme_bw()
p2 <- ggplot(data=dfy_long, aes(x=bins, y=reps)) +
  geom_boxplot()+labs(title = "Y dimension", x = "bin", y = "target counts") +
  theme_bw()
p3 <- ggplot(data=dfz_long, aes(x=bins, y=reps)) +
  geom_boxplot()+labs(title = "Z dimension", x = "bin", y = "target counts") +
  theme_bw()


grid.arrange(p1, p2, p3, ncol = 1, nrow = 3)
p1 + p2 + p3 + plot_layout(ncol = 1) # Stacks plots vertically
