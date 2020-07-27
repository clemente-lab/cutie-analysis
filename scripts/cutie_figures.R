t###
# CUTIE
###
library(ggplot2)
library(reshape)
library(scales)
path = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/'

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      SE = sqrt((sum(x[[col]])/length(x[[col]])) * (1 - (sum(x[[col]])/length(x[[col]])))/length(x[[col]])))
      # sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# figure ?
plist = c("88_NP_50_0.05","88_NP_50_0.9","88_FP_50_0.9","88_FP_50_0.2","88_FN_50_0.9","88_FN_50_0.5")

for (p in plist){
  df = read.csv(file.path(path,paste("cutie_simulations/inputs/",p,".txt", sep='')), header=T, sep='\t', row.names=1)
  ggplot(df, aes(x = X, y = Y, colors = 'dodgerblue3'))+# , fill = Collection_method))  + 
    geom_point(colour = 'cornflowerblue', size=3, aes(x = X)) + 
    theme_classic()+
    theme(text = element_text(size=15))+
    #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
    labs(x = "X", y = "Y") 
  ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,".pdf",sep='')), height=4, width=4, units='in')
}

# figure 1e
df_1e = read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_cookdcompare_FP_50.csv"), 
                   header=T,  sep=",", row.names = 1)

df_1e <- data_summary(df_1e, varname="indicator", 
                    groupnames=c("Method", "corr_strength"))

l = as.numeric(df_1e$corr_strength)
ggplot(df_1e, aes(x=corr_strength, y=indicator, group=Method, color=Method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#9BBB59','#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ #theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("CUTIE vs Cook's D on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F1_legend.pdf"), height=4, width=5, units='in')

ggplot(df_1e, aes(x=corr_strength, y=indicator, group=Method, color=Method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#9BBB59','#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("CUTIE vs Cook's D on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F1e.pdf"), height=4, width=4, units='in')



# figure 1f
df_1f = read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_cookdcompare_FN_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_1f <- data_summary(df_1f, varname="indicator", 
                      groupnames=c("Method", "corr_strength"))

l = as.numeric(df_1f$corr_strength)

ggplot(df_1f, aes(x=corr_strength, y=indicator, group=Method, color=Method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#9BBB59','#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=11)) + ggtitle("CUTIE vs Cook's D on FN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F1f.pdf"), height=4, width=4, units='in')



# figure 1g
df_1g= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_cookdcompare_NP_50.csv"), 
                 header=T,  sep=",", row.names = 1)

df_1g <- data_summary(df_1g, varname="indicator", 
                      groupnames=c("Method", "corr_strength"))

l = as.numeric(df_1g$corr_strength)
ggplot(df_1g, aes(x=corr_strength, y=indicator, group=Method, color=Method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#9BBB59','#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=11)) + ggtitle("CUTIE vs Cook's D on TP/TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F1g.pdf"), height=4, width=4, units='in')




# S2
df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True") + #theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ 
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2_legend.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2i.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_spearman_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Spearman/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2ii.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_kendall_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Kendall/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2iii.pdf"), height=4, width=4, units='in')



# S3
df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_FP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S3i.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_spearman_False_FP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Spearman/CUTIE on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S3ii.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_kendall_False_FP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Kendall/CUTIE on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S3iii.pdf"), height=4, width=4, units='in')




# S4
df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_25.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S4i.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S4ii.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_100.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S4iii.pdf"), height=4, width=4, units='in')



# S5
df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/r_1_pearson_False_NP_25.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S5i.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/r_1_pearson_False_NP_50.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on TP and TN Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S5ii.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/r_1_pearson_False_NP_100.csv"),#p_1_pearson_cookdcompare_NP_50.csv"), 
                header=T,  sep=",", row.names = 1)

df_S2 <- data_summary(df_S2, varname="indicator", 
                      groupnames=c("Significance", "corr_strength"))

l = as.numeric(df_S2$corr_strength)
ggplot(df_S2, aes(x=corr_strength, y=indicator, group=Significance, color=Significance)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion of Correlations classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion of Correlations classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=12)) + ggtitle("Pearson/CUTIE on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
#        axis.text.x = element_text(angle=90, hjust=1)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S5iii.pdf"), height=4, width=4, units='in')









# figure 2a
df_2a = read.csv(file.path(path,"cutie_figures/outputs/exp06/Fig2_rawdata_df.txt"), 
                header=T,  sep="\t", row.names = 1)
df_2a = read.csv(file.path(path,"cutie_figures/outputs/exp06/condensed_results.txt"), 
                 header=T,  sep="\t")#, row.names = 1)

# p_nomc_1_pearson_False_lungc
# p_nomc_1_pearson_True_lungc
options(scipen=100000)
df = df_2a[which(df_2a$analysis_id=='p_nomc_1_pearson_False_lungc' | df_2a$analysis_id=='p_nomc_1_pearson_True_lungc'), ]
df$analysis_id = c('CUTIE','Cook\'s D')
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=c('CUTIE','Cook\'s D'))
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Method", y = "Number of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12))+
  ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2ai.pdf"), height=3, width=3, units='in')

options(scipen=100000)
df = df_2a[which(df_2a$analysis_id=='p_nomc_1_pearson_False_hdac' | df_2a$analysis_id=='p_nomc_1_pearson_True_hdac'), ]
df$analysis_id = c('CUTIE','Cook\'s D')
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=c('CUTIE','Cook\'s D'))
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Method", y = "Number of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12))+
  ggtitle("Gene Expression") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2aii.pdf"), height=3, width=3, units='in')

options(scipen=100000)
df = df_2a[which(df_2a$analysis_id=='p_nomc_1_pearson_False_who' | df_2a$analysis_id=='p_nomc_1_pearson_True_who'), ]
df$analysis_id = c('CUTIE','Cook\'s D')
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=c('CUTIE','Cook\'s D'))
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Method", y = "Number of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12))+
  ggtitle("WHO") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2aiii.pdf"), height=3, width=3, units='in')











# fig 2b lung
df_lungc = as.data.frame(t(read.csv(file.path(path,"cutie_lungc/inputs/otu_table.MSQ34_L6.txt"),
                    header=T, sep='\t', row.names=1, skip = 1)))
i = 162
j = 156
X = 'Bacteroidales_S24.7'
Y = 'u.c._Prevotellaceae'
df = as.data.frame(df_lungc[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$cd >= 1,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$cd >= 1,"Influential according to Cook's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2b_legend.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2bi.pdf"), height=4, width=4, units='in')

i = 317
j = 157
X = 'Shuttleworthia'
Y = 'Prevotella'
df = as.data.frame(df_lungc[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$Shuttleworthia >= 0.0005,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Shuttleworthia >= 0.0005,"Influential according to CUTIE", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color='Status'))+#, aes(color=Status))  + 
  geom_point(size=3)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2b_legend2.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2bii.pdf"), height=4, width=4, units='in')


# fig 2b gene
df_hdac = as.data.frame(t(read.csv(file.path(path,"cutie_MINE/inputs/GSE15222_series_matrix_x1000.txt"),
                      header=T, sep='\t', skip = 62, row.names=1)))

i = 532
j = 87
X = 'GI_12056470.S' # GPRC5A, G protein-coupled receptor class C group 5 member A
Y = 'GI_10834981.S' # IGFBP5, insulin like growth factor binding protein 5 
df = as.data.frame(df_hdac[,c(i+1, j+1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$cd >= 1,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$cd >= 1,"Influential according to Cook's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2biii.pdf"), height=4, width=4, units='in')

i = 47
j = 24
X = 'GI_10337584.S' # DEFA5, defensin alpha 5
Y = 'GI_10190657.S' # putative asosciation with teratozoospermia
df = as.data.frame(df_hdac[,c(i+1, j+1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$GI_10337584.S >= 30,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$GI_10337584.S >= 30,"Influential according to CUTIE", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2biv.pdf"), height=4, width=4, units='in')



# fig 2b who
df_who = read.csv(file.path(path,"cutie_MINE/inputs/WHOfinal.txt"),
                     header=T, sep='\t', row.names=1, skip = 0)

i = 178
j = 46
X = 'Breast_cancer_number_of_female_deaths'
Y = 'Number_of_nursing_and_midwifery_personnel' 
df = as.data.frame(df_who[,c(i+1, j+1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df = df[complete.cases(df), ]
df$cd = cooks.distance(model)
df$color <- ifelse(df$cd >= 1,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$cd >= 1,"Influential according to Cook's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2bv.pdf"), height=4, width=4, units='in')


df_who = read.csv(file.path(path,"cutie_MINE/inputs/WHOfinal.txt"),
                  header=T, sep='\t', row.names=1, skip = 0)

i = 171
j = 79
X = 'Annual_freshwater_withdrawals_total' # GPRC5A, G protein-coupled receptor class C group 5 member A
Y = 'Measles_immunization_coverage_among_one_year_olds_for_lowest_educational_level_of_mother' # IGFBP5, insulin like growth factor binding protein 5 
df = as.data.frame(df_who[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df = df[complete.cases(df), ]
df$cd = cooks.distance(model)
df$color <- ifelse(df$Annual_freshwater_withdrawals_total >= 3000,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Annual_freshwater_withdrawals_total >= 3000,"Influential according to CUTIE", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2bvi.pdf"), height=4, width=4, units='in')










# S6


df_who = read.csv(file.path(path,"cutie_MINE/inputs/WHOfinal.txt"),
                  header=T, sep='\t', row.names=1, skip = 0)

# 242_165.pdf
# 347_43.pdf

i = 242
j = 165
X = 'Imports_Unit_Value' # GPRC5A, G protein-coupled receptor class C group 5 member A
Y = 'Agriculture_Contribution_to_Economy' # IGFBP5, insulin like growth factor binding protein 5 
df = as.data.frame(df_who[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Imports_Unit_Value >= 300,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Imports_Unit_Value >= 300,"Influential and reverse-sign according to CUTIE", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) # + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6_legend.pdf"), height=4, width=8, units='in')

df_who = read.csv(file.path(path,"cutie_MINE/inputs/WHOfinal.txt"),
                  header=T, sep='\t', row.names=1, skip = 0)

# 242_165.pdf
# 347_43.pdf

i = 242
j = 165
X = 'Imports_Unit_Value' # GPRC5A, G protein-coupled receptor class C group 5 member A
Y = 'Agriculture_Contribution_to_Economy' # IGFBP5, insulin like growth factor binding protein 5 
df = as.data.frame(df_who[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Imports_Unit_Value >= 300,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Imports_Unit_Value >= 300,"Influential and reverse-sign according to CUTIE's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6a.pdf"), height=4, width=4, units='in')


options(scipen=100000000000000)
i = 347
j = 43
X = 'Trade_Balance_between_Goods_and_Services' # GPRC5A, G protein-coupled receptor class C group 5 member A
Y = 'Number_of_Dentistry_Personnel' # IGFBP5, insulin like growth factor binding protein 5 
df = as.data.frame(df_who[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Trade_Balance_between_Goods_and_Services < -100000000000,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Trade_Balance_between_Goods_and_Services < -10000000000 ,"Influential and reverse-sign according to CUTIE's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6b.pdf"), height=4, width=4, units='in')



# fig S7


df = read.csv(file.path(path,"cutie_figures/outputs/exp06/p_nomc_1_pearson_False_who/data_processing/counter_samp_resample1.txt"),
                  header=T, sep='\t',  skip = 0)

X = 'Sample_Number'
Y = 'Number_of_False_Positives'
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Number_of_False_Positives >= 600,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Number_of_False_Positives >= 600,"Influential and reverse-sign according to CUTIE", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S7a.pdf"), height=4, width=4, units='in')

df = read.csv(file.path(path,"cutie_figures/outputs/exp06/p_nomc_1_pearson_False_who/data_processing/counter_var_resample1.txt"),
              header=T, sep='\t',  skip = 0)

X = 'Variable_Number'
Y = 'Number_of_False_Positives'
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Number_of_False_Positives >= 70,"dodgerblue3", "darkorange2")
df$Status <- ifelse(df$Number_of_False_Positives >= 70,"Influential and reverse-sign according to CUTIE's D", "Non-influential")

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("dodgerblue3", "darkorange2"))+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S7b.pdf"), height=4, width=4, units='in')




# fig 2c
df_2c = read.csv(file.path(path,"cutie_figures/outputs/exp06/condensed_results_prop_nomc.txt"), 
                 header=T,  sep="\t")#, row.names = 1)

options(scipen=100000)
df = df_2c
# df = df_2c[which(df_2c$analysis_id=='p_nomc_1_pearson_False_lungc' | df_2c$analysis_id=='p_nomc_1_spearman_False_lungc' | df_2c$analysis_id=='p_nomc_1_kendall_False_lungc'), ]
labels = rep(c('P','S','K'),3)
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=c("Microbiome","Gene","WHO")))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12))+
  facet_wrap(~ dataset)# +
  # ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2cI.pdf"), height=3, width=5, units='in')


# fig 2c
df_2c = read.csv(file.path(path,"cutie_figures/outputs/exp06/condensed_results_prop_fdr.txt"), 
                 header=T,  sep="\t")#, row.names = 1)

options(scipen=100000)
df = df_2c
# df = df_2c[which(df_2c$analysis_id=='p_nomc_1_pearson_False_lungc' | df_2c$analysis_id=='p_nomc_1_spearman_False_lungc' | df_2c$analysis_id=='p_nomc_1_kendall_False_lungc'), ]
labels = rep(c('P','S','K'),3)
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=c("Microbiome","Gene","WHO")))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12))+
  facet_wrap(~ dataset)# +
# ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2cII.pdf"), height=3, width=5, units='in')


df = df_2c[which(df_2c$analysis_id=='p_nomc_1_pearson_False_hdac' | df_2c$analysis_id=='p_nomc_1_spearman_False_hdac' | df_2c$analysis_id=='p_nomc_1_kendall_False_hdac'), ]
labels = c('Pearson','Spearman','Kendall')
df$analysis_id = labels
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Gene Expression", y = "Number of Correlations") + scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=12)) #+
  # ggtitle("Gene Expression") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2cii.pdf"), height=3, width=4, units='in')


df = read.csv('/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_simulations/outputs/exp20/24_NP_50_0.7.txt', sep='\t')
model = lm(df$Y ~ df$X)#@, data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$cd >= 1,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$cd >= 1,"Influential according to Cook's D", "Non-influential")




ggplot(df, aes(x = X, y = Y, colors = 'dodgerblue3'))+# , fill = Collection_method))  + 
  geom_point(colour = 'cornflowerblue', size=3, aes(x = X)) + 
  theme_classic()+
  theme(text = element_text(size=15))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "X", y = "Y") 
ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,".pdf",sep='')), height=4, width=4, units='in')

df = as.data.frame(df_lungc[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)





