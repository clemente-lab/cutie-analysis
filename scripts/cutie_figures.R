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


# scatterplots for figure 2
plist = c("88_NP_50_0.05","88_NP_50_0.9","88_FP_50_0.9","88_FP_50_0.2","88_FN_50_0.9","88_FN_50_0.5")

for (p in plist){
  df = read.csv(file.path(path,paste("cutie_simulations/inputs/",p,".txt", sep='')), header=T, sep='\t', row.names=1)
  df$color <- ifelse(df$X >= 2,"darkorange2", "dodgerblue3")
  X = 'X'
  Y = 'Y'
  ggplot(df, aes_string(x = X, y = Y, color = 'color'))+
    geom_point(size=3, aes(x = X)) + theme_classic()+
    scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
    theme(text = element_text(size=16))+
    labs(x = "X", y = "Y") + theme(legend.position = "none")
  ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,"_2.pdf",sep='')), height=4, width=4, units='in')
}

for (p in plist){
  df = read.csv(file.path(path,paste("cutie_simulations/inputs/",p,".txt", sep='')), header=T, sep='\t', row.names=1)
  df$color <- ifelse(df$X >= 3,"darkorange2", "dodgerblue3")
  X = 'X'
  Y = 'Y'
  ggplot(df, aes_string(x = X, y = Y, color = 'color'))+
    geom_point(size=3, aes(x = X)) + theme_classic()+
    scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
    theme(text = element_text(size=16))+
    labs(x = "X", y = "Y") + theme(legend.position = "none")
  ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,"_3.pdf",sep='')), height=4, width=4, units='in')
}

for (p in plist){
  df = read.csv(file.path(path,paste("cutie_simulations/inputs/",p,".txt", sep='')), header=T, sep='\t', row.names=1)
  df$color <- ifelse(df$X >= 5,"darkorange2", "dodgerblue3")
  X = 'X'
  Y = 'Y'
  ggplot(df, aes_string(x = X, y = Y, color = 'color'))+
    geom_point(size=3, aes(x = X)) + theme_classic()+
    scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
    theme(text = element_text(size=16))+
    labs(x = "X", y = "Y") + theme(legend.position = "none")
  ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,"_5.pdf",sep='')), height=4, width=4, units='in')
}



# figure 2a
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ 
  theme(text = element_text(size=16)) + ggtitle("CUTIE vs Cook's D on FP Simulations") + theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F1_legend.pdf"), height=4, width=5, units='in')

ggplot(df_1e, aes(x=corr_strength, y=indicator, group=Method, color=Method)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=indicator-SE, ymax=indicator+SE), width=.05,
                position=position_dodge(0.05)) + 
  theme_classic()+
  scale_color_manual(values=c('#9BBB59','#4F81BD','#C0504D'))+
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2a.pdf"), height=4, width=4, units='in')



# figure 2b
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")++ 
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2b.pdf"), height=4, width=4, units='in')



# figure 2c
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2c.pdf"), height=4, width=4, units='in')



# S2
df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_50.csv"),
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True") +
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ 
  theme(text = element_text(size=16))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2_legend.pdf"), height=4, width=4, units='in')

df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_pearson_False_NP_50.csv"),
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2i.pdf"), height=4, width=4, units='in')



df_S2= read.csv(file.path(path,"cutie_simulations/outputs/exp20/p_1_spearman_False_NP_50.csv"), 
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16))
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
  labs(x ="Correlation Strength", y = "Proportion Classified as True")+# theme(legend.position = "none")+ 
  scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
  scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
  theme(text = element_text(size=16))
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S2iii.pdf"), height=4, width=4, units='in')

#### Power curves for DFFITS DSR CD

slist = c("dffits",
          "dsr")
plist = c("NP","FN","FP")
path2 = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp'

for (s in slist){
  for (p in plist){
    df_1g = read.csv(file.path(path2,paste("outputs/jobs26/p_1_pearson_",s,"_",p,"_50.csv", sep='')),
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
      scale_color_manual(values=c('#4F81BD','#C0504D','#9BBB59','#8064A2'))+
      labs(x ="Correlation Strength", y = "Proportion Classified as True")+
      scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
      scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
      theme(text = element_text(size=16)) 
    ggsave(file.path(path2,paste("outputs/jobs26/S2_",s,"_",p,".pdf", sep='')), height=4, width=4, units='in')
  }
}

slist = c("cookd")
plist = c("NP","FN","FP")
path2 = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp'

for (s in slist){
  for (p in plist){
    df_1g = read.csv(file.path(path2,paste("outputs/jobs26/p_1_pearson_",s,"_",p,"_50.csv", sep='')),
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
      scale_color_manual(values=c('#9BBB59','#8064A2','#4F81BD','#C0504D'))+
      labs(x ="Correlation Strength", y = "Proportion Classified as True")+# theme(legend.position = "none")+ 
      scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
      scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
      theme(text = element_text(size=16))
    ggsave(file.path(path2,paste("outputs/jobs26/S2_",s,"_",p,".pdf", sep='')), height=4, width=4, units='in')
  }
}


### new S3


slist = c("pearson",
          "spearman",
          "kendall")
plist = c("NP","FN","FP")
path2 = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp'

for (s in slist){
  for (p in plist){
    df_S2 = read.csv(file.path(path2,paste("outputs/jobs27/p_1_",s,"_False_",p,"_50.csv", sep='')),
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
      labs(x ="Correlation Strength", y = "Proportion Classified as True")+ 
      scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
      scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
      theme(text = element_text(size=16)) 
    ggsave(file.path(path2,paste("outputs/jobs27/S3_",s,"_",p,".pdf", sep='')), height=4, width=4, units='in')
  }
}





### new S4

rlist = c("p","r")
slist = c("25",
          "50",
          "100")
plist = c("NP","FN","FP")
path2 = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp'

for (r in rlist){
  for (s in slist){
    for (p in plist){
      df_S2 = read.csv(file.path(path2,paste("outputs/jobs37/",r,"_1_pearson_False_",p,"_",s,".csv", sep='')),
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
        labs(x ="Correlation Strength", y = "Proportion Classified as True")+# theme(legend.position = "none")+ 
        scale_y_continuous(name="Proportion Classified as True", limits=c(0, 1))+
        scale_x_continuous(breaks=pretty(l))+ theme(legend.position = "none")+
        theme(text = element_text(size=16)) 
      ggsave(file.path(path2,paste("outputs/jobs37/S4_",r,"_",s,"_",p,".pdf", sep='')), height=4, width=4, units='in')
    }
  }
}


# figure 3
df_2a = read.csv(file.path(path,"cutie_figures/outputs/exp06/Fig2_rawdata_df.txt"), 
                header=T,  sep="\t", row.names = 1)
df_2a = read.csv(file.path(path,"cutie_figures/outputs/exp06/condensed_results.txt"), 
                 header=T,  sep="\t")#, row.names = 1)

options(scipen=100000)
df = df_2a[which(df_2a$analysis_id=='p_nomc_1_pearson_False_lungc' | df_2a$analysis_id=='p_nomc_1_pearson_True_lungc'), ]
df$analysis_id = c('CUTIE','Cook\'s D')
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=c('CUTIE','Cook\'s D'))
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Method", y = "Number of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16))+
  ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2ai.pdf"), height=3, width=3, units='in')

options(scipen=100000)
df = df_2a[which(df_2a$analysis_id=='p_nomc_1_pearson_False_hdac' | df_2a$analysis_id=='p_nomc_1_pearson_True_hdac'), ]
df$analysis_id = c('CUTIE','Cook\'s D')
df = melt(df, id=c("analysis_id"))
df = df[which(df$variable != 'rsTP'),]
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=c('CUTIE','Cook\'s D'))
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Method", y = "Number of Correlations")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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

df2 = subset(df, u.c._Prevotellaceae < 0.00075)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE))

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2b_legend.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none") + 
  geom_abline(aes(intercept =  -1.721873e-06  , slope =  1.379536e-03  ),
               color = "darkorange2")+
  geom_abline(aes(intercept =  3.265424e-07  , slope =   6.203977e-04   ),
               color = "dodgerblue3")
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
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2b_legend2.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6a.pdf"), height=4, width=4, units='in')

df2 = subset(df, Imports_Unit_Value < 300)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE))

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none") +
  geom_abline(aes(intercept =  -3.07334211, slope = 0.08394777   ),
            color = "darkorange2") +
  geom_abline(aes(intercept =  20.776668, slope = -0.116673   ),
              color = "dodgerblue3")

ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6a_line.pdf"), height=4, width=4, units='in')



ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none") +
  geom_abline(aes(intercept =  -3.07334211, slope = 0.08394777   ),
              color = "darkorange2") 
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6a_bline.pdf"), height=4, width=4, units='in')






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
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6b.pdf"), height=4, width=4, units='in')


df2 = subset(df, Number_of_Dentistry_Personnel < 3e+05)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE))

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none") +
  geom_abline(aes(intercept =  1.130699e+04, slope = -4.882384e-07    ),
              color = "darkorange2") +
  geom_abline(aes(intercept =  6.750003e+03  , slope = 3.713964e-07   ),
              color = "dodgerblue3")

ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6b_line.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'color'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) + theme(legend.position = "none") +
  geom_abline(aes(intercept =  1.130699e+04, slope = -4.882384e-07    ),
              color = "darkorange2") 

ggsave(file.path(path,"cutie_simulations/outputs/exp20/S6b_bline.pdf"), height=4, width=4, units='in')








geom_abline(aes(intercept =  3.265424e-07  , slope =   6.203977e-04   ),
            color = "dodgerblue3")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/F2bi.pdf"), height=4, width=4, units='in')

df2 = subset(df, u.c._Prevotellaceae < 0.00075)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE))




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
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16))+
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
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16))+
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
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16))+
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
  theme(text = element_text(size=16)) #+
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
  theme(text = element_text(size=16))+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "X", y = "Y") 
ggsave(file.path(path,paste("cutie_simulations/outputs/exp20/Fig_",p,".pdf",sep='')), height=4, width=4, units='in')

df = as.data.frame(df_lungc[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)



# fig nat com
df_spatial = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/spatial/df_f1.csv"),
                                    header=T, row.names=1, skip = 0))
i = 41
j = 15
X = 'PIK3CD'
Y = 'CDR1as' # 'CDR1as(ciRS-7)'
df = as.data.frame(df_spatial[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$CDR1as >= 15000,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$CDR1as >= 15000,"Influential according to CUTIE", "Non-influential")

df2 = subset(df, CDR1as <= 15000)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept = 2154, slope = 8.85),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 4077.411372, slope = 1.563407 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/spatial_41_15_FP.pdf"), height=4, width=4, units='in')


# fig nat com
df_spatial = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/spatial/df_f1.csv"),
                                    header=T, row.names=1, skip = 0))
i = 41
j = 20
X = 'PIK3CD'
Y = 'circSLC8A1' # 'CDR1as(ciRS-7)'
df = as.data.frame(df_spatial[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$circSLC8A1 >= 1600,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$circSLC8A1 >= 1600,"Influential according to CUTIE", "Non-influential")

df2 = subset(df, circSLC8A1 <= 1600)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  389.640538, slope = 0.227836),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 261.98589, slope = 0.46703 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/spatial_41_20_FP.pdf"), height=4, width=4, units='in')




# df covid
df_covid = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/df_covid.txt"),
                                    header=T, row.names=1, skip = 0, sep='\t'))
i = 4
j = 6
X = 'GDP.expense.on.Health.care..15....'
Y = 'Total.cases.due.to.COVID.19..16.' # 'CDR1as(ciRS-7)'
df = as.data.frame(df_covid[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$Total.cases.due.to.COVID.19..16. >= 300000,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Total.cases.due.to.COVID.19..16. >= 300000,"Influential according to CUTIE", "Non-influential")

df2 = subset(df, Total.cases.due.to.COVID.19..16. <= 300000)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  2866.239, slope = 7013.486),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 37728.214, slope = 2716.543 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_4.pdf"), height=4, width=4, units='in')


ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  2866.239, slope = 7013.486),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") 

ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_4b.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_4c.pdf"), height=4, width=4, units='in')



# df covid
df_covid = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/df_covid.txt"),
                                  header=T, row.names=1, skip = 0, sep='\t'))
i = 5
j = 6
X = 'Critical.care.beds.Per.Capita.per.100.000.inhabitants.of.the.country...14.'
Y = 'Total.cases.due.to.COVID.19..16.' # 'CDR1as(ciRS-7)'
df = as.data.frame(df_covid[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$cd = cooks.distance(model)
df$color <- ifelse(df$Total.cases.due.to.COVID.19..16. >= 300000,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Total.cases.due.to.COVID.19..16. >= 300000,"Influential according to CUTIE", "Non-influential")

df2 = subset(df, Total.cases.due.to.COVID.19..16. <= 300000)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  2866.239, slope = 7013.486),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 37728.214, slope = 2716.543 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_5.pdf"), height=4, width=4, units='in')

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  2866.239, slope = 7013.486),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") 

ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_5b.pdf"), height=4, width=4, units='in')


ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) 
ggsave(file.path(path,"cutie_exp/outputs/jobs24/covid_6_5c.pdf"), height=4, width=4, units='in')



# df mennonites
df_mennonites = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/mennonites/df_mennonites.txt"),
                                  header=T, row.names=1, skip = 0, sep='\t'))
i = 152
j = 98
X = "TGF.beta1.values" 
Y = "Mouse.IgA2"    # 'CDR1as(ciRS-7)'
df = as.data.frame(df_mennonites[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$Mouse.IgA2 >= 15,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$Mouse.IgA2 >= 15,"Influential according to CUTIE", "Non-influential")
coef(model)

df2 = subset(df, Mouse.IgA2 <= 15)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  1.245398058, slope = 0.002249325),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 0.52548587, slope = 0.00318395 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/mennonites_151_97.pdf"), height=4, width=4, units='in')




# df mennonites
df_mennonites = as.data.frame(read.csv(file.path(path,"cutie_exp/inputs/mennonites/df_mennonites.txt"),
                                       header=T, row.names=1, skip = 0, sep='\t'))
i = 88
j = 65
X = "Cat.IgA1" 
Y = "biftotal"    # 'CDR1as(ciRS-7)'
df = as.data.frame(df_mennonites[,c(i + 1, j + 1)])
colnames(df) <- c(X, Y)
model = lm(unlist(df[Y], use.names=FALSE)~unlist(df[X], use.names=FALSE), data=df)
df$color <- ifelse(df$biftotal >= 8.5e+07,"darkorange2", "dodgerblue3")
df$Status <- ifelse(df$biftotal >= 8.5e+07,"Influential according to CUTIE", "Non-influential")
coef(model)

df2 = subset(df, biftotal <= 8.5e+07)
model = lm(unlist(df2[Y], use.names=FALSE)~unlist(df2[X], use.names=FALSE), data=df2)
coef(model)

ggplot(df, aes_string(x = X, y = Y, color = 'Status'))  + 
  geom_point(size=3)+#, color=df$color)+# shape=1 + 
  theme_classic()+scale_color_manual(values=c("darkorange2", "dodgerblue3"))+
  theme(text = element_text(size=16))+theme(legend.position = "none")+
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = X, y = Y) +
  geom_abline(#data = dfc, 
    aes(intercept =  10512428.4, slope = 409109.5),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "darkorange2") +
  geom_abline(#data = dfc, 
    aes(intercept = 18221970.0, slope = -153429.1 ),
    #              aes(intercept = `(Intercept)`, slope = dfc),  
    color = "dodgerblue3")

ggsave(file.path(path,"cutie_exp/outputs/jobs24/mennonites_87_65.pdf"), height=4, width=4, units='in')


# boxplot for meta analysis



slist = c("Pearson",
          "Spearman",
          "Kendall")
plist = c("FDR","FOR","P")
path2 = '/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp'
for (s in slist){
  # plist=c('P')
  for (p in plist){
    df = read.csv(file.path(path2,paste("outputs/jobs20/",p,"_",s,".csv", sep='')), header=T, row.names=1)
    # Basic box plot
    ggplot(df, aes_string(x="Type", y=p, fill="Type")) + 
      geom_boxplot()+ #fill='gray'
      # labs(x="Datatype", y = p)+#title="Positive Rate across Datatypes",+
      labs(x='',y='') +
      theme_classic() + 
      scale_y_continuous(limits=c(0,1)) + 
      theme(legend.position = "none", text = element_text(size=20)) + 
      geom_jitter(shape=16, position=position_jitter(0.1))
    
    # Change fill colors manually :
    # Continuous colors
    # bp + scale_fill_brewer(palette="Blues") + theme_classic()
    ggsave(file.path(path2,paste("outputs/jobs20/Fig_",p,"_",s,"_2.pdf",sep='')), height=4, width=4.5, units='in')
  }
}



list=c('LungC','Cell','PLOS','Micro') # prop_micro.txt
c('LungTX','LungPT','CRC','HGStool','HGOral','IBD','CA','Statin','Mennonites') # prop_mixed
c('NC','LiverM','LiverF','CovidLong','Spatial','HDAC') # prop_gene
c('WHO','Covid','Airplane','Baseball') #prop_other

# fig barplots
df_s8 = read.csv(file.path(path,"cutie_exp/outputs/jobs25/prop_micro.txt"), 
                 header=T,  sep="\t")#, row.names = 1)
options(scipen=100000)
df = df_s8
labels = rep(c('P','S','K'),length(list))
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=list))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16)) +
  facet_wrap(~ dataset, nrow=1)# +
# ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S88i.pdf"), height=3, width=length(list)*2, units='in')


list=c('LungTX','LungPT','CRC','HGStool','HGOral','IBD','CA','Statin','Mennonites') # prop_mixed
c('NC','LiverM','LiverF','CovidLong','Spatial','HDAC') # prop_gene
c('WHO','Covid','Airplane','Baseball') #prop_other

# fig barplots
df_s8 = read.csv(file.path(path,"cutie_exp/outputs/jobs25/prop_mixed.txt"), 
                 header=T,  sep="\t")#, row.names = 1)
options(scipen=100000)
df = df_s8
labels = rep(c('P','S','K'),length(list))
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=list))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16)) +
  facet_wrap(~ dataset, nrow=1)# +
# ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S88ii.pdf"), height=3, width=length(list)*2, units='in')

list=c('NC','LiverM','LiverF','CovidLong','Spatial','HDAC') # prop_gene
c('WHO','Covid','Airplane','Baseball') #prop_other

# fig barplots
df_s8 = read.csv(file.path(path,"cutie_exp/outputs/jobs25/prop_gene.txt"), 
                 header=T,  sep="\t")#, row.names = 1)
options(scipen=100000)
df = df_s8
labels = rep(c('P','S','K'),length(list))
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=list))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16)) +
  facet_wrap(~ dataset, nrow=1)# +
# ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S88iii.pdf"), height=3, width=length(list)*2, units='in')


list=c('WHO','Covid','Airplane','Baseball') #prop_other

# fig barplots
df_s8 = read.csv(file.path(path,"cutie_exp/outputs/jobs25/prop_other.txt"), 
                 header=T,  sep="\t")#, row.names = 1)
options(scipen=100000)
df = df_s8
labels = rep(c('P','S','K'),length(list))
df$analysis_id = labels
df = melt(df, id=c("analysis_id", "dataset"))
df$Class = factor(df$variable, levels=c('TP', 'FP', 'FN', 'TN'))# 'rsTP', 'FP', 'FN', 'TN'))
df$analysis_id = factor(df$analysis_id, levels=labels)
class_colors = c('#66b3ff',  '#ff9999', '#228B22','#8064A2')# '#ADD8E6', '#ff9999', '#228B22','#8064A2')
df = transform(df, analysis_id=factor(analysis_id,levels=c("P","S","K")))
df = transform(df, dataset=factor(dataset,levels=list))
ggplot(df, aes(fill=Class, y=value, x=analysis_id)) + scale_fill_manual(values=class_colors) +
  geom_bar(position= position_stack(reverse = TRUE), stat="identity")  + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Metric", y = "Proportion")+ scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme(text = element_text(size=16)) +
  facet_wrap(~ dataset, nrow=1)# +
# ggtitle("Microbiome") + theme(plot.title = element_text(hjust = 0.5)) #+ theme(legend.position = "none")
ggsave(file.path(path,"cutie_simulations/outputs/exp20/S88iv.pdf"), height=3, width=length(list)*2, units='in')



### UPSETR

df <- read.csv('~/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs06/nat_Spearman_nomc/data_processing/summary_df_resample_1.txt', sep='\t')
p = df$pvalues
df$q = p.adjust(p, method = "BY")
df <- df[order(df$q),]

library(UpSetR)
input <- c(
  "Pearson"=1466,
  "Spearman"=34,
  "Kendall"=54,
  "Pearson&Spearman"=3,
  "Pearson&Kendall"=7,
  "Spearman&Kendall"=542,
  "Pearson&Spearman&Kendall"=108
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/statin_FP.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
dev.off()

input <- c(
  "Pearson"=1079,
  "Spearman"=63,
  "Kendall"=53,
  "Pearson&Spearman"=3,
  "Pearson&Kendall"=6,
  "Spearman&Kendall"=669,
  "Pearson&Spearman&Kendall"=32
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/statin_FN.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
dev.off()


input <- c(
  "Pearson"=71,
  "Spearman"=11,
  "Kendall"=19,
  "Pearson&Spearman"=0,
  "Pearson&Kendall"=0,
  "Spearman&Kendall"=21,
  "Pearson&Spearman&Kendall"=0
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/ibd_FP.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      empty.intersections=T,
      show.numbers=F)
dev.off()

input <- c(
  "Pearson"=115,
  "Spearman"=11,
  "Kendall"=13,
  "Pearson&Spearman"=0,
  "Pearson&Kendall"=0,
  "Spearman&Kendall"=32,
  "Pearson&Spearman&Kendall"=0
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/ibd_FN.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      empty.intersections=T,
      set_size.show=F,
      show.numbers=F)
dev.off()


input <- c(
  "Pearson"=69,
  "Spearman"=10,
  "Kendall"=10,
  "Pearson&Spearman"=8,
  "Pearson&Kendall"=11,
  "Spearman&Kendall"=30,
  "Pearson&Spearman&Kendall"=20
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/spatial_FP.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
dev.off()


input <- c(
  "Pearson"=70,
  "Spearman"=32,
  "Kendall"=19,
  "Pearson&Spearman"=6,
  "Pearson&Kendall"=8,
  "Spearman&Kendall"=65,
  "Pearson&Spearman&Kendall"=30
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/spatial_FN.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
dev.off()


input <- c(
  "Pearson"=3797,
  "Spearman"=654,
  "Kendall"=701,
  "Pearson&Spearman"=136,
  "Pearson&Kendall"=147,
  "Spearman&Kendall"=1343,
  "Pearson&Spearman&Kendall"=319
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/who_FP.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
dev.off()


input <- c(
  "Pearson"=6239,
  "Spearman"=891,
  "Kendall"=640,
  "Pearson&Spearman"=234,
  "Pearson&Kendall"=158,
  "Spearman&Kendall"=1975,
  "Pearson&Spearman&Kendall"=660
)
pdf(file = "/Users/KevinBu/Desktop/clemente_lab/Projects/cutie_paper/cutie_exp/outputs/jobs23/who_FN.pdf", width = 6, height = 6, onefile = F)
upset(fromExpression(input), 
      order.by = "degree", 
      text.scale = c(2.5,2.5,1.25,1.25,2,1), 
      sets=c('Kendall','Spearman','Pearson'),
      keep.order=T,
      set_size.show=F,
      show.numbers=F)
# show.numbers=F)#, set_size.angles=90, set_size.show=T,set_size.numbers_size = 5)
dev.off()


