
library(tidyr)
library(dplyr)
library(tibble)
library(Biobase)
library(ALL)
library(hgu95av2.db)
library(ggplot2)

data(ALL)
## data
ex_data <- exprs(ALL)[1:30,]
ph_data <- pData(ALL)[,c("cod", "sex", "BT")]

## remove missing | duplicated genes
ph_data <- ph_data[complete.cases(ph_data),]
feature_names <- rownames(ex_data)
gene_names <- unlist(as.list(hgu95av2SYMBOL[feature_names]))
idx <- which(is.na(gene_names) | duplicated(gene_names))
ex_data <- as.data.frame(ex_data[-idx,])
rownames(ex_data) <- gene_names[-idx]
ex_data[1:3,1:3]


ex_data_mlt <- ex_data %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  mutate(bt = ph_data[name,"BT"])

### boxplot
ex_data_mlt %>% 
  group_by(rowname) %>% 
  ggplot(aes(x=bt, y=value, group=bt)) +
  facet_wrap(~rowname, ncol=9, scales="free") +
  geom_boxplot() +
  theme(
    axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust=0.5)
  )

## NA
ph_data$BT
tmp <- ex_data %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) 
ph_data[tmp$name,]$BT


ex_data_mlt <- ex_data %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  mutate(bt = ph_data[name,"BT"]) %>% 
  drop_na()
ex_data_mlt %>% complete.cases()


### 평균 비교
ex_summary <- ex_data_mlt %>% 
  group_by(bt, rowname) %>% 
  summarize(m=mean(value))


ggplot(ex_summary, aes(x=bt, y=m, group=bt)) +
  facet_wrap(~rowname, ncol=9, scales="free") +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust=0.5)
  )

### scale 
ggplot(ex_summary, aes(x=bt, y=m, group=bt)) +
  facet_wrap(~rowname, ncol=9) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust=0.5)
  )

### testing (t-test)
ex_data_mlt <- ex_data %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  mutate(bt = ph_data[name,"sex"])

ex_data_mlt %>% head
t.test(value~bt, data=ex_data_mlt)

test_results <- ex_data_mlt %>% 
  group_by(rowname) %>% 
  summarize(
    tval=t.test(value~bt)$statistic,
    pval=t.test(value~bt)$p.value,
    )


#### test all
data(ALL)
## data
ex_data <- exprs(ALL)
ph_data <- pData(ALL)[,c("cod", "sex", "BT")]

## remove missing | duplicated genes
ph_data <- ph_data[complete.cases(ph_data),]
feature_names <- rownames(ex_data)
gene_names <- unlist(as.list(hgu95av2SYMBOL[feature_names]))
idx <- which(is.na(gene_names) | duplicated(gene_names))
ex_data <- as.data.frame(ex_data[-idx,])
rownames(ex_data) <- gene_names[-idx]
dim(ex_data)

ex_data_mlt <- ex_data %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  mutate(bt = ph_data[name,"sex"]) %>% 
  drop_na()

test_results <- ex_data_mlt %>% 
  group_by(rowname) %>% 
  summarize(
    tval=t.test(value~bt)$statistic,
    pval=t.test(value~bt)$p.value,
  )

dim(test_results)
head(test_results)
sig_results <- test_results %>% 
  filter(pval<0.01) 

sel_genes <- ex_data_mlt %>% 
  filter(rowname %in% sig_results$rowname)

sel_genes %>% 
  group_by(rowname) %>% 
  ggplot(aes(x=bt, y=value, group=bt)) +
  facet_wrap(~rowname, ncol=9, scales="free") +
  geom_boxplot() +
  theme(
    axis.text.x = element_text(angle = 90, size=8, hjust = 1, vjust=0.5)
  )



### testing (anova)
lm(data=ex_data_mlt, formula = value~bt)





### ngs

#library(Rsamtools)
#library(ShortRead)
#library(GenomicAlignments)

refseq <- readDNAStringSet("ngs/dmpR_GESSv4.fasta")
pfm <- consensusMatrix(refseq)
pfm_atgc <- t(pfm[1:4,])
nt <- unlist(strsplit(as.character(refseq), split=""))

load("./ngs/R1_target_freq.Rdata")
ref_freq1 <- pfm_atgc*target_freq

load("./ngs/R2_target_freq.Rdata")
ref_freq2 <- pfm_atgc*target_freq

load("./ngs/R3_target_freq.Rdata")
ref_freq3 <- pfm_atgc*target_freq

load("./ngs/R4_target_freq.Rdata")
ref_freq4 <- pfm_atgc*target_freq

mydata <- data.frame(r1=rowSums(ref_freq1), 
                     r2=rowSums(ref_freq2), 
                     r3=rowSums(ref_freq3), 
                     r4=rowSums(ref_freq4))

mydata2 <- mydata %>% 
  apply(1, scale) %>% 
  t %>% as.data.frame %>% 
  mutate(pos=nrow(.):1, nt=nt) 

library(NbClust)
nc <- NbClust(mydata, min.nc=4, max.nc=7, method="kmeans")
barplot(table(nc$Best.n[1,]),
        xlab="Numer of Clusters", ylab="Number of Criteria",
        main="Number of Clusters Chosen by 26 Criteria")
table(nc$Best.n[1,])

cl <- kmeans(mydata, centers = 4, iter.max = 10000)
mydf3 <- mydata2 %>% 
  mutate(cl=factor(cl$cluster)) %>% 
  group_split(cl)

mydf5 <- mydf4_group %>% 
  lapply(function(x){ mutate(x, int_pos=c(nrow(x):1)) }) %>% do.call(rbind, .)
head(mydf5)
colnames(mydf5)[1:4] <- c("R2", "R3", "R4", "R5")
mydf5_melt <- melt(mydf5, id.vars=c("pos", "cl", "int_pos", "nt", "wtnt", "R4nt"))

tiff("Figure.tiff", width = 15, height = 5.5, units = 'in', res = 300, compression = 'lzw')
ggplot(mydf5_melt, aes(x=int_pos, y=variable, fill=value)) +
  scale_fill_viridis(name="FR",
                     option = 'C',
                     direction = 1,
                     na.value = "grey93") +
  geom_tile(color="black", size=0.5) +
  #geom_tile() +
  #facet_grid(cl~.) +
  theme_ipsum_rc() +
  annotate(geom = "text", x=1:length(aalab), y = -0.7, label =aalab, size = 4) +
  #annotate(geom = "text", x=1:length(aapos), y = -1.2, label =aapos, size = 4) +
  #annotate(geom = "text", x=1:length(aalab), y = -1.6, label =aalab, size = 4) +
  scale_x_continuous(breaks=1:length(ntlab), labels=ntlab) +
  coord_cartesian(ylim = c(0, 5), expand = FALSE, clip = "off") +
  theme(panel.spacing = unit(0.5, "lines"),
        plot.margin = unit(c(1,1,4,1), "lines"), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

dev.off()




subdat.mlt2 <- newdata.mlt2 %>% mutate(Inducer_conc_2 = factor(Inducer_conc_1))
levels(subdat.mlt2$Inducer_1) <- c("2AP", "2HB", "2NP", "3AP", "3HB", "3NP", "4AP", "4HB", "4NP", "benzene", "phenol")
ylimit <- c(-40000, 500000)
br.lab[1] <- "0"
tiff("FigureS52.tiff", width = 14, height = 10, units = 'in', res = 300, compression = 'lzw')
ggplot(subdat.mlt2, aes(x=Inducer_conc_1, y=value, color=main_date)) +
  geom_point(size=3, shape=16, color="gray") +
  geom_smooth(method="loess", level=0.95, span=1) + ## confidence interval
  scale_x_continuous(trans="log", breaks = br, labels=br.lab) +
  #facet_grid(Inducer_2~device_name, margins=FALSE, scales="fixed") +
  facet_grid(device_name~Inducer_1, margins=FALSE, scales="fixed") +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(angle = 90, size=10, hjust = 1, vjust=0.5),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        strip.text = element_text(size=12),
        panel.grid.major = element_line(colour = "#eeeeee"),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "gray", fill = "white"),
        #legend.position="none",
        panel.border = element_rect(color="grey", fill = NA)) +
  scale_y_continuous(limits = ylimit)+
  ggtitle("") + 
  labs(x=expression("Concentration ("*mu*"M)"), y="RFU", color="Date")
dev.off()
## 1400 x 1000




