## Day1

x <- 10
y <- 20

# Alt + -
# <- 

2 + 2
((2-1)^2 + (1-3)^2 )^(1/2)
2 + 2; 2 - 2

sqrt((4+3)*(2+1))


x <- c(1, 3, 5, 7, 9)
y <- c(2, 4, 6, 8, 10)
x
y
z <- x + y

?mean
example("mean")
example()
help(package="UsingR")

data(package="UsingR")
mydat <- SAT
?SAT
data(SAT)

#The number of whale beachings in Texas during the 1990s
whale <- c(74, 122, 235, 111, 292, 111, 211, 133, 156, 79)
#Object `whale` is a data vector == (univariate) data set

# The size 
length(whale)
sum(whale)
sum(whale)/length(whale)
mean(whale)


whale^2
whale - mean(whale)
whale^2 - mean(whale)
sqrt(whale)

x <- 1
x <- c(1, 2)
x <- c(x, 3)

hip_cost <- c(10500, 45000, 74100, NA, 83500)
sum(hip_cost)
sum(hip_cost, na.rm=TRUE)
?sum
?mean
mean(hip_cost, na.rm=TRUE)

?precip
length(precip)
class(precip)
data(precip)

names(precip)
o <- order(names(precip))
names(precip)[o]
sort(names(precip))

test_scores <- c(100, 90, 80)
names(test_scores) <- c("Alice", "Bob", "Shirley")
test_scores[1]
test_scores[2]

precip[1]
precip[2:10]
precip[c(1,3,5)]
precip[-1]
precip["Seattle Tacoma"]
precip[2] <- 10


1:10
seq(1, 10, 1)
seq(1, by=1, to=10)
?seq
seq(1, by=0.1, to=10)
seq(1, 100, length.out=10)


rep(1, 4)
?rep
rep(1:4, 5)
rep(1:4, each=5)

### 3.1.2.1 exercise

x <- rep("a", 5)
x <- seq(1, 99, by=2)
x <- rep(1:3, each=3)
x <- rep(seq(1,3), each=3)
x <- c(1:5, 4:1)

ch <- c("Lincoln", "said", "and")

paste("X", 1:10)
paste("X", 1:10, sep="")
paste(c("A", "B"), 1:10, sep="")

paste("The", "quick", "brown", "fox")
paste(c("The", "quick", "brown", "fox"))
paste(c("The", "quick", "brown", "fox"), collapse=" ")


x <- 1:10
paste(x)
paste(x, collapse=":")


x <- c("Red", "Blue", "Yellow", "Green", "Blue", "Green")
y <- factor(x)
y
y[1] <- "Gold"
y

levels(y) <- c(levels(y), "Gold")
y[1] <- "Gold"
y

is.na(1)
is.numeric(1)
is.logical(TRUE)

x <- 1
is.numeric(x)
if(is.numeric(x)){
  # do something..1
}else{
  # do something...2
}

pi < 3

precip < 30
precip[which(precip < 30)]


sum(c(TRUE, TRUE))

x <- 1:100
x < 10 | x >90


head(rivers)
tail(rivers)

?rivers

head(rivers)


my_mean <- function(x){
  z <- sum(x, na.rm=T)/length(x)
  return(z)
}

x <- c(1:10, NA)

my_mean(x)


### Day2 2020 10 14

ex <- exprs(ALL)
ph <- pData(ALL)
dim(ph)
ph[1:10,1:5]
ph2 <- phenoData(ALL)
varMetadata(ph2)
pData(ph2)

library(GEOquery)
gds <- getGEO(filename=system.file("extdata/GDS507.soft.gz",package="GEOquery"))
class(gds)
methods(class="GDS")
Table(gds)
Columns(gds)

gsm <- getGEO(filename=system.file("extdata/GSM11805.txt.gz",package="GEOquery"))
Meta(gsm)


gse2553 <- getGEO('GSE2553',GSEMatrix=TRUE)
gse2553
class(gse2553)
length(gse2553)
class(gse2553[[1]])

mydat <- gse2553[[1]]
ex <- exprs(mydat)
dim(ex)  

## phenotype data
ph <- pData(mydat)
dim(ph)  
ph[1:10,1:5]
colnames(ph)
varMetadata(phenoData(mydat))

## featureData
fe <- fData(mydat)
dim(fe)
fe[1:10,1:5]
colnames(fe)
featureData(mydat)
varMetadata(featureData(mydat))
  

eset <- GDS2eSet(gds, do.log2=TRUE)
eset

## ALL
ex <- exprs(ALL)
dim(ex)
ex[1:10,1:5]

ph <- pData(ALL)
dim(ph)
colnames(ph)
varMetadata(phenoData(ALL))

featureData(ALL)
ids <- featureNames(ALL)[1:10]

library(hgu95av2.db)
?hgu95av2ENTREZID

xx <- as.list(hgu95av2ENTREZID)
xx[ids]
ALL

xx <- as.list(hgu95av2SYMBOL)
xx[ids]

subsx <- ex[1:10,]
rownames(subsx) <- xx[ids]
plot(subsx["MAPK3",])
hist(subsx["MAPK3",], br=30)


####
library(Biostrings)


x0 <- sample(DNA_BASES, 10, replace = T)
x1 <- paste(x0, collapse = "")
class(x1)
x1 <- DNAString(x1)
class(x1)
length(x1)
nchar(x1)
x2 <- toString(x1)
complement(x1)
reverseComplement(x1)
translate(x1)

x1[1:3]
subseq(x1, 3, 5)

gccont <- letterFrequency(x1, c("G", "C"), as.prob=TRUE)
sum(gccont)

## exercise
x0 <- sample(DNA_BASES, 30, replace = T)
x1 <- paste(x0, collapse = "")
x2 <- paste("ATG", x1, "TAG", sep="")
x0
table(x0)
x3 <- DNAString(x2)
class(x3)
sum(letterFrequency(x3, c("G", "C"), as.prob=T))


x0 <- c("CTC-NACCAGTAT", "TTGA", "TACCTAGAG")
x1 <- DNAStringSet(x0)
## 추가 삭제
length(x1)
nchar(x1)
x4 <- c(x1, DNAStringSet(x3))
x4[-1]
x4
subseq(x4, 2, 4)
x4[1]
x4[[1]]

## loot
for(i in 1:10){
  cat(i, "\n")
}

x0 <- rep("", 10)
for(i in 1:length(x0)){
  tmp <- paste(sample(DNA_BASES, 30, replace = T), collapse="")
  x0[i] <- paste("ATG", tmp, "TAG", sep="")
}
x0
x1 <- DNAStringSet(x0)


random_dna <- function(len){
  #require(Biostrings)
  tmp <- paste(sample(DNA_BASES, len, replace = T), collapse="")
  x0 <- paste("ATG", tmp, "TAG", sep="")
  #x0 <- DNAString(x0)
  return(x0)
}

random_dna(len=30)
random_dna(len=40)


x0 <- rep("", 10)
for(i in 1:length(x0)){
  x0[i] <- random_dna(len=30)
}
x0


x0 <- replicate(10, random_dna(30))
x1 <- DNAStringSet(x0)

x1 <- random_dna(200)
x2 <- DNAString(x1)


Views(x2, 1, width=20)
successiveViews(x2, width=rep(20, 2))
successiveViews(x2, width=rep(20, 10))


x1 <- DNAStringSet(x0)
names(x1) <- paste("Seq", 1:length(x1), sep="")
x1
writeXStringSet(x1, "myfastaseq.fasta", format="fasta")
myseq <- readDNAStringSet("myfastaseq.fasta", format="fasta")

###

library(rentrez)
help(package="rentrez")
browseVignettes("rentrez")

entrez_dbs()
entrez_db_summary("nuccore")

acc <- c("NC_001477", "NC_001474", "NC_001475", "NC_002640")
all_recs <- entrez_fetch(db="nuccore", id=acc[1], rettype="fasta")
?entrez_fetch


r_search <- entrez_search(db="pubmed", term="R Language")
all_recs <- entrez_fetch(db="pubmed", id=r_search$ids, rettype="fasta")
write(all_recs, file="mypub.txt")













