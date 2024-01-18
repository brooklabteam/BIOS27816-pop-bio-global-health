# Intro to phylogenetic modeling - Global epidemiology 2024
# Influenza phylogeny
# Script written by Gwen Kettenburg

# Okay! you have your phylogeny file, what can you do with it? With R, you can flip it, change the 
# outgroup, highlight clades, change the size of text and bootstrap numbers, change the top colors and symbols, 
# add vector images to the tips, and so much more! FigTree is a GUI interface and can be used to quickly look at 
# the tree before doing more with it, and you can edit in it as well, but R is much more powerful and you have 
# many more options. 

# First, let's load some packages in R that are specific to tree editing and visualization

rm(list=ls()) #clear your environment

# ggtree runs from biocunductor, so you may need to install it before getting ggtree
## install.packages("BiocManager")
#library(BiocManager)

library(tidyverse) #package for cleaning things up
library(ggplot2) #package for making nice looking plots
library(ggtree) #package for manipulating phylogenetic trees
library(ape) #another package for manipulating phylogenetic trees
library(ggnewscale) #package for adding scales on the plots
library(cowplot) #package for stitching multiple plots together into one


# Second, load the tree files into R by setting working directory to where the file is located. 

# You can go up top to the options, click on Session -> Set Working Directory -> Choose Directory

# You can also set a path to your file, example below is what I do when working with my R scripts, it makes it
# easier for anyone working with your file down the line. 

#homewd= "/Users/gwenddolenkettenburg/Desktop/"
#setwd("~/Desktop/intro_phylo_paris/RAxML_output")


# We're going to make two trees to compare, one that is for hemagglutinin protein of influenza which is involved with cell entry, and
# PB2 protein which encodes polymerase, a conserved protein involved with viral replication

pb2.tree <- read.tree("pb2_h1n1_seq_newick") #loading the tree files into R
ha.tree <- read.tree("ha_h1n1_seq_newick") #loading the tree files into R

root.pb2.tree <- root(pb2.tree, which(pb2.tree$tip.label == "NC_006505_Infectious_salmon_anemia_virus")) #setting the root as Isavirus PB2
root.ha.tree <- root(ha.tree, which(ha.tree$tip.label == "AY601904_Infectious_salmon_anemia_virus")) #setting the root as Isavirus HA

plot(root.pb2.tree) #check how it looks quickly before doing anything else
plot(root.ha.tree) #check how it looks quickly before doing anything else

#They look different! Any ideas why? 

#load tree data prepared from elsewhere, so we can use it with the tree if needed
dat.ha <- read.csv(("ha_h1n1_france_human_swine_metadata.csv"), header = T, stringsAsFactors = F)
dat.pb2 <- read.csv(("pb2_h1n1_france_human_swine_metadata.csv"), header = T, stringsAsFactors = F)

#check that your metadata matches your tree data
nrow(dat.ha) #8
length(root.ha.tree$tip.label) #8
nrow(dat.pb2) #8
length(root.pb2.tree$tip.label) #8

# If you have a very large dataset of sequences from NCBI, it may be good to have a separate csv file of specific parameters
# to associate with the tree, like if you want to show sampling location, host, novel sequence, year sampled, and more. This
# is especially important in making viral phylogenies. I'll show you what those look like, but in the folder look for PB2 and HA metadata 
# files and take a look. We will not go over this today but reach out if you would like to learn how to add that.

#Okay now we have everything in R

#Let's explore what we can do with ggtree

#When adding features to plots, it is important to save versions of it that you can add to, p is traditionally used

p.pb2 <- ggtree(root.pb2.tree) + ggtitle("PB2") #base and most simple version of our phylogeny, adding a title
p.pb2 #check it out

p.ha <- ggtree(root.ha.tree) + ggtitle("HA") #base and most simple version of our phylogeny, adding a title
p.ha #check it out

#Can't infer much right now


#Let's add the tip labels back, and add some points to the tip labels too so we can make them colorful later on

p1.pb2 <- p.pb2 + geom_nodepoint() + geom_tiplab(size=3) + xlim(0,3) # adding blank tip points and original labels to change later
p1.pb2 # check it out

p1.ha <- p.ha + geom_nodepoint() + geom_tiplab(size=3) + xlim (0,3) # adding blank tip points and original labels to change later
p1.ha # check it out

#It's easy to add some personalization here, putting colors and shapes inside the parenthesis for the commands to make nodes and tip labels

# You pick the colors and shapes! Let's also make the outgroup text black so we can see it more easily, but first we need to manually define
# a color scheme, let's make the HA tree orange and the PB2 tree blue, and do different shades depending on whether the host is swine or human

colz.ha<-c("AHC68973_A_Lyon_1_12_2011_H1N1"="orange",
            "ACS91401_A_Paris_2709_2009_H1N1"="orange",
            "AGC13493_A_swine_Cotes_dArmor_110466_2010_H1N1"="orange4",
            "AHC68951_A_StEtienne_1139_2010_H1N1"="orange",
            "AGK62667_A_swine_France_CotesdArmor_0388_2009_H1N1"="orange4",
            "AKJ80511_A_swine_France_Finistere_0011_2011_H1N1"="orange4",
            "AID48624_A_swine_Marseille_2260_1980_H1N1"="orange4",
            "AY601904_Infectious_salmon_anemia_virus"="black")

colz.pb2<-c("AKJ80506_A_swine_France_Finistere_0011_2011_H1N1"="deepskyblue4",
            "AKJ83080_A_swine_France_CotesdArmor_0388_2009_H1N1"="deepskyblue4",
            "AID48619_A_swine_Marseille_2260_1980_H1N1"="deepskyblue4",
            "AHC68947_A_StEtienne_1139_2010_H1N1"="deepskyblue",
            "ACS91418_A_Paris_2709_2009_H1N1"="deepskyblue",
            "AHC68969_A_Lyon_1_12_2011_H1N1"="deepskyblue",
            "AGF37528_A_swine_Cotes_dArmor_110466_2010_H1N1"="deepskyblue4",
            "NC_006505_Infectious_salmon_anemia_virus"="black")

#Now, we can make the phylogeny nice looking

#Here we highlight the node points, we highlight the branch tip points, and color the actual word labels
p2.ha <- p1.ha + geom_nodepoint(color="lightgreen", shape=16, size=4) + geom_tippoint(color="indianred1", shape=8, size=1) +
  geom_tiplab(color=colz.ha, size=3)
p2.ha # check it out

p2.pb2 <- p1.pb2 + geom_nodepoint(color="lightgreen", shape=16, size=4) + geom_tippoint(color="indianred1", shape=8, size=1) +
  geom_tiplab(color=colz.pb2, size=3)
p2.pb2 # check it out

#Now the above is how phylogenies are normally shown, in this shape...but you can make round phylogenies too


# Let's change the shape to a fan, first we need a version of the tree without the labels, since they will be the wrong direction
p3.ha <- p.ha + geom_nodepoint(color="lightgreen", shape=16, size=5) + geom_tippoint(color="indianred1", shape=8, size=1)
p3.ha #check it out

p3.pb2 <- p.pb2 + geom_nodepoint(color="lightgreen", shape=16, size=5) + geom_tippoint(color="indianred1", shape=8, size=1)
p3.pb2 #check it out

#ggtree feature opens the tree specified tree view into a fan shape
#and rotates it 180 degrees, we are re-adding the labels on top of the fan shape
p4.ha <- open_tree(p3.ha, 180) + geom_tiplab()
p4.ha # check it out

p4.pb2 <- open_tree(p3.pb2, 180) + geom_tiplab()
p4.pb2 # check it out

#Sure it looks nice but it's harder to read so let's go back to playing with the first tree with the colored labels

# We can add some clade labels! Let's look at the tree again and think, what could be some clades?

# First, we can add a scale bar for the amino acid substitution rate, now we can actually analyze some things
p5.ha <- p2.ha+geom_treescale(0.1) 
p5.ha # check it out

p5.pb2 <- p2.pb2+geom_treescale(0.1)
p5.pb2 # check it out

#Hmmm it's hard to tell where the clades are in the PB2 tree since it's such a conserved region, but there seems to be
#at least one swine clade

# ggtree labels the nodes with numbers, so let's get those so we know where we are labeling the clades

# look at it this way, which is a table form and fairly easy to read
x <- as_tibble(root.ha.tree)
print(x, n=15) # see all the rows, 15 rows gives the whole data set
y <- as_tibble(root.pb2.tree)
print(y, n=15)

# Here we want to use the parent category for telling ggtree what to do, clades will all share a number


#Two ways you can highlight the clades we have specified

#Method number 1 - use a vertical bar and label
p6.ha <- p5.ha + geom_cladelab(node=15, label="H1N1 swine",  
                    offset = 0.8, textcolor='orange4')
p6.ha #check it out

p6.pb2 <- p5.pb2 + geom_cladelab(node=11, label="H1N1 swine",  
                               offset = 0.8, textcolor='deepskyblue4')
p6.pb2 #check it out


#Method number 2 - highlight the whole clade
p7.ha <- p5.ha + geom_hilight(node=15, fill="yellow", alpha=0.3, extend=1)
p7.ha #check it out

p7.pb2 <- p5.pb2 + geom_hilight(node=11, fill="yellow", alpha=0.3, extend=1)
p7.pb2 #check it out


#Let's put both trees side by side
h1n1<-plot_grid(p7.ha,p7.pb2,rel_widths = c(1,1), nrow = 1)
h1n1

#What kind of differences do you see?

#For further reading, which is a wonderful and easy resource: https://yulab-smu.top/treedata-book/index.html

