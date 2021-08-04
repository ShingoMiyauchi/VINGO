#  Visual macro-synteny with VINGO

#                                       v1.0    Shingo Miyauchi 3Aug21
#-------------------
#   Descriptions
#-------------------
# < Info >
# The script generates a genome-wide plot with genome scaffolds, genes, TEs, synteny for the overview of macrosynteny among species compared. Visually Integrated Numerous Genres of Omics (VINGO) combines output from the visual integrative omics platfrom (a.k.a ShingoTools) such as SHIN+GO, TINGO, PRINGO, and SynGO. The figure made here was used in Looney et al. (2021).   

# < Input files >
# 1) Gene coordinates with genes, TEs in scaffold 1 to 10 from SynGO
# SynGO_Scaff1to10_GenomeFeature.csv

# 2) Scaffold size calculated from SynGO
# SynGO_Scaffolds_Size.csv

# 3) Synteny identified from SynGO
# SynGO_synteny_locations.csv

# 4) Scaffold of interest
# Glocon1_compared_with_7fungi.txt

# 5) CAZymes present in genomes
# CAZymes_CAZy_Total.csv

# < Overview of Methods >
# 1) Preparation 
# 2) Make a coordinate table 
# 3) Make a scaffold size table 
# 4) Make a table for synteny
# 5) Set up colours for variables
# 6) Make GRanges objects
# 7) Visualise macro synteny

# < R version + packages required >
# We used R version 3.6.3
# dplyr 1.0.2
# karyoploteR 1.12.4
# scales 1.1.0
#-------------------
# Global trun off for factorisation
options(stringsAsFactors=FALSE)

# Load R packages
library(dplyr)
library(karyoploteR)
library(scales)

#================
# 1) Preparation 
#================
# Load selected speciese with scaffolds  
interest.df <- read.delim(list.files(pattern="Glocon1_compared_with_7fungi.txt"))

# Insert scaffold name manually 
scaffold_name <- "Glocon1_scaffold1"
  
# Make a unique fungus + scaffold ID
uniqID <- paste0(interest.df$FungalID %>% gsub("\\d+|_", "", .), 
                 "_",  interest.df$Scaffold) %>% print()

# Load total CAZyme list
cazy.df <- read.csv("CAZymes_CAZy_Total.csv")

# Make a key column for merging with other data frames
cazy.df$key <- paste0(cazy.df$JGI_ID, ".",  cazy.df$ProteinID)
  
# Load genome coordinate + genes + TEs 
cord.df <- read.csv(list.files(pattern= "SynGO_Scaff1to10_GenomeFeature.csv"))

#============================
# 2) Make a coordinate table 
#============================
# Make a key column with fungus + proteinID 
cord.df$key <- paste0(cord.df$fungus, ".", cord.df$proteinID)

# Merge cazy.df with cord.df
m.df <- merge(cord.df, cazy.df[, c("key","Annot")], all.x=T)

# Make a proper column
m.df$CAZyTotal <- m.df$Annot

# Replace CAZy with "" in SSP column (because it is purely duplicated)
m.df$CAZyTotal[m.df$type == "SSP"] <- ""

# Replace all NA with ""
m.df[is.na(m.df)] <- ""

# Change the order of colums - important for transformation into Grange object
cord.df <- m.df[, c("scaffold", "start", "end", "strand",
                       "proteinID", "type", "annot","CAZyTotal", "fungus")]

# Make a column wiht unique fungalID + scaffold ID
cord.df$Scaffold <- paste0(cord.df$fungus %>% gsub("\\d+|_", "", .), cord.df$scaffold %>% gsub("scaffold", "", .))

# Extract only scaffolds of interest
scaf.df <- cord.df[cord.df$Scaffold %in% uniqID, ]

# Clean scaf.df
scaf.df <- scaf.df[, c("Scaffold", "start", "end", "type", "annot", "CAZyTotal")]

# Remove unwanted R objects
rm(cord.df, m.df, cazy.df)

#====================================
# 3) Make a scaffold size table 
#====================================
# Load scaffold size information
scaf.size.df <- read.csv(list.files(pattern="SynGO_Scaffolds_Size.csv"))

### Correct scaffold IDs - JGI genome labels are sometimes inconsistent
temp.ls <- lapply(scaf.size.df$FungalID %>% unique(), function(X) scaf.size.df[scaf.size.df$FungalID == X,])

# Make sure scaffolds are in large to small order
temp.ls <- lapply(temp.ls, function(X) arrange(X, desc(Length)))

# Convert scaffold IDs into "scaffold_#" format
temp.ls <- lapply(temp.ls, function(X) 
  data.frame(FungalID= X$FungalID, 
             Scaffold= paste0("scaffold_", rep(1:length(X$ScaffoldID))), 
             Start= rep(1),
             End= X$Length))

# Put them back to a single data frame 
scaf.size.df <- do.call(rbind, temp.ls)

# Make a unique ID column with fungiID + scaffoldID
scaf.size.df$Scaffold <- paste0(gsub("\\d+|_", "", scaf.size.df$FungalID),  gsub("scaffold", "", scaf.size.df$Scaffold))

# Make scaffold size table for conversion into GRanges 
scaf.size.df <- scaf.size.df[,c("Scaffold", "Start", "End", "FungalID")]

# Extract only for scaffolds of interest
scaf.size.df <- scaf.size.df[scaf.size.df$Scaffold %in% uniqID, ] 

# Arrange by evolutionary order based on uniqID from interest.df
scaf.size.df$Scaffold <- factor(scaf.size.df$Scaffold, levels= uniqID)
scaf.size.df <- scaf.size.df %>% arrange(Scaffold)

# Remove temp.ls
rm(temp.ls)

#===========================
# 4) Make a table for synteny 
#===========================
# Load macrosynteny information
syn.df <- read.csv(list.files(pattern="SynGO_synteny_locations.csv"))

#------------------------------
# Select scaffolds of interest with interest.df
#------------------------------
funID <- interest.df$FungalID
scafID <- interest.df$Scaffold
pairID <- unique(syn.df$Pair)

# Extract pairIDs in syn.df
ls <- list()
for(i in 1:(length(funID)-1)) ls[[i]] <- pairID[grep(funID[i], pairID)] %>% grep(funID[i+1], ., value = T) %>% print()

# Extract synteny by selected pairIDs
temp.ls <- lapply(ls, function(X) syn.df[syn.df$Pair %in% X, ])

# name temp.ls 
names(temp.ls) <- unlist(ls)

#--------------------------------
# Extract fungi + scaffolds of interest from temp.ls
#--------------------------------
start.syn.df.ls <- list()
end.syn.df.ls <- list()

for (k in 1:length(temp.ls)) {
  # Put separate two and put them vertically
  start.df <- temp.ls[[k]][,c("Fungus.1", "Scaffold.1", "Start.1", "End.1")]
  end.df <- temp.ls[[k]][,c("Fungus.2", "Scaffold.2", "Start.2", "End.2")] 
  
  # Make a unique ID column with fungus ID + scaffold ID
  start.df$ID <- paste0(start.df[, "Fungus.1"] %>% gsub("\\d+|_", "", .), "_", start.df[, "Scaffold.1"] %>% gsub("^\\s+|\\s+$", "", .))
  
  end.df$ID <- paste0(end.df[, "Fungus.2"] %>% gsub("\\d+|_", "", .), "_", end.df[, "Scaffold.2"] %>% gsub("^\\s+|\\s+$", "", .))
  
  # IMPORTANT - Get row IDs containing ID of interest + insect row IDs
  row.id <- intersect(which(start.df$ID %in% uniqID), 
                      which(end.df$ID %in% uniqID))
  
  # Separate into start and end link data frames for plotting
  start.syn.df.ls[[k]] <- start.df[row.id, c("ID", "Start.1", "End.1")] 
  end.syn.df.ls[[k]] <- end.df[row.id, c("ID", "Start.2", "End.2")]
}

# Name temp.ls 
names(start.syn.df.ls) <- names(temp.ls)
names(end.syn.df.ls) <- names(temp.ls)

# Remove intermediate R objects 
rm(temp.ls, ls, start.df, end.df, interest.df, syn.df, funID, i, k, row.id, uniqID, scafID, pairID)

#============================
# 5) Set up colours for secretome group
#============================
# Select colours to use in plots
iro <- c("#E41A1C", "#377EB8", "#4DAF4A", "#DDD23B", "#268785", "#D0104C", "#F75C2F", "darkorange") %>% sapply(function(X) alpha(X, 0.6)) %>% print()

# Check the colours
plot(rep(1,length(iro)),col=iro,pch=19,cex=3)

#=========================
# 6) Make GRanges objects
#=========================
# Extract scaffold size & Create a GRanges object for scaffold size
scaf.size.gr <- toGRanges(scaf.size.df)

# Separate TEs and genes + Convert into Grange object
# 1) TEs
scaf.te.gr <- scaf.df[scaf.df$type == "TE/Repeat sequences", ] %>% toGRanges()
 
# 2) Genes
scaf.gene.gr <- scaf.df[scaf.df$type != "TE/Repeat sequences", ] %>% toGRanges()

#  Synteny links
start.gr <- do.call(rbind, start.syn.df.ls) %>% toGRanges()
end.gr <- do.call(rbind, end.syn.df.ls) %>% toGRanges()

#==========================
# 7) Visualise macro synteny
#==========================
# Function for visual macro synteny
plot.fun <- function(GeneTerm, Scaf.length) {
  
  # Add scaffolds 
  kp <- plotKaryotype(genome= scaf.size.gr,
                      ideogram.plotter = NULL,
                      # Set up panels to use upper and lower panels for plots
                      plot.type = 2, 
                      col="grey30")
  
  # Visualise scaffold as a line
  kpAddCytobandsAsLine(kp, lwd=1)
  
  # Add length measure
  kpAddBaseNumbers(kp, tick.dist = Scaf.length, 
                   minor.tick.dist = Scaf.length/1e+01, 
                   add.units=T, col="grey60") 
  
  # Add synteny link
  kpPlotLinks(kp, data=start.gr, data2=end.gr, 
              col= alpha("lightblue1", 0.15), 
              border = alpha("lightblue1", 0.15), 
              data.panel=2, 
              r0= 0)
  
  # Gene density map
  kpPlotDensity(kp, data=scaf.gene.gr,
                window.size = Scaf.length / 1e+02,
                col= alpha("gray85",0.5),
                border= NA,
                r0= 0, r1= 0.6)
  
  # CAZyme total map 
  if (GeneTerm == "CAZymeTotal") {
    kpPlotDensity(kp, 
                  data=scaf.gene.gr[grep(".+", scaf.gene.gr$CAZyTotal),], 
                  window.size = Scaf.length / 1e+02,
                  col= iro[3],
                  border=NA,
                  r0= 0.1, 
                  r1= 0.7)
  }
  
  # Label genes coding for CAZymes 
  if (GeneTerm == "CAZymeTotal") {
    kpPlotMarkers(kp, chr=scaf.df[grep(".+", scaf.df$CAZyTotal), "Scaffold"], 
                  x= scaf.df[grep(".+", scaf.df$CAZyTotal), "start"], 
                  labels=scaf.df[grep(".+", scaf.df$CAZyTotal), "CAZyTotal"],
                  line.color = alpha("gray80",0.7), 
                  label.color = alpha("gray40",0.8),
                  text.orientation = "vertical",
                  ignore.chromosome.ends= F,
                  #label.dist = 0.000001,
                  # arrow bottom, middle, top
                  marker.parts=c(0, 0.2, 0.1),
                  cex= 0.45, 
                  r0= 0.3,
                  r1= 0.8)
  }
  
  # Total TE coveage in lower panel
  kpPlotDensity(kp, data=scaf.te.gr,
                window.size = Scaf.length / 1e+02,
                col= alpha("gray80",0.5),
                border= NA,
                r0= 0,
                r1= 1,
                data.panel = 2)
}

#------------
# Save plots
#------------
pdf (file=paste0(getwd(),"/" , "CAZymeTotal","_", scaffold_name, "_MacroSynteny",  ".pdf"),  height = 8, width = 7)
plot.fun("CAZymeTotal", 1e+06)  
dev.off()
