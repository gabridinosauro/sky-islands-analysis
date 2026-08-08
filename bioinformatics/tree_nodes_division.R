library(ape)
library(DECIPHER)
library(phangorn)

setwd("/Users/gabri/Desktop/mouse/SILVA_tree")


#1 Visualize Alignment -----------
### First thing I Want to check how does the alignment looks like
### read the allignment in, to visualize it 
aln <- Biostrings::readDNAStringSet("allignment_trimmed.fasta")
#aln <- read.dna("allignment_trimmed.fasta", format = "fasta")
# Check dimensions
length(aln)  # Should be ~47000
hdrs <- names(aln) ## extract headers
head(hdrs, 20)   # quick peek at the first 20 headers
noASV_hdrs <- hdrs[!grepl("ASV", hdrs, ignore.case = TRUE)]
ASVs <- paste0("ASV", 1:50) ### jest get the first more or less 50 most abudant ASVs
set.seed(123)  # optional, for reproducibility
rand_noASV <- sample(noASV_hdrs, 10) ### take 10 random Archaeae
# logical index: which sequence names are ASVs
idx <- names(aln) %in% c(ASVs, rand_noASV)
names(aln)[idx]
# visualize only ASVs
BrowseSeqs(aln[idx], highlight = 1)




#2 Tree - 3 species ligilactobacillus ----
#### Now read the tree 
tt  = read.tree("asv_tree_no_fastest.nwk")
length(tt$tip.label)
"ASV2" %in% tt$tip.label # not there
"ASV1" %in% hdrs # not there!!!!
### something is wrong here....

taxtab = as.data.frame(readRDS("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/all_data_run/run_dada2/Taxtab.RDS"))
taxtab_lc <- taxtab[!is.na(taxtab$Species) &
                      taxtab$Species == "Lactobacillus_sp007570935(RS_GCF_007570935_1", ]
unique(taxtab$Phylum)
# vector of tips you want to keep
tips_to_keep <- rownames(taxtab_lc)[rownames(taxtab_lc) %in% tt$tip.label]
# prune the tree to ONLY those tips
tt_min <- keep.tip(tt, tips_to_keep)
plot(tt_min)
# species 1
taxtab_lc <- taxtab[!is.na(taxtab$Species) &
                      taxtab$Species == "Lactobacillus_sp007570935(RS_GCF_007570935_1", ]
tips_lc <- intersect(rownames(taxtab_lc), tt$tip.label)
# species 2
taxtab_la <- taxtab[!is.na(taxtab$Species) &
                      taxtab$Species == "Ligilactobacillus_animalis(RS_GCF_001434535_1", ]
tips_la <- intersect(rownames(taxtab_la), tt$tip.label)
# species 3
taxtab_lap <- taxtab[!is.na(taxtab$Species) &
                       taxtab$Species == "Ligilactobacillus_apodemi(RS_GCF_001434405_1", ]
tips_lap <- intersect(rownames(taxtab_lap), tt$tip.label)
# combine (use unique to avoid duplicates)
tips_all <- unique(c(tips_lc, tips_la, tips_lap))
# prune tree
tt_min <- keep.tip(tt, tips_all)
# --- create tip color vector BEFORE plotting
tip_cols <- rep("black", length(tt_min$tip.label))
names(tip_cols) <- tt_min$tip.label
# assign colors only if those tip names exist in the pruned tree
if(length(tips_lc) > 0) tip_cols[intersect(tips_lc, names(tip_cols))] <- "red"
if(length(tips_la) > 0) tip_cols[intersect(tips_la, names(tip_cols))] <- "blue"
if(length(tips_lap) > 0) tip_cols[intersect(tips_lap, names(tip_cols))] <- "darkgreen"
# plot
plot(tt_min, tip.color = tip_cols)







# 3 collpase tree by sipport ---------------
### I have to collapse the nodes
collapseBySupport <- function(tree, cutoff = 0.5) {
  support <- as.numeric(tree$node.label)
  bad_nodes <- which(is.na(support) | support < cutoff)
  edge_idx <- which(tree$edge[,2] %in% (Ntip(tree) + bad_nodes))
  tree$edge.length[edge_idx] <- 0
  tree <- di2multi(tree, tol = 1e-8)
  collapsed_tree <-  collapse.singles(tree)
  return(collapsed_tree)
}


collapsed_tree <- collapseBySupport(tt_min, cutoff = 0.5)
plot(collapsed_tree,tip.color = tip_cols)
### looks the same ok
collapsed_tree <- collapseBySupport(tt, cutoff = 0.5)





#  4 root the tree using archaea -------------
is.rooted(collapsed_tree)
labs = tt$tip.label
noASV_hdrs <- labs[!grepl("ASV", labs, ignore.case = TRUE)]
tt_rooted <- root(
  collapsed_tree,
  outgroup = noASV_hdrs,
  resolve.root = TRUE
)
ASVs <- paste0("ASV", 1:200)
set.seed(123)  # optional, for reproducibility
rand_noASV <- sample(noASV_hdrs, 10)
rand_noASV

# logical index: which sequence names are ASVs
idx <- names(aln) %in% c(ASVs, rand_noASV)
#keep <- sample(tt_rooted$tip.label, 500) ### random looks weird, but I think the problem is that some low abudancne ASVs are actually chimeras and ruin thetip.labeltip.labeltip.label
keep <- names(aln)[idx]
subtree <- drop.tip(tt_rooted,
                    setdiff(tt_rooted$tip.label, keep))

plot(subtree, cex=0.3)
head(taxtab)
subtree$tip.label

fam <- factor(taxtab[subtree$tip.label, "Phylum"])   # <- FORCE factor
library(Polychrome)

cols <- setNames(
  palette36.colors(length(levels(fam))),
  levels(fam)
)

tipcols <- cols[fam]                         # color for each tip
tipcols
plot(subtree, tip.color = tipcols, cex = 0.3)
legend("bottomleft",
       legend = levels(fam),
       col = cols,
       pch = 3,
       cex = 0.6,                  # shrink text
       pt.cex = 0.7,               # shrink symbols
       x.intersp = 0.4,            # tighten horizontal spacing
       y.intersp = 0.6,            # tighten vertical spacing
       xjust = 1,                  # RIGHT-align text inside box
       yjust = 1,
       bty = "o",                  # draw box
       box.lwd = 0.7,              # thinner box line
       inset = 0.02)

### remove archaea -----------
tt_no_archaea <- drop.tip(tt_rooted, noASV_hdrs)
saveRDS(tt_no_archaea,"/Users/gabri/Desktop/mouse/SILVA_tree/tt_no_archaea.RDS")
idx <- names(aln) %in% c(ASVs)
keep <- names(aln)[idx]
subtree <- drop.tip(tt_rooted,
                    setdiff(tt_rooted$tip.label, keep))

fam <- factor(taxtab[subtree$tip.label, "Phylum"])   # <- FORCE factor
library(Polychrome)

cols <- setNames(
  palette36.colors(length(levels(fam))),
  levels(fam)
)

tipcols <- cols[fam]                         # color for each tip
tipcols
plot(subtree, tip.color = tipcols, cex = 0.3)
legend("bottomright",
       legend = levels(fam),
       col = cols,
       pch = 3,
       cex = 0.6,                  # shrink text
       pt.cex = 0.7,               # shrink symbols
       x.intersp = 0.4,            # tighten horizontal spacing
       y.intersp = 0.6,            # tighten vertical spacing
       xjust = 1,                  # RIGHT-align text inside box
       yjust = 1,
       bty = "o",                  # draw box
       box.lwd = 0.7,              # thinner box line
       inset = 0.02)




### Ok great now
# 5 calculate division to 1-4th of length -------------------------------------
### the logic here is to divide the tree and take only those nodes 1/4 away from the top to the leaf...
### The problem is that different clades seem at a differet distance from the root... so I need to normalize the distance for each node...
### And take those 1/4 away from the tip they countain more far away from the root
## 5.1 calculate distance from root -------------
d <- node.depth.edgelength(tt_no_archaea) ### calculate distance from root

## 5.2 calculate for each node which leaf it contains -------------------------
Ntip <- length(tt_no_archaea$tip.label) # number of tips
edges <- tt_no_archaea$edge# edge matrix: parent -> child
children <- split(edges[,2], edges[,1])# children list for each node
# function to get descendant tips (fast, recursive)
get_desc_tips <- function(node) {
  kids <- children[[as.character(node)]]
  if (is.null(kids)) return(character(0))
  
  tips <- integer(0)
  for (k in kids) {
    if (k <= Ntip) {
      tips <- c(tips, k)
    } else {
      tips <- c(tips, get_desc_tips(k))
    }
  }
  tips
}
# internal nodes
internal_nodes <- (Ntip + 1):(Ntip + tt_no_archaea$Nnode)
# compute descendant tips (FAST)
pb <- txtProgressBar(min = 0, max = length(internal_nodes), style = 3)
descendant_tips <- lapply(seq_along(internal_nodes), function(i) {
  setTxtProgressBar(pb, i)
  get_desc_tips(internal_nodes[i])
})
close(pb)

#### Let's do a check
### select three tips and calcualte their ancestry:
# pick 3 tips at random
set.seed(1)
tips_to_check <- sample(1:Ntip, 3)
#tips_to_check
get_ancestors <- function(tip, edges) {
  ancestors <- integer(0)
  current <- tip
  
  repeat {
    parent <- edges[edges[,2] == current, 1]
    if (length(parent) == 0) break
    ancestors <- c(ancestors, parent)
    current <- parent
  }
  
#  ancestors
}
ancestors_list <- lapply(tips_to_check, get_ancestors, edges = edges)
names(ancestors_list) <- paste0("tip_", tips_to_check)
#ancestors_list
node <- 47837 ### try a combiantion of thee
tip  <- 43307
# find index of this node in internal_nodes
idx <- which(internal_nodes == node)
tip %in% descendant_tips[[idx]]
### seem to work... great
#saveRDS(descendant_tips,"descendant_tips.RDS")
descendant_tips = readRDS("descendant_tips.RDS")

# convert indices to labels
descendant_tips <- lapply(descendant_tips, function(x)
  tt_no_archaea$tip.label[x]
)
names(descendant_tips) <- internal_nodes



## 5.3 now calculate max distances distances -----------------------
### this code calculate for each node how far away their most distant leaf is from the root 
tip_index <- setNames(seq_len(Ntip), tt_no_archaea$tip.label)
pb <- txtProgressBar(min = 0, max = length(descendant_tips), style = 3)
descendant_idx <- lapply(seq_along(descendant_tips), function(i) {
  setTxtProgressBar(pb, i)
  tip_index[ descendant_tips[[i]] ]
})
close(pb)
node_score <- vapply(descendant_idx, function(idx)
  if (length(idx) == 0) NA_real_ else max(d[idx]),
  numeric(1)
)


## 5.4 Now calculate the distance of each node from the leaf ----------
node_root_distance <- d[internal_nodes]

## 5.5 calcualate position of node relative to tip furthest away -----
node_score[1:10]
ape::Ntip(tt_no_archaea) ## #ok correct

nodes_df <- data.frame(
  node_id = internal_nodes,
  node_score = node_score,
  node_distance = d[internal_nodes]
)

#### Ok now calcualte percet distance

nodes_df$fraction  = nodes_df$node_distance / nodes_df$node_score
nodes_df$one_forth  = ifelse(nodes_df$fraction > 0.75, TRUE, FALSE)
nodes_df$one_eight  = ifelse(nodes_df$fraction > 0.875, TRUE, FALSE)
nodes_df$one_sixteeen  = ifelse(nodes_df$fraction > 0.9375, TRUE, FALSE)

head(nodes_df)
table(nodes_df$one_forth)
table(nodes_df$one_eight)
table(nodes_df$one_sixteeen)

### Ok now plot
# Assign node IDs as node labels to the full tree
tt_no_archaea_alt_lab = tt_no_archaea
n_tips <- length(tt_no_archaea_alt_lab$tip.label)
length(tt_no_archaea_alt_lab$node.label) == nrow(nodes_df)
tt_no_archaea_alt_lab$node.label <- nodes_df$node_id
subtree_alt_lab <- drop.tip(tt_no_archaea_alt_lab,
                    setdiff(tt_no_archaea_alt_lab$tip.label, keep))
  
  # Then color the nodes
tree_node_ids <- as.numeric(subtree_alt_lab$node.label)
matched_idx <- match(tree_node_ids, nodes_df$node_id)
node_colors <- ifelse(nodes_df$one_forth[matched_idx], "red", "black")

# Plot
plot(subtree_alt_lab, tip.color = tipcols, cex = 0.3)
nodelabels(pch = 16, col = node_colors, cex = 0.5)


### cool, it works.
### save the list of nodes
#saveRDS(nodes_df, "nodes_df.RDS")
# 6 how many tips and which tips are above each node ----
### this is to do the filtering later in the day...
# If TRUE, add it as a column
nodes_df$descendants <- descendant_idx
head(nodes_df)
nodes_df$descendants[53]
nodes_df$n_descendants <- lengths(nodes_df$descendants)

#saveRDS(nodes_df, "nodes_df.RDS")

### ok supercool
# 7 collapse the nodes ------------------------------------------------------
### this last part collapsing the tree into all 1/4 nodes
### so the logic sort of creating phylogenetically aware OTUs so to say
### to do it, I need to collapse all nodes above the last excluded node.
head(nodes_df)
parent_lookup <- setNames(edges[,1], edges[,2])# Create parent lookup from edges
nodes_df$parent <- parent_lookup[as.character(nodes_df$node_id)]# Assign father node
### HERE change to collapse different values -------------------------------
# Create a lookup for one_forth values by node_id
one_forth_lookup <- setNames(nodes_df$one_forth, nodes_df$node_id)
one_eight_lookup <- setNames(nodes_df$one_eight, nodes_df$node_id)
one_soxt_lookup <- setNames(nodes_df$one_sixteeen, nodes_df$node_id)

# Find transition nodes: where one_forth is TRUE but parent's one_forth is FALSE
nodes_df$is_transition <- nodes_df$one_forth & 
  !one_forth_lookup[as.character(nodes_df$parent)]
nodes_df$is_transition[is.na(nodes_df$parent)] <- FALSE

nodes_df$is_transition_1_8 <- nodes_df$one_eight & 
  !one_eight_lookup[as.character(nodes_df$parent)]
nodes_df$is_transition_1_8[is.na(nodes_df$parent)] <- FALSE

nodes_df$is_transition_1_16 <- nodes_df$one_sixteeen & 
  !one_soxt_lookup[as.character(nodes_df$parent)]
nodes_df$is_transition_1_16[is.na(nodes_df$parent)] <- FALSE

table(nodes_df$is_transition_1_16 )
table(nodes_df$is_transition_1_8 )

# Get the transition nodes
transition_nodes <- nodes_df[nodes_df$is_transition, ]

# View results
head(transition_nodes)


subtree_alt_lab <- drop.tip(tt_no_archaea_alt_lab,
                    setdiff(tt_no_archaea_alt_lab$tip.label, keep))
  
  # Then color the nodes
tree_node_ids <- as.numeric(subtree_alt_lab$node.label)
matched_idx <- match(tree_node_ids, nodes_df$node_id)
node_colors <- ifelse(nodes_df$is_transition[matched_idx], "red", "black")

# Plot
plot(subtree_alt_lab, tip.color = tipcols, cex = 0.3)
nodelabels(pch = 16, col = node_colors, cex = 0.5)

## cool some nodes missing so only few visualizd, but we have ur winners....
### save the final dataframe a
saveRDS(nodes_df, "nodes_df.RDS")
nodes_df = readRDS("/Users/gabri/Desktop/mouse/SILVA_tree/nodes_df.RDS")
