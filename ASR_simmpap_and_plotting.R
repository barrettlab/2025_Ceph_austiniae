

# tip data
# species	leafy	myco
# Aphyllorchis_montana_KU551262	leafless	FM
# Cephalanthera_austiniae_2010_Ce_stiniae	leafless	FM
# Cephalanthera_damasonium_MH590345	leafy	PM
# Cephalanthera_falcata_OR161968	leafy	PM
# Cephalanthera_humilis_KU551265	leafless	FM
# Cephalanthera_longibracteata_MH590346	leafy	PM
# Cephalanthera_longifolia_KU551263	leafy	PM
# Cephalanthera_nanchuanica_OR947462	leafy	PM
# Cephalanthera_rubra_MH590347	leafy	PM
# Diplandrorchis_sinica_MZ014629	leafless	FM
# Diplandrorchis_sinica_OP310037	leafless	FM
# Epipactis_albensis_MH590348	leafy	PM
# Epipactis_atrorubens_MH590349	leafy	PM
# Epipactis_gigantea_MH590350	leafy	AU
# Epipactis_helleborine_MH590351	leafy	PM
# Epipactis_helleborine_MK593537	leafy	PM
# Epipactis_helleborine_MK608776	leafy	PM
# Epipactis_mairei_KU551264	leafy	AU
# Epipactis_mairei_MG925367	leafy	AU
# Epipactis_microphylla_MH590352	leafy	PM
# Epipactis_palustris_MH590353	leafy	AU
# Epipactis_purpurata_MH590354	leafy	PM
# Epipactis_thunbergii_MN200387	leafy	AU
# Epipactis_veratrifolia_KU551267	leafy	AU
# Limodorum_abortivum_MH590355	leafless	PM
# Neottia_acuminata_KU551268	leafless	FM
# Neottia_camtschatea_KU551266	leafless	FM
# Neottia_cordata_MH590356	leafy	PM
# Neottia_fugongensis_KU551270	leafy	AU
# Neottia_japonica_MH321183	leafy	AU
# Neottia_japonica_MH321184	leafy	AU
# Neottia_listeroides_KU551272	leafless	FM
# Neottia_nidus_avis_JF325876	leafless	FM
# Neottia_ovata_KU551271	leafy	PM
# Neottia_pinetorum_KU551269	leafy	AU
# Neottia_suzukii_MH321185	leafy	AU
# Nervilia_fordii_ON515491	leafy	AU
# Palmorchis_pabstii_MH590357	leafy	AU
# Sobralia_callosa_KM032623	leafy	AU
# Spiranthes_sinensis_MK936427	leafy	AU


# tip data, binary format
# species	leafy	myco
# Aphyllorchis_montana_KU551262	1	2
# Cephalanthera_austiniae_2010_Ce_stiniae	1	2
# Cephalanthera_damasonium_MH590345	0	1
# Cephalanthera_falcata_OR161968	0	1
# Cephalanthera_humilis_KU551265	1	2
# Cephalanthera_longibracteata_MH590346	0	1
# Cephalanthera_longifolia_KU551263	0	1
# Cephalanthera_nanchuanica_OR947462	0	1
# Cephalanthera_rubra_MH590347	0	1
# Diplandrorchis_sinica_MZ014629	1	2
# Diplandrorchis_sinica_OP310037	1	2
# Epipactis_albensis_MH590348	0	1
# Epipactis_atrorubens_MH590349	0	1
# Epipactis_gigantea_MH590350	0	0
# Epipactis_helleborine_MH590351	0	1
# Epipactis_helleborine_MK593537	0	1
# Epipactis_helleborine_MK608776	0	1
# Epipactis_mairei_KU551264	0	0
# Epipactis_mairei_MG925367	0	0
# Epipactis_microphylla_MH590352	0	1
# Epipactis_palustris_MH590353	0	0
# Epipactis_purpurata_MH590354	0	1
# Epipactis_thunbergii_MN200387	0	0
# Epipactis_veratrifolia_KU551267	0	0
# Limodorum_abortivum_MH590355	1	1
# Neottia_acuminata_KU551268	1	2
# Neottia_camtschatea_KU551266	1	2
# Neottia_cordata_MH590356	0	1
# Neottia_fugongensis_KU551270	0	0
# Neottia_japonica_MH321183	0	0
# Neottia_japonica_MH321184	0	0
# Neottia_listeroides_KU551272	1	2
# Neottia_nidus_avis_JF325876	1	2
# Neottia_ovata_KU551271	0	1
# Neottia_pinetorum_KU551269	0	0
# Neottia_suzukii_MH321185	0	0
# Nervilia_fordii_ON515491	0	0
# Palmorchis_pabstii_MH590357	0	0
# Sobralia_callosa_KM032623	0	0
# Spiranthes_sinensis_MK936427	0	0



# Load libraries
library(ape)
library(phytools)

# 1. Read and root the tree
tree <- read.tree("mixfinder.contree")
tree <- root(tree, outgroup = "Spiranthes_sinensis_MK936427", resolve.root = TRUE)
tree <- ladderize(tree)

# 2. Make tree ultrametric with chronos
tree <- chronos(tree, lambda = 1)  # Use λ=1 for moderate smoothing
if (!inherits(tree, "phylo")) stop("chronos() failed and returned a non-phylo object")

# 3. Read in tip data
tipdata <- read.csv("tipdata.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(tipdata) <- tipdata$species

# 4. Match tree and trait data
common_tips <- intersect(tree$tip.label, rownames(tipdata))
tree <- drop.tip(tree, setdiff(tree$tip.label, common_tips))
tipdata <- tipdata[tree$tip.label, ]  # reorder to match tree

# 5. Format traits as named vectors
leafy <- setNames(as.character(tipdata$leafy), rownames(tipdata))
myco  <- setNames(as.character(tipdata$myco), rownames(tipdata))

# 6. Define irreversible transition matrices

# Only allowed:			leafy to leafless (0->1)

# Only Allowed:			AU to PM (0->1); PM->FM (1->2)

Q_leafy <- matrix(c(0, 1,
                    0, 0), nrow = 2, byrow = TRUE,
                  dimnames = list(c("0", "1"), c("0", "1")))

Q_myco <- matrix(c(0, 1, 0,
                   0, 0, 1,
                   0, 0, 0), nrow = 3, byrow = TRUE,
                 dimnames = list(c("0", "1", "2"), c("0", "1", "2")))

# 7. Run stochastic character mapping
set.seed(123)

sim_leafy <- make.simmap(tree, leafy, model = Q_leafy, nsim = 100,
                         rate = "fixed", pi = c("0" = 1, "1" = 0))

sim_myco <- make.simmap(tree, myco, model = Q_myco, nsim = 100,
                        rate = "fixed", pi = c("0" = 1, "1" = 0, "2" = 0))

summary_leafy <- summary(sim_leafy)
summary_myco  <- summary(sim_myco)

# 8. Define color schemes
leafy_colors <- c("0" = "blue", "1" = "red")  # 0 = leafy, 1 = leafless
myco_colors  <- c("0" = "blue", "1" = "orange", "2" = "red")  # 0 = AU, 1 = PM, 2 = FM

# 9. Plot posterior summaries from all simulations
par(mfrow = c(2, 1))

# Leafiness
plot(summary_leafy, colors = leafy_colors, fsize = 1,
     ftype = "i", xlim = c(0, 1.2 * max(nodeHeights(tree))))
title("Leafiness (0 = leafy, 1 = leafless)")
add.simmap.legend(colors = leafy_colors, prompt = TRUE, x = 0.01, y = Ntip(tree) * 0.95)

# Mycorrhizal mode
plot(summary_myco, colors = myco_colors, fsize = 1,
     ftype = "i", xlim = c(0, 1.2 * max(nodeHeights(tree))))
title("Mycorrhizal Mode (0 = AU, 1 = PM, 2 = FM)")
add.simmap.legend(colors = myco_colors, prompt = TRUE, x = 0.01, y = Ntip(tree) * 0.95)


