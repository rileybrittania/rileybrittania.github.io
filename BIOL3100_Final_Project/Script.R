library(bio3d)
library(tidyverse)
library(easystats)
library(skimr)
library(MASS)
protein_df <- read.csv('Protein_Structures_Comparison.csv') %>% 
  mutate(Method = trimws(Method))

##3-to-1 AA Conversion (for pasting into AlphaFold)
###Step-by-Step
pdb_1CA2 <- read.pdb('PDBs/1CA2.pdb')
aa1_1CA2 <- aa321(pdb_1CA2$seqres)
aa1_string_1CA2 <- paste(aa1_1CA2, collapse = "")
print(aa1_string_1CA2)
##Using Tidyverse, 1 file at a time
aa1_1CA2 <- read.pdb('PDBs/1CA2.pdb')$seqres %>%
  aa321() %>%
  paste(collapse = "") %>%
  print()
##Using Tidyverse, all files (WHAT I DID)
pdb_files <- list.files(path = "PDBs/", pattern = "\\.pdb$", full.names = TRUE)
all_sequences <- pdb_files %>%
  set_names(str_remove(basename(.), "\\.pdb$")) %>%
  map(~ read.pdb(.)$seqres %>% 
        aa321() %>% 
        paste(collapse = ""))
print(all_sequences)

##Got the number of amino acids per sequence through Alpha Fold 3, code too complex
##Used pymol to superimpose proteins and obtain RMSD values
###using "align to molecule" function
###removed waters from x-ray files, as well as other atoms/DNA/ligands
##Look at screenshots for pymol code 
##create and save images of overlaps for website (maybe only do a few, all if time)
##Save everything into objects before transferring to markdown!

##Plot for Assignment 4, AA Sequence Length vs X-Ray Diffraction Resolution
protein_df <- read.csv('Protein_Structures_Comparison.csv')
View(protein_df)
unique(protein_df$PDB_Class)
names(protein_df)
library(ggpubr)

assignment_4_plot <- protein_df %>% 
  filter(!is.na(Resolution)) %>% 
  ggplot(aes(x = Number_of_AA,
             y = Resolution)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F, color = 'orange') +
  labs(x = "Number of Amino Acids",
       y = "Resolution (in Å)",
       title = "Number of Amino Acids vs X-Ray Diffraction Resolution") +
  stat_regline_equation(label.x.npc = "left", label.y.npc = "top") +
  stat_cor(aes(label = paste(after_stat(rr.label))),
           label.x.npc = "left", label.y.npc = 0.85, hjust = 0)
ggsave('Assignment_4_Plot.png', plot = assignment_4_plot,  
       dpi = 300, height = 5, width = 10)


###Sequence Length vs X-Ray Diffraction Resolution
protein_df %>% 
  filter(!is.na(Resolution)) %>% 
  ggplot(aes(x = Number_of_AA,
             y = Resolution)) +
  geom_point() +
  labs(x = "Number of Amino Acids",
       y = "Resolution (in Å)",
       title = "Sequence Length vs X-Ray Diffraction Resolution") 

###Sequence Length vs RMSD_2
protein_df %>% 
  filter(!is.na(RMSD_2)) %>% 
  ggplot(aes(x = Number_of_AA,
             y = RMSD_2)) +
  geom_point() +
  labs(x = "Number of Amino Acids",
       y = "RMSD from Alignment with AlphaFold2 (in Å)",
       title = "Sequence Length vs RMSD using AlphaFold2") 

###Sequence Length vs RMSD_3
protein_df %>% 
  filter(!is.na(RMSD_3)) %>% 
  ggplot(aes(x = Number_of_AA,
             y = RMSD_3)) +
  geom_point() +
  labs(x = "Number of Amino Acids",
       y = "RMSD from Alignment with AlphaFold3 (in Å)",
       title = "Sequence Length vs RMSD using AlphaFold3") 

###pTM vs RMSD_2
protein_df %>% 
  filter(!is.na(highest_pTM) & !is.na(RMSD_2)) %>% 
  ggplot(aes(x = highest_pTM,
             y = RMSD_2)) +
  geom_point() +
  labs(x = "pTM Score",
       y = "RMSD from Alignment with AlphaFold2 (in Å)",
       title = "pTM from Highest Rank vs RMSD using AlphaFold2")

###pTM vs RMSD_3
protein_df %>% 
  filter(!is.na(highest_pTM) & !is.na(RMSD_3)) %>% 
  ggplot(aes(x = highest_pTM,
             y = RMSD_3)) +
  geom_point() +
  labs(x = "pTM Score",
       y = "RMSD from Alignment with AlphaFold3 (in Å)",
       title = "pTM from Highest Rank vs RMSD using AlphaFold3")

###Protein vs Sequence Length (colored by PDB Class)
protein_df %>% 
  ggplot(aes(x = Protein_Name,
             y = Number_of_AA, 
             fill = PDB_Class)) +
  geom_bar(stat = 'identity') +
  labs(x = "Protein_Name",
       y = "Number of Amino Acids",
       title = "Protein vs Sequence Length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

###pTM vs Sequence Length (separated by method)
protein_df %>% 
  filter(!is.na(highest_pTM)) %>% 
  ggplot(aes(x = Number_of_AA,
             y = highest_pTM, 
             color = Method)) +
  geom_point() +
  labs(x = "Number of Amino Acids",
       y = "pTM Score",
       title = "Sequence Length vs pTM from Highest Rank")

##STATISTICAL ANALYSES
###Linear Models
mod_aa_xray <- glm(dat = protein_df,
                formula = Resolution ~ Number_of_AA)
summary(mod_aa_xray)

mod_aa_rmsd2 <- glm(dat = protein_df,
                formula = RMSD_2 ~ Number_of_AA)
summary(mod_aa_rmsd2)

mod_aa_rmsd3 <- glm(dat = protein_df,
                 formula = RMSD_3 ~ Number_of_AA)
summary(mod_aa_rmsd3)

mod_ptm_rmsd2 <- glm(dat = protein_df,
                     formula = RMSD_2 ~ highest_pTM)
summary(mod_ptm_rmsd2)

mod_ptm_rmsd3 <- glm(dat = protein_df,
                     formula = RMSD_3 ~ highest_pTM)
summary(mod_ptm_rmsd3)

mod_aa_ptm <- glm(dat = protein_df,
                  formula = Number_of_AA ~ highest_pTM)

###Comparison of Linear Models
compare_performance(mod_aa_xray, mod_aa_rmsd2, mod_aa_rmsd3, 
                    mod_ptm_rmsd2, mod_ptm_rmsd3, mod_aa_ptm) %>% plot()

### Finding the best models for RMSD2 and RMSD 3 with MASS
mod_all_rmsd2 <- glm(dat = protein_df,
                     formula = RMSD_2 ~ Number_of_AA * highest_pTM)
step_rmsd2 = stepAIC(mod_all_rmsd2)
step_rmsd2$formula

mod_all_rmsd3 <- glm(dat = protein_df,
                     formula = RMSD_3 ~ Number_of_AA * highest_pTM)
step_rmsd3 = stepAIC(mod_all_rmsd3)
step_rmsd3$formula

### Comparison of Linear Models With Best Models
compare_performance(mod_aa_xray, mod_aa_rmsd2, mod_aa_rmsd3, 
                    mod_ptm_rmsd2, mod_ptm_rmsd3, mod_aa_ptm,
                    mod_xray_rmsd2, mod_xray_rmsd3,
                    mod_all_rmsd2, mod_all_rmsd3) %>% plot()

### Average RMSD & pTM per Method
Average_RMSD_2 <- protein_df %>% 
  filter(!is.na(RMSD_2)) %>% 
  group_by(Method) %>%
  summarize(Mean_RMSD_AlphaFold2 = mean(RMSD_2))
Average_RMSD_2

Average_RMSD_3 <- protein_df %>% 
  filter(!is.na(RMSD_3)) %>% 
  group_by(Method) %>%
  summarize(Mean_RMSD_AlphaFold3 = mean(RMSD_3))
Average_RMSD_3

Average_pTM <- protein_df %>% 
  filter(!is.na(highest_pTM)) %>% 
  group_by(Method) %>%
  summarize(Mean_pTM = mean(highest_pTM))
Average_pTM
