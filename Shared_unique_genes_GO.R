library(ggvenn)
library(VennDiagram)
library(openxlsx)
library(topGO)
library(gprofiler2)
library("qpcR")    
library(openxlsx)
library("readxl")

all=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_RBE.xlsx', sheet='Genes'))
logFC=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_RBE.xlsx', sheet='LogFC'))
p_adj=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_RBE.xlsx', sheet='p_adj'))

d_m=na.omit(all['Male.Down'])
d_f=na.omit(all['Female.Down'])
u_m=na.omit(all['Male.Up'])
u_f=na.omit(all['Female.Up'])

d_a=na.omit(all['Young.Down'])
d_oa=na.omit(all['Old.Down'])
u_a=na.omit(all['Young.Up'])
u_oa=na.omit(all['Old.Up'])

a <- list('Young' = c(d_a[[1]],u_a[[1]]),
          'Young' = c(d_oa[[1]],u_oa[[1]]))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Genes at T0")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_Young_Adult_T0.jpg");
dev.off ()

a <- list('Male' = c(d_m[[1]],u_m[[1]]),
          'Female' = c(d_f[[1]],u_f[[1]]))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Genes at T0")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_male_female_T0.jpg");
dev.off ()


GO_d_m <- gost(query = d_m[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_f <- gost(query = d_f[[1]], 
                organism = "hsapiens", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_m <- gost(query = u_m[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_f <- gost(query = u_f[[1]], 
                organism = "hsapiens", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_a <- gost(query = d_a[[1]], 
                           organism = "hsapiens", ordered_query = TRUE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_oa <- gost(query = d_oa[[1]], 
                           organism = "hsapiens", ordered_query = TRUE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_a <- gost(query = u_a[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_oa <- gost(query = u_oa[[1]], 
                organism = "hsapiens", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

# intersections
shared_down_f_m=intersect(d_f[[1]], d_m[[1]])
ind_d_m=match(shared_down_f_m,d_m[[1]])
ind_d_f=match(shared_down_f_m,d_f[[1]])

GO_shared_down_f_m <- gost(query = shared_down_f_m, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('DownReg Female' = d_f[[1]],
          'DownReg Male' = d_m[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_down_f_m.jpg");
dev.off ()

shared_up_f_m=intersect(u_f[[1]], u_m[[1]])
ind_u_m=match(shared_up_f_m,u_m[[1]])
ind_u_f=match(shared_up_f_m, u_f[[1]])
GO_shared_up_f_m <- gost(query = shared_up_f_m, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('UpReg Female' = u_f[[1]],
          'UpReg Male' = u_m[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_up_f_m.jpg");
dev.off ()

shared_down_a_oa=intersect(d_a[[1]], d_oa[[1]])
ind_d_a=match(shared_down_a_oa, d_a[[1]])
ind_d_oa=match(shared_down_a_oa,d_oa[[1]])
GO_shared_down_a_oa <- gost(query = shared_down_a_oa, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('DownReg Young' = d_a[[1]],
          'DownReg Old' = d_oa[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_down_a_oa.jpg");
dev.off ()

shared_up_a_oa=intersect(u_a[[1]], u_oa[[1]])
ind_u_a=match(shared_up_a_oa, u_a[[1]])
ind_u_oa=match(shared_up_a_oa, u_oa[[1]])
GO_shared_up_a_oa <- gost(query = shared_up_a_oa, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('UpReg Young' = u_a[[1]],
          'UpReg Old' = u_oa[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_up_a_oa.jpg");
dev.off ()

# Non-shared values 
unique_down_m <- setdiff(d_m[[1]],d_f[[1]])
unique_down_f <- setdiff(d_f[[1]],d_m[[1]])
unique_up_m <- setdiff(u_m[[1]],u_f[[1]])
unique_up_f <- setdiff(u_f[[1]],u_m[[1]])
unique_down_a <- setdiff(d_a[[1]],d_oa[[1]])
unique_down_oa <- setdiff(d_oa[[1]],d_a[[1]])
unique_up_a <- setdiff(u_a[[1]],u_oa[[1]])
unique_up_oa <- setdiff(u_oa[[1]],u_a[[1]])

ind_unique_d_m=match(unique_down_m, d_m[[1]])
ind_unique_d_f=match(unique_down_f, d_f[[1]])
ind_unique_u_m=match(unique_up_m, u_m[[1]])
ind_unique_u_f=match(unique_up_f, u_f[[1]])

ind_unique_d_a=match(unique_down_a, d_a[[1]])
ind_unique_d_oa=match(unique_down_oa, d_oa[[1]])
ind_unique_u_a=match(unique_up_a, u_a[[1]])
ind_unique_u_oa=match(unique_up_oa, u_oa[[1]])


GO_unique_down_m <- gost(query = unique_down_m, 
                        organism = "hsapiens", ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_down_f <- gost(query = unique_down_f, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Down_f_m=intersect(GO_unique_down_m$result$term_id,GO_unique_down_f$result$term_id)
term_ids <- GO_unique_down_f[["result"]][["term_id"]]
matching_indices <- which(Shared_G0ID_Down_f_m %in% term_ids)
Shared_G0name_Down_f_m=intersect(GO_unique_down_m$result$term_name,GO_unique_down_f$result$term_name)
Shared_G0name_Down_f_m=GO_unique_down_f$result$term_name[matching_indices]

a <- list('DownReg Female' = GO_unique_down_m$result$term_id,
          'DownReg Male' = GO_unique_down_f$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/Shared_G0ID_Down_f_m.jpg");
dev.off ()

GO_unique_up_m <- gost(query = unique_up_m, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_up_f <- gost(query = unique_up_f, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Up_f_m=intersect(GO_unique_up_m$result$term_id,GO_unique_up_f$result$term_id)
term_ids <- GO_unique_up_f[["result"]][["term_id"]]
matching_indices <- which(Shared_G0ID_Up_f_m %in% term_ids)
Shared_G0name_Up_f_m=intersect(GO_unique_up_m$result$term_name,GO_unique_up_f$result$term_name)
Shared_G0name_Up_f_m=GO_unique_up_f$result$term_name[matching_indices]


a <- list('UpReg Female' = GO_unique_up_m$result$term_id,
          'UpReg Male' = GO_unique_up_f$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/Shared_G0ID_Up_f_m.jpg");
dev.off ()

GO_unique_down_a <- gost(query = unique_down_a, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_down_oa <- gost(query = unique_down_oa, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Down_a_oa=intersect(GO_unique_down_a$result$term_id,GO_unique_down_oa$result$term_id)
Shared_G0name_Down_a_oa=intersect(GO_unique_down_a$result$term_name,GO_unique_down_oa$result$term_name)
a <- list('DownReg Young' = GO_unique_down_a$result$term_id,
          'DownReg Old' = GO_unique_down_oa$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/Shared_G0ID_Down_a_oa.jpg");
dev.off ()

GO_unique_up_a <- gost(query = unique_up_a, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_up_oa <- gost(query = unique_up_oa, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Up_a_oa=intersect(GO_unique_up_a$result$term_id,GO_unique_up_oa$result$term_id)
Shared_G0name_Up_a_oa=intersect(GO_unique_up_a$result$term_name,GO_unique_up_oa$result$term_name)
a <- list('UpReg Young' = GO_unique_up_a$result$term_id,
          'UpReg Old' = GO_unique_up_oa$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/Shared_G0ID_Up_a_oa.jpg");
dev.off ()

# Create an Excel file and add sheets with the lists
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_GO_T0_RBE.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Shared_G0_Down_f_m')
writeData(wb, sheet = 'Shared_G0_Down_f_m', x = cbind(Shared_G0ID_Down_f_m,Shared_G0name_Down_f_m))
addWorksheet(wb, sheetName = 'Shared_GO_Up_f_m')
writeData(wb, sheet = 'Shared_GO_Up_f_m', x = cbind(Shared_G0ID_Up_f_m,Shared_G0name_Up_f_m))
addWorksheet(wb, sheetName = 'Shared_G0_Down_a_oa')
writeData(wb, sheet = 'Shared_G0_Down_a_oa', x = cbind(Shared_G0ID_Down_a_oa,Shared_G0name_Down_a_oa))
addWorksheet(wb, sheetName = 'Shared_G0_Up_a_oa')
writeData(wb, sheet = 'Shared_G0_Up_a_oa', x = cbind(Shared_G0ID_Up_a_oa,Shared_G0name_Up_a_oa))
#l_d_a=logFC[ind_d_a,'Young.Down']
addWorksheet(wb, sheetName = 'Shar_Gene_Down_f_m')
writeData(wb, sheet = 'Shar_Gene_Down_f_m', x = cbind(GO_shared_down_f_m$result$term_id,GO_shared_down_f_m$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Up_f_m')
writeData(wb, sheet = 'Shar_Gene_Up_f_m', x = cbind(GO_shared_up_f_m$result$term_id,GO_shared_up_f_m$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Down_a_oa')
writeData(wb, sheet = 'Shar_Gene_Down_a_oa', x = cbind(GO_shared_down_a_oa$result$term_id,GO_shared_down_a_oa$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Up_a_oa')
writeData(wb, sheet = 'Shar_Gene_Up_a_oa', x = cbind(GO_shared_up_a_oa$result$term_id,GO_shared_up_a_oa$result$term_name))

addWorksheet(wb, sheetName = 'G0_Down_a')
writeData(wb, sheet = 'G0_Down_a', x = cbind(GO_d_a$result$term_id,GO_d_a$result$term_name))
addWorksheet(wb, sheetName = 'G0_Down_oa')
writeData(wb, sheet = 'G0_Down_oa', x = cbind(GO_d_oa$result$term_id,GO_d_oa$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_a')
writeData(wb, sheet = 'G0_Up_a', x = cbind(GO_u_a$result$term_id,GO_u_a$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_oa')
writeData(wb, sheet = 'G0_Up_oa', x = cbind(GO_u_oa$result$term_id,GO_u_oa$result$term_name))

addWorksheet(wb, sheetName = 'G0_Down_f')
writeData(wb, sheet = 'G0_Down_f', x = cbind(GO_d_f$result$term_id,GO_d_f$result$term_name))
addWorksheet(wb, sheetName = 'G0_Down_m')
writeData(wb, sheet = 'G0_Down_m', x = cbind(GO_d_m$result$term_id,GO_d_m$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_f')
writeData(wb, sheet = 'G0_Up_f', x = cbind(GO_u_f$result$term_id,GO_u_f$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_m')
writeData(wb, sheet = 'G0_Up_m', x = cbind(GO_u_m$result$term_id,GO_u_m$result$term_name))

#GENES 

Shared_down_up_m <- qpcR:::cbind.na(shared_down_f_m,logFC[ind_d_m,'Male.Down'],p_adj[ind_d_m,'Male.Down'],shared_up_f_m,logFC[ind_u_m,'Male.Up'],p_adj[ind_u_m,'Male.Up']) # Bind as columns
Shared_down_up_f <- qpcR:::cbind.na(shared_down_f_m,logFC[ind_d_f,'Female.Down'],p_adj[ind_d_f,'Female.Down'],shared_up_f_m,logFC[ind_u_f,'Female.Up'],p_adj[ind_u_f,'Female.Up']) # Bind as columns
Shared_down_up_a <- qpcR:::cbind.na(shared_down_a_oa,logFC[ind_d_a,'Young.Down'],p_adj[ind_d_a,'Young.Down'],shared_up_a_oa,logFC[ind_u_a,'Young.Up'],p_adj[ind_u_a,'Young.Up'])
Shared_down_up_oa <- qpcR:::cbind.na(shared_down_a_oa,logFC[ind_d_oa,'Old.Down'],p_adj[ind_d_oa,'Old.Down'],shared_up_a_oa,logFC[ind_u_oa,'Old.Up'],p_adj[ind_u_oa,'Old.Up'])

Unique_down_up_m <- qpcR:::cbind.na(unique_down_m,logFC[ind_unique_d_m,'Male.Down'],p_adj[ind_unique_d_m,'Male.Down'],unique_up_m,logFC[ind_unique_u_m,'Male.Up'],p_adj[ind_unique_u_m,'Male.Up'])
Unique_down_up_f <- qpcR:::cbind.na(unique_down_f,logFC[ind_unique_d_f,'Female.Down'],p_adj[ind_unique_d_f,'Female.Down'],unique_up_f,logFC[ind_unique_u_f,'Female.Up'],p_adj[ind_unique_u_f,'Female.Up'])
Unique_down_up_a <- qpcR:::cbind.na(unique_down_a,logFC[ind_unique_d_a,'Young.Down'],p_adj[ind_unique_d_a,'Young.Down'],unique_up_a,logFC[ind_unique_u_a,'Young.Up'],p_adj[ind_unique_u_a,'Young.Up'])
Unique_down_up_oa <- qpcR:::cbind.na(unique_down_oa,logFC[ind_unique_d_oa,'Old.Down'],p_adj[ind_unique_d_oa,'Old.Down'],unique_up_oa,logFC[ind_unique_u_oa,'Old.Up'],p_adj[ind_unique_u_oa,'Old.Up'])

addWorksheet(wb, sheetName = 'Shared_down_up_m')
writeData(wb, sheet = 'Shared_down_up_m', x = Shared_down_up_m)
addWorksheet(wb, sheetName = 'Shared_down_up_f')
writeData(wb, sheet = 'Shared_down_up_f', x = Shared_down_up_f)
addWorksheet(wb, sheetName = 'Shared_down_up_a')
writeData(wb, sheet = 'Shared_down_up_a', x = Shared_down_up_a)
addWorksheet(wb, sheetName = 'Shared_down_up_oa')
writeData(wb, sheet = 'Shared_down_up_oa', x = Shared_down_up_oa)

addWorksheet(wb, sheetName = 'Unique_down_up_m')
writeData(wb, sheet = 'Unique_down_up_m', x = Unique_down_up_m)
addWorksheet(wb, sheetName = 'Unique_down_up_f')
writeData(wb, sheet = 'Unique_down_up_f', x = Unique_down_up_f)
addWorksheet(wb, sheetName = 'Unique_down_up_a')
writeData(wb, sheet = 'Unique_down_up_a', x = Unique_down_up_a)
addWorksheet(wb, sheetName = 'Unique_down_up_oa')
writeData(wb, sheet = 'Unique_down_up_oa', x = Unique_down_up_oa)

# Save the Excel file
#output_file="Up_down_categories_GO_gender_T0_adj_age.xlsx"
#output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_GO_gender_T0.xlsx'
saveWorkbook(wb, output_file)

## T1

all=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_T1_RBE.xlsx', sheet='Genes'))
logFC=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_T1_RBE.xlsx', sheet='LogFC'))
p=data.frame(read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_T1_RBE.xlsx', sheet='p_adj'))

d_m=na.omit(all['Male.Down'])
d_f=na.omit(all['Female.Down'])
u_m=na.omit(all['Male.Up'])
u_f=na.omit(all['Female.Up'])
d_a=na.omit(all['Young.Down'])
d_oa=na.omit(all['Old.Down'])
u_a=na.omit(all['Young.Up'])
u_oa=na.omit(all['Old.Up'])


a <- list('Young' = c(d_a[[1]],u_a[[1]]),
          'Old' = c(d_oa[[1]],u_oa[[1]]))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Genes at T1")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_Young_Adult_T1.jpg");
dev.off ()

a <- list('Male' = c(d_m[[1]],u_m[[1]]),
          'Female' = c(d_f[[1]],u_f[[1]]))
venn <- ggvenn(a)
# Add a title to the Venn diagram
venn <- venn + ggtitle("DE Genes at T1")
# Print the Venn diagram
print(venn)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_male_female_T1.jpg");
dev.off ()


GO_d_m <- gost(query = d_m[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_f <- gost(query = d_f[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_m <- gost(query = u_m[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_f <- gost(query = u_f[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_a <- gost(query = d_a[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_d_oa <- gost(query = d_oa[[1]], 
                organism = "hsapiens", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_a <- gost(query = u_a[[1]], 
               organism = "hsapiens", ordered_query = TRUE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "annotated", custom_bg = NULL, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

GO_u_oa <- gost(query = u_oa[[1]], 
                organism = "hsapiens", ordered_query = TRUE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

# intersections
shared_down_f_m=intersect(d_f[[1]], d_m[[1]])
ind_d_m=match(shared_down_f_m,d_m[[1]])
ind_d_f=match(shared_down_f_m,d_f[[1]])

GO_shared_down_f_m <- gost(query = shared_down_f_m, 
                           organism = "hsapiens", ordered_query = TRUE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = "fdr", 
                           domain_scope = "annotated", custom_bg = NULL, 
                           numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('DownReg Female' = d_f[[1]],
          'DownReg Male' = d_m[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_down_f_m_T1.jpg");
dev.off ()

shared_up_f_m=intersect(u_f[[1]], u_m[[1]])
ind_u_m=match(shared_up_f_m,u_m[[1]])
ind_u_f=match(shared_up_f_m, u_f[[1]])
GO_shared_up_f_m <- gost(query = shared_up_f_m, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('UpReg Female' = u_f[[1]],
          'UpReg Male' = u_m[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/shared_up_f_m_T1.jpg");
dev.off ()

shared_down_a_oa=intersect(d_a[[1]], d_oa[[1]])
ind_d_a=match(shared_down_a_oa, d_a[[1]])
ind_d_oa=match(shared_down_a_oa,d_oa[[1]])
GO_shared_down_a_oa <- gost(query = shared_down_a_oa, 
                            organism = "hsapiens", ordered_query = TRUE, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                            measure_underrepresentation = FALSE, evcodes = FALSE, 
                            user_threshold = 0.05, correction_method = "fdr", 
                            domain_scope = "annotated", custom_bg = NULL, 
                            numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('DownReg Young' = d_a[[1]],
          'DownReg Old' = d_oa[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_down_a_oa_T1.jpg");
dev.off ()

shared_up_a_oa=intersect(u_a[[1]], u_oa[[1]])
ind_u_a=match(shared_up_a_oa, u_a[[1]])
ind_u_oa=match(shared_up_a_oa, u_oa[[1]])
GO_shared_up_a_oa <- gost(query = shared_up_a_oa, 
                          organism = "hsapiens", ordered_query = TRUE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = FALSE, 
                          user_threshold = 0.05, correction_method = "fdr", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

a <- list('UpReg Young' = u_a[[1]],
          'UpReg Old' = u_oa[[1]])
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/shared_up_a_oa_T1.jpg");
dev.off ()

# Non-shared values 
unique_down_m <- setdiff(d_m[[1]],d_f[[1]])
unique_down_f <- setdiff(d_f[[1]],d_m[[1]])
unique_up_m <- setdiff(u_m[[1]],u_f[[1]])
unique_up_f <- setdiff(u_f[[1]],u_m[[1]])
unique_down_a <- setdiff(d_a[[1]],d_oa[[1]])
unique_down_oa <- setdiff(d_oa[[1]],d_a[[1]])
unique_up_a <- setdiff(u_a[[1]],u_oa[[1]])
unique_up_oa <- setdiff(u_oa[[1]],u_a[[1]])

ind_unique_d_m=match(unique_down_m, d_m[[1]])
ind_unique_d_f=match(unique_down_f, d_f[[1]])
ind_unique_u_m=match(unique_up_m, u_m[[1]])
ind_unique_u_f=match(unique_up_f, u_f[[1]])

ind_unique_d_a=match(unique_down_a, d_a[[1]])
ind_unique_d_oa=match(unique_down_oa, d_oa[[1]])
ind_unique_u_a=match(unique_up_a, u_a[[1]])
ind_unique_u_oa=match(unique_up_oa, u_oa[[1]])


GO_unique_down_m <- gost(query = unique_down_m, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_down_f <- gost(query = unique_down_f, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Down_f_m=intersect(GO_unique_down_m$result$term_id,GO_unique_down_f$result$term_id)
term_ids <- GO_unique_down_f[["result"]][["term_id"]]
matching_indices <- which(Shared_G0ID_Down_f_m %in% term_ids)
Shared_G0name_Down_f_m=intersect(GO_unique_down_m$result$term_name,GO_unique_down_f$result$term_name)
Shared_G0name_Down_f_m=GO_unique_down_f$result$term_name[matching_indices]

a <- list('DownReg Female' = GO_unique_down_m$result$term_id,
          'DownReg Male' = GO_unique_down_f$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/Shared_G0ID_Down_f_m_T1.jpg");
dev.off ()

GO_unique_up_m <- gost(query = unique_up_m, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_up_f <- gost(query = unique_up_f, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Up_f_m=intersect(GO_unique_up_m$result$term_id,GO_unique_up_f$result$term_id)
term_ids <- GO_unique_up_f[["result"]][["term_id"]]
matching_indices <- which(Shared_G0ID_Up_f_m %in% term_ids)
Shared_G0name_Up_f_m=intersect(GO_unique_up_m$result$term_name,GO_unique_up_f$result$term_name)
Shared_G0name_Up_f_m=GO_unique_up_f$result$term_name[matching_indices]


a <- list('UpReg Female' = GO_unique_up_m$result$term_id,
          'UpReg Male' = GO_unique_up_f$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/GENDER/Shared_G0ID_Up_f_m_T1.jpg");
dev.off ()

GO_unique_down_a <- gost(query = unique_down_a, 
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = FALSE, 
                         user_threshold = 0.05, correction_method = "fdr", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_down_oa <- gost(query = unique_down_oa, 
                          organism = "hsapiens", ordered_query = TRUE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = FALSE, 
                          user_threshold = 0.05, correction_method = "fdr", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Down_a_oa=intersect(GO_unique_down_a$result$term_id,GO_unique_down_oa$result$term_id)
Shared_G0name_Down_a_oa=intersect(GO_unique_down_a$result$term_name,GO_unique_down_oa$result$term_name)
a <- list('DownReg Young' = GO_unique_down_a$result$term_id,
          'DownReg Old' = GO_unique_down_oa$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/Shared_G0ID_Down_a_oa_T1.jpg");
dev.off ()

GO_unique_up_a <- gost(query = unique_up_a, 
                       organism = "hsapiens", ordered_query = TRUE, 
                       multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                       measure_underrepresentation = FALSE, evcodes = FALSE, 
                       user_threshold = 0.05, correction_method = "fdr", 
                       domain_scope = "annotated", custom_bg = NULL, 
                       numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_up_oa <- gost(query = unique_up_oa, 
                        organism = "hsapiens", ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

Shared_G0ID_Up_a_oa=intersect(GO_unique_up_a$result$term_id,GO_unique_up_oa$result$term_id)
Shared_G0name_Up_a_oa=intersect(GO_unique_up_a$result$term_name,GO_unique_up_oa$result$term_name)
a <- list('UpReg Young' = GO_unique_up_a$result$term_id,
          'UpReg Old' = GO_unique_up_oa$result$term_id)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_Gender_comparison/Shared/AGE/Shared_G0ID_Up_a_oa_T1.jpg");
dev.off ()

# Create an Excel file and add sheets with the lists
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_categories_GO_T1_RBE.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Shared_G0_Down_f_m')
writeData(wb, sheet = 'Shared_G0_Down_f_m', x = cbind(Shared_G0ID_Down_f_m,Shared_G0name_Down_f_m))
addWorksheet(wb, sheetName = 'Shared_GO_Up_f_m')
writeData(wb, sheet = 'Shared_GO_Up_f_m', x = cbind(Shared_G0ID_Up_f_m,Shared_G0name_Up_f_m))
addWorksheet(wb, sheetName = 'Shared_G0_Down_a_oa')
writeData(wb, sheet = 'Shared_G0_Down_a_oa', x = cbind(Shared_G0ID_Down_a_oa,Shared_G0name_Down_a_oa))
addWorksheet(wb, sheetName = 'Shared_G0_Up_a_oa')
writeData(wb, sheet = 'Shared_G0_Up_a_oa', x = cbind(Shared_G0ID_Up_a_oa,Shared_G0name_Up_a_oa))
#l_d_a=logFC[ind_d_a,'Young.Down']
addWorksheet(wb, sheetName = 'Shar_Gene_Down_f_m')
writeData(wb, sheet = 'Shar_Gene_Down_f_m', x = cbind(GO_shared_down_f_m$result$term_id,GO_shared_down_f_m$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Up_f_m')
writeData(wb, sheet = 'Shar_Gene_Up_f_m', x = cbind(GO_shared_up_f_m$result$term_id,GO_shared_up_f_m$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Down_a_oa')
writeData(wb, sheet = 'Shar_Gene_Down_a_oa', x = cbind(GO_shared_down_a_oa$result$term_id,GO_shared_down_a_oa$result$term_name))
addWorksheet(wb, sheetName = 'Shar_Gene_Up_a_oa')
writeData(wb, sheet = 'Shar_Gene_Up_a_oa', x = cbind(GO_shared_up_a_oa$result$term_id,GO_shared_up_a_oa$result$term_name))

addWorksheet(wb, sheetName = 'G0_Down_a')
writeData(wb, sheet = 'G0_Down_a', x = cbind(GO_d_a$result$term_id,GO_d_a$result$term_name))
addWorksheet(wb, sheetName = 'G0_Down_oa')
writeData(wb, sheet = 'G0_Down_oa', x = cbind(GO_d_oa$result$term_id,GO_d_oa$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_a')
writeData(wb, sheet = 'G0_Up_a', x = cbind(GO_u_a$result$term_id,GO_u_a$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_oa')
writeData(wb, sheet = 'G0_Up_oa', x = cbind(GO_u_oa$result$term_id,GO_u_oa$result$term_name))

addWorksheet(wb, sheetName = 'G0_Down_f')
writeData(wb, sheet = 'G0_Down_f', x = cbind(GO_d_f$result$term_id,GO_d_f$result$term_name))
addWorksheet(wb, sheetName = 'G0_Down_m')
writeData(wb, sheet = 'G0_Down_m', x = cbind(GO_d_m$result$term_id,GO_d_m$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_f')
writeData(wb, sheet = 'G0_Up_f', x = cbind(GO_u_f$result$term_id,GO_u_f$result$term_name))
addWorksheet(wb, sheetName = 'G0_Up_m')
writeData(wb, sheet = 'G0_Up_m', x = cbind(GO_u_m$result$term_id,GO_u_m$result$term_name))

#GENES 

Shared_down_up_m <- qpcR:::cbind.na(shared_down_f_m,logFC[ind_d_m,'Male.Down'],p[ind_d_m,'Male.Down'],shared_up_f_m,logFC[ind_u_m,'Male.Up'],p[ind_u_m,'Male.Up']) # Bind as columns
Shared_down_up_f <- qpcR:::cbind.na(shared_down_f_m,logFC[ind_d_f,'Female.Down'],p[ind_d_f,'Female.Down'],shared_up_f_m,logFC[ind_u_f,'Female.Up'],p[ind_u_f,'Female.Up']) # Bind as columns
Shared_down_up_a <- qpcR:::cbind.na(shared_down_a_oa,logFC[ind_d_a,'Young.Down'],p[ind_d_a,'Young.Down'],shared_up_a_oa,logFC[ind_u_a,'Young.Up'],p[ind_u_a,'Young.Up'])
Shared_down_up_oa <- qpcR:::cbind.na(shared_down_a_oa,logFC[ind_d_oa,'Old.Down'],p[ind_d_oa,'Old.Down'],shared_up_a_oa,logFC[ind_u_oa,'Old.Up'],p[ind_u_oa,'Old.Up'])

Unique_down_up_m <- qpcR:::cbind.na(unique_down_m,logFC[ind_unique_d_m,'Male.Down'],p[ind_unique_d_m,'Male.Down'],unique_up_m,logFC[ind_unique_u_m,'Male.Up'],p[ind_unique_u_m,'Male.Up'])
Unique_down_up_f <- qpcR:::cbind.na(unique_down_f,logFC[ind_unique_d_f,'Female.Down'],p[ind_unique_d_f,'Female.Down'],unique_up_f,logFC[ind_unique_u_f,'Female.Up'],p[ind_unique_u_f,'Female.Up'])
Unique_down_up_a <- qpcR:::cbind.na(unique_down_a,logFC[ind_unique_d_a,'Young.Down'],p[ind_unique_d_a,'Young.Down'],unique_up_a,logFC[ind_unique_u_a,'Young.Up'],p[ind_unique_u_a,'Young.Up'])
Unique_down_up_oa <- qpcR:::cbind.na(unique_down_oa,logFC[ind_unique_d_oa,'Old.Down'],p[ind_unique_d_oa,'Old.Down'],unique_up_oa,logFC[ind_unique_u_oa,'Old.Up'],p[ind_unique_u_oa,'Old.Up'])

addWorksheet(wb, sheetName = 'Shared_down_up_m')
writeData(wb, sheet = 'Shared_down_up_m', x = Shared_down_up_m)
addWorksheet(wb, sheetName = 'Shared_down_up_f')
writeData(wb, sheet = 'Shared_down_up_f', x = Shared_down_up_f)
addWorksheet(wb, sheetName = 'Shared_down_up_a')
writeData(wb, sheet = 'Shared_down_up_a', x = Shared_down_up_a)
addWorksheet(wb, sheetName = 'Shared_down_up_oa')
writeData(wb, sheet = 'Shared_down_up_oa', x = Shared_down_up_oa)

addWorksheet(wb, sheetName = 'Unique_down_up_m')
writeData(wb, sheet = 'Unique_down_up_m', x = Unique_down_up_m)
addWorksheet(wb, sheetName = 'Unique_down_up_f')
writeData(wb, sheet = 'Unique_down_up_f', x = Unique_down_up_f)
addWorksheet(wb, sheetName = 'Unique_down_up_a')
writeData(wb, sheet = 'Unique_down_up_a', x = Unique_down_up_a)
addWorksheet(wb, sheetName = 'Unique_down_up_oa')
writeData(wb, sheet = 'Unique_down_up_oa', x = Unique_down_up_oa)

# Save the Excel file
saveWorkbook(wb, output_file)

#___________________Ages and Genders ________________________________________________________________________________


library(ggvenn)
library(VennDiagram)
library(openxlsx)
library(topGO)
library(gprofiler2)
library("qpcR") 
library("readxl") 
p_adj=0.05
logFC=0
transcriptomics_HS_T0 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Results_DE_analysis_genes.xlsx', sheet= 'Heat stress - Stroke T0')
colnames(transcriptomics_HS_T0)=c('SYMBOL','LogFC','FDR')
logFC_trans=transcriptomics_HS_T0$LogFC
FDR_trans=transcriptomics_HS_T0$FDR
proteomic_genes = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP.xlsx', sheet= 'Stress-T0 (Genes)')
proteomic_genes_HS_T0 <- subset(proteomic_genes, sapply(proteomic_genes$adj.P.Val, function(x) x < p_adj))
proteomic_proteins = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP.xlsx', sheet= 'Stress-T0 (Proteins)')
proteomic_proteins_HS_T0 <- subset(proteomic_proteins, sapply(proteomic_proteins$adj.P.Val, function(x) x < p_adj))
logFC_prot=proteomic_proteins_HS_T0$logFC
FDR_prot=proteomic_proteins_HS_T0$adj.P.Val
proteomic_genes_1 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP.xlsx', sheet= 'Stress-T1 (Genes)')
proteomic_genes_HS_T1 <- subset(proteomic_genes_1, sapply(proteomic_genes_1$adj.P.Val, function(x) x < p_adj))
proteomic_proteins_1 = read_excel('/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/TopGenes_DEP.xlsx', sheet= 'Stress-T1 (Proteins)')
proteomic_proteins_HS_T1 <- subset(proteomic_proteins_1, sapply(proteomic_proteins_1$adj.P.Val, function(x) x < p_adj))
logFC_prot_1=proteomic_proteins_HS_T1$logFC
FDR_prot_1=proteomic_proteins_HS_T1$adj.P.Val

up_transcriptomics=subset(transcriptomics_HS_T0$SYMBOL, sapply(transcriptomics_HS_T0$LogFC, function(x) x > 0))
down_transcriptomics=subset(transcriptomics_HS_T0$SYMBOL, sapply(transcriptomics_HS_T0$LogFC, function(x) x < -0))
up_proteomic_genes=subset(proteomic_genes_HS_T0$SYMBOL, sapply(proteomic_genes_HS_T0$logFC, function(x) x > logFC))
down_proteomic_genes=subset(proteomic_genes_HS_T0$SYMBOL, sapply(proteomic_genes_HS_T0$logFC, function(x) x < -logFC))
up_proteomic_proteins=subset(proteomic_proteins_HS_T0$SYMBOL, sapply(proteomic_proteins_HS_T0$logFC, function(x) x > logFC))
down_proteomic_proteins=subset(proteomic_proteins_HS_T0$SYMBOL, sapply(proteomic_proteins_HS_T0$logFC, function(x) x < -logFC))

up_proteomic_genes_HT1=subset(proteomic_genes_HS_T1$SYMBOL, sapply(proteomic_genes_HS_T1$logFC, function(x) x > logFC))
down_proteomic_genes_HT1=subset(proteomic_genes_HS_T1$SYMBOL, sapply(proteomic_genes_HS_T1$logFC, function(x) x < -logFC))
up_proteomic_proteins_HT1=subset(proteomic_proteins_HS_T1$SYMBOL, sapply(proteomic_proteins_HS_T1$logFC, function(x) x > logFC))
down_proteomic_proteins_HT1=subset(proteomic_proteins_HS_T1$SYMBOL, sapply(proteomic_proteins_HS_T1$logFC, function(x) x < -logFC))

#Check in proteins shared from T0 to T1
shared_T0_T1=intersect(up_proteomic_proteins_HT1,up_proteomic_proteins)
unique_T1 <- setdiff(up_proteomic_proteins_HT1,up_proteomic_proteins)
unique_T0 <- setdiff(up_proteomic_proteins,up_proteomic_proteins_HT1)
ind_unique_T1=match(unique_T1, proteomic_genes_HS_T1$SYMBOL)
ind_unique_T0=match(unique_T0, proteomic_genes_HS_T0$SYMBOL)
ind_shared_T0_T1_inT1=match(shared_T0_T1, proteomic_genes_HS_T1$SYMBOL)
ind_shared_T0_T1_inT0=match(shared_T0_T1, proteomic_genes_HS_T0$SYMBOL)

Shared_T0_T1_in_T0 <- data.frame(qpcR:::cbind.na(shared_T0_T1,logFC_prot[ind_shared_T0_T1_inT0],FDR_prot[ind_shared_T0_T1_inT0])) # Bind as columns
colnames(Shared_T0_T1_in_T0)=c('SYMBOLS','LogFC','FDR')
Shared_T0_T1_in_T1 <- data.frame(qpcR:::cbind.na(shared_T0_T1,logFC_prot_1[ind_shared_T0_T1_inT1],FDR_prot_1[ind_shared_T0_T1_inT1])) # Bind as columns
colnames(Shared_T0_T1_in_T1)=c('SYMBOLS','LogFC','FDR')
Unique_T0 <- data.frame(qpcR:::cbind.na(unique_T0,logFC_prot[ind_unique_T0],FDR_prot[ind_unique_T0])) # Bind as columns
colnames(Unique_T0)=c('SYMBOLS','LogFC','FDR')
Unique_T1 <- data.frame(qpcR:::cbind.na(unique_T1,logFC_prot_1[ind_unique_T1],FDR_prot_1[ind_unique_T1])) # Bind as columns
colnames(Unique_T1)=c('SYMBOLS','LogFC','FDR')
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Shared_T0_T1_up.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Shared_in_T0')
writeData(wb, sheet = 'Shared_in_T0', x = Shared_T0_T1_in_T0)
addWorksheet(wb, sheetName = 'Shared_in_T1')
writeData(wb, sheet = 'Shared_in_T1', x = Shared_T0_T1_in_T1)
addWorksheet(wb, sheetName = 'Unique_T0')
writeData(wb, sheet = 'Unique_T0', x = Unique_T0)
addWorksheet(wb, sheetName = 'Unique_T1')
writeData(wb, sheet = 'Unique_T1', x = Unique_T1)
saveWorkbook(wb, output_file)

#Check in proteins shared from T0 to T1
shared_T0_T1=intersect(down_proteomic_proteins_HT1,down_proteomic_proteins)
unique_T1 <- setdiff(down_proteomic_proteins_HT1,down_proteomic_proteins)
unique_T0 <- setdiff(down_proteomic_proteins,down_proteomic_proteins_HT1)
ind_unique_T1=match(unique_T1, proteomic_genes_HS_T1$SYMBOL)
ind_unique_T0=match(unique_T0, proteomic_genes_HS_T0$SYMBOL)
ind_shared_T0_T1_inT1=match(shared_T0_T1, proteomic_genes_HS_T1$SYMBOL)
ind_shared_T0_T1_inT0=match(shared_T0_T1, proteomic_genes_HS_T0$SYMBOL)

Shared_T0_T1_in_T0 <- data.frame(qpcR:::cbind.na(shared_T0_T1,logFC_prot[ind_shared_T0_T1_inT0],FDR_prot[ind_shared_T0_T1_inT0])) # Bind as columns
colnames(Shared_T0_T1_in_T0)=c('SYMBOLS','LogFC','FDR')
Shared_T0_T1_in_T1 <- data.frame(qpcR:::cbind.na(shared_T0_T1,logFC_prot_1[ind_shared_T0_T1_inT1],FDR_prot_1[ind_shared_T0_T1_inT1])) # Bind as columns
colnames(Shared_T0_T1_in_T1)=c('SYMBOLS','LogFC','FDR')
Unique_T0 <- data.frame(qpcR:::cbind.na(unique_T0,logFC_prot[ind_unique_T0],FDR_prot[ind_unique_T0])) # Bind as columns
colnames(Unique_T0)=c('SYMBOLS','LogFC','FDR')
Unique_T1 <- data.frame(qpcR:::cbind.na(unique_T1,logFC_prot_1[ind_unique_T1],FDR_prot_1[ind_unique_T1])) # Bind as columns
colnames(Unique_T1)=c('SYMBOLS','LogFC','FDR')
output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Proteomics_C_data/Shared_T0_T1_down.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Shared_in_T0')
writeData(wb, sheet = 'Shared_in_T0', x = Shared_T0_T1_in_T0)
addWorksheet(wb, sheetName = 'Shared_in_T1')
writeData(wb, sheet = 'Shared_in_T1', x = Shared_T0_T1_in_T1)
addWorksheet(wb, sheetName = 'Unique_T0')
writeData(wb, sheet = 'Unique_T0', x = Unique_T0)
addWorksheet(wb, sheetName = 'Unique_T1')
writeData(wb, sheet = 'Unique_T1', x = Unique_T1)
saveWorkbook(wb, output_file)

# # CHECK DIRECTION

# proteomic_genes=proteomic_genes_HS_T0[,1:2]
# colnames(proteomic_genes)=c('SYMBOL','prot_logFC')
# transcriptomics=transcriptomics_HS_T0[,1:2]
# colnames(transcriptomics)=c('SYMBOL','trans_logFC')
# merged_data <- merge(proteomic_genes, transcriptomics, by = "SYMBOL")
# 
# # Assuming you have a data frame named merged_data
# 
# # Create a new column to store the direction match information
# merged_data$direction_match <- NA
# 
# # Check if the directions of fold changes match and assign 0, 1, or 2
# for (i in 1:nrow(merged_data)) {
#   protein_fc <- merged_data$prot_logFC[i]
#   transcript_fc <- merged_data$trans_logFC[i]
#   
#   if (protein_fc > 0 && transcript_fc > 0) {
#     direction_match <- 1  # Both upregulated
#   } else if (protein_fc < 0 && transcript_fc < 0) {
#     direction_match <- 2  # Both downregulated
#   } else {
#     direction_match <- 0  # Not concordantly regulated
#   }
#   
#   merged_data$direction_match[i] <- direction_match
# }
# 
# # Print the updated data frame with the direction match information
# Both_downregulated=sum(merged_data$direction_match==2)
# Both_downregulated
# Both_upregulated=sum(merged_data$direction_match==1)
# Both_upregulated
# Not_concordantly=sum(merged_data$direction_match==0)
# Not_concordantly
# Both_downregulated+Both_upregulated+Not_concordantly



a <- list('Proteins' = proteomic_proteins_HS_T0$SYMBOL,
          'Genes' = transcriptomics_HS_T0$SYMBOL)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Shared_genes_overall.jpg");
dev.off ()
a <- list('Upregulated: proteins' = up_proteomic_proteins,
          'Upregulated: genes' = up_transcriptomics)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Shared_upregulated.jpg");
dev.off ()
a <- list('Downregulated: proteins' = down_proteomic_proteins,
          'Downregulated: genes' = down_transcriptomics)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Shared_downregulated.jpg");
dev.off ()
shared_upregulated_p=intersect(up_transcriptomics,up_proteomic_proteins)
universe_size=length(transcriptomics_HS_T0$SYMBOL)+length(proteomic_proteins_HS_T0$SYMBOL)#length(up_transcriptomics)+length(up_proteomic_proteins)
overlap_size=length(intersect(proteomic_proteins_HS_T0$SYMBOL,transcriptomics_HS_T0$SYMBOL))
#overlap_size=1000  
result <- phyper(q = overlap_size - 1,
                 m = length(transcriptomics_HS_T0$SYMBOL),
                 n = universe_size  - length(transcriptomics_HS_T0$SYMBOL),
                 k = length(proteomic_proteins_HS_T0$SYMBOL),
                 lower.tail = FALSE)
result
result <- phyper(q = overlap_size - 1,
                 m = length(transcriptomics_HS_T0$SYMBOL),
                 n = universe_size  - length(transcriptomics_HS_T0$SYMBOL),
                 k = length(proteomic_proteins_HS_T0$SYMBOL),
                 lower.tail = TRUE)
result
######
result <- phyper(q = overlap_size - 1,
                 m = length(proteomic_proteins_HS_T0$SYMBOL),
                 n = universe_size  - length(proteomic_proteins_HS_T0$SYMBOL),
                 k = length(transcriptomics_HS_T0$SYMBOL),
                 lower.tail = FALSE)
result
result <- phyper(q = overlap_size - 1,
                 m = length(proteomic_proteins_HS_T0$SYMBOL),
                 n = universe_size  - length(proteomic_proteins_HS_T0$SYMBOL),
                 k = length(transcriptomics_HS_T0$SYMBOL),
                 lower.tail = TRUE)
result

contingency_table <- matrix(c(
  length(intersect(transcriptomics_HS_T0, proteomic_proteins_HS_T0)),  # Genes in both sets
  length(setdiff(transcriptomics_HS_T0, proteomic_proteins_HS_T0)),   # Genes in transcriptomics only
  length(setdiff(proteomic_proteins_HS_T0, transcriptomics_HS_T0)),   # Genes in proteomics only
  length(union(transcriptomics_HS_T0, proteomic_proteins_HS_T0))    # Total unique genes
), ncol = 2)

# Add row and column names for clarity
rownames(contingency_table) <- c("Overlap", "Transcriptomics Only")
colnames(contingency_table) <- c("Overlap", "Proteomics Only")

# Perform Fisher's Exact Test
fisher_result <- fisher.test(contingency_table)

# Print the result
print(fisher_result)

a <- list('Upregulated: proteins' = up_proteomic_genes,
          'Upregulated: genes' = up_transcriptomics)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Shared_upregulated.jpg");
dev.off ()
a <- list('Downregulated: proteins' = down_proteomic_genes,
          'Downregulated: genes' = down_transcriptomics)
ggvenn(a)
dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Shared_downregulated.jpg");
dev.off ()

shared_upregulated_p=intersect(up_transcriptomics,up_proteomic_genes)
shared_upregulated_p=intersect(up_transcriptomics,up_proteomic_proteins)
ind_shared_upregulated_p=match(shared_upregulated_p, transcriptomics_HS_T0$SYMBOL)

GO_shared_upregulated_g <- gost(query = shared_upregulated_p, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "fdr", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

shared_downregulated_g=intersect(down_transcriptomics,down_proteomic_genes)
shared_downregulated_p=intersect(down_transcriptomics,down_proteomic_proteins)
ind_shared_downregulated_p=match(shared_downregulated_p, transcriptomics_HS_T0$SYMBOL)

GO_shared_downregulated_g <- gost(query = shared_downregulated_p, 
                        organism = "hsapiens", ordered_query = TRUE, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = FALSE, 
                        user_threshold = 0.05, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = NULL, 
                        numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

unique_transcript_up <- setdiff(up_transcriptomics,up_proteomic_proteins)
unique_proteins_up <- setdiff(up_proteomic_proteins,up_transcriptomics)
unique_transcript_down <- setdiff(down_transcriptomics,down_proteomic_proteins)
unique_proteins_down <- setdiff(down_proteomic_proteins,down_transcriptomics)
ind_unique_transcript_up=match(unique_transcript_up, transcriptomics_HS_T0$SYMBOL)
ind_unique_proteins_up=match(unique_proteins_up, proteomic_proteins_HS_T0$SYMBOL)
ind_unique_transcript_down=match(unique_transcript_down, transcriptomics_HS_T0$SYMBOL)
ind_unique_proteins_down=match(unique_proteins_down, proteomic_proteins_HS_T0$SYMBOL)

Shared_prot_trans_up <- data.frame(qpcR:::cbind.na(shared_upregulated_p,logFC_trans[ind_shared_upregulated_p],FDR_trans[ind_shared_upregulated_p])) # Bind as columns
Shared_prot_trans_down <- data.frame(qpcR:::cbind.na(shared_downregulated_p,logFC_trans[ind_shared_downregulated_p],FDR_trans[ind_shared_downregulated_p])) # Bind as columns
Unique_transcript_up <- data.frame(qpcR:::cbind.na(unique_transcript_up,logFC_trans[ind_unique_transcript_up],FDR_trans[ind_unique_transcript_up])) # Bind as columns
Unique_transcript_down <- data.frame(qpcR:::cbind.na(unique_transcript_down,logFC_trans[ind_unique_transcript_down],FDR_trans[ind_unique_transcript_down])) # Bind as columns
Unique_proteins_up <- data.frame(qpcR:::cbind.na(unique_proteins_up,logFC_prot[ind_unique_proteins_up],FDR_prot[ind_unique_proteins_up]))
Unique_proteins_down <- data.frame(qpcR:::cbind.na(unique_proteins_down,logFC_prot[ind_unique_proteins_down],FDR_prot[ind_unique_proteins_down]))

GO_unique_upregulated_trans <- gost(query = unique_transcript_up, 
                                      organism = "hsapiens", ordered_query = TRUE, 
                                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                                      user_threshold = 0.05, correction_method = "fdr", 
                                      domain_scope = "annotated", custom_bg = NULL, 
                                      numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_downregulated_trans <- gost(query = unique_transcript_down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "fdr", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_upregulated_prot <- gost(query = unique_proteins_up, 
                                    organism = "hsapiens", ordered_query = TRUE, 
                                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                                    user_threshold = 0.05, correction_method = "fdr", 
                                    domain_scope = "annotated", custom_bg = NULL, 
                                    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)
GO_unique_downregulated_prot <- gost(query = unique_proteins_down, 
                                      organism = "hsapiens", ordered_query = TRUE, 
                                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                      measure_underrepresentation = FALSE, evcodes = FALSE, 
                                      user_threshold = 0.05, correction_method = "fdr", 
                                      domain_scope = "annotated", custom_bg = NULL, 
                                      numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = FALSE)

output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Up_down_transcripts_proteins_GO_Genes.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Shared_trans_prot_up_GO')
writeData(wb, sheet = 'Shared_trans_prot_up_GO', x = cbind(GO_shared_upregulated_g$result$term_id,GO_shared_upregulated_g$result$term_name))
addWorksheet(wb, sheetName = 'Shared_trans_prot_down_GO')
writeData(wb, sheet = 'Shared_trans_prot_down_GO', x = cbind(GO_shared_downregulated_g$result$term_id,GO_shared_downregulated_g$result$term_name))
addWorksheet(wb, sheetName = 'Unique_transcripts_g_up_GO')
writeData(wb, sheet = 'Unique_transcripts_g_up_GO', x = cbind(GO_unique_upregulated_trans$result$term_id,GO_unique_upregulated_trans$result$term_name))
addWorksheet(wb, sheetName = 'Unique_transcripts_g_down_GO')
writeData(wb, sheet = 'Unique_transcripts_g_down_GO', x = cbind(GO_unique_downregulated_trans$result$term_id,GO_unique_downregulated_trans$result$term_name))
addWorksheet(wb, sheetName = 'Unique_proteins_g_up_GO')
writeData(wb, sheet = 'Unique_proteins_g_up_GO', x = cbind(GO_unique_upregulated_prot$result$term_id,GO_unique_upregulated_prot$result$term_name))
addWorksheet(wb, sheetName = 'Unique_proteins_g_down_GO')
writeData(wb, sheet = 'Unique_proteins_g_down_GO', x = cbind(GO_unique_downregulated_prot$result$term_id,GO_unique_downregulated_prot$result$term_name))


addWorksheet(wb, sheetName = 'Shared_prot_trans_up')
writeData(wb, sheet = 'Shared_prot_trans_up', x = Shared_prot_trans_up)
addWorksheet(wb, sheetName = 'Shared_prot_trans_down')
writeData(wb, sheet = 'Shared_prot_trans_down', x = Shared_prot_trans_down)

addWorksheet(wb, sheetName = 'Unique_transcript_up')
writeData(wb, sheet = 'Unique_transcript_up', x = Unique_transcript_up)
addWorksheet(wb, sheetName = 'Unique_transcript_down')
writeData(wb, sheet = 'Unique_transcript_down', x = Unique_transcript_down)
addWorksheet(wb, sheetName = 'Unique_proteins_up')
writeData(wb, sheet = 'Unique_proteins_up', x = Unique_proteins_up)
addWorksheet(wb, sheetName = 'Unique_proteins_down')
writeData(wb, sheet = 'Unique_proteins_down', x = Unique_proteins_down)

saveWorkbook(wb, output_file)

# G=data.frame(cbind(palmieri_final@featureData@data[["SYMBOL"]],palmieri_final@assayData[["exprs"]]))
# aggregate_strings <- function(x) {
#   if (is.character(x)) {
#     if (length(x) > 1) {
#       return(x[1])
#     } else {
#       return(x)
#     }
#   } else {
#     return(mean(x))
#   }
# }
# G_mean <- aggregate(. ~ V1, data = G, FUN = aggregate_strings)
# rownames(G_mean)=G_mean$V1
# ind_down=match(shared_downregulated_p, G_mean$V1)
# ind_up=match(shared_upregulated_p, G_mean$V1)
# all_ind=c(ind_down,ind_up)
# G_mean <- subset(G_mean, select = -V1)
# exp_val=G_mean[all_ind,]
# 
# ind_down=match(shared_downregulated_p, palmieri_final@featureData@data[["SYMBOL"]])
# ind_up=match(shared_upregulated_p, palmieri_final@featureData@data[["SYMBOL"]])
# 
# Prueba=data.frame(palmieri_final@featureData@data[["SYMBOL"]][ind_down])
# symbols=rbind(shared_downregulated_p,shared_upregulated_p)
# exp_val=cbind(palmieri_final@assayData[["exprs"]][ind_down],palmieri_final@assayData[["exprs"]][ind_up])


library(dendextend)
library(Biobase)
library(multiClust)
library(preprocessCore)
library(ctc)
library(gplots)
library(dendextend)
library(graphics)
library(grDevices)
library(amap)
library(survival)
library(ggpubr)
library(factoextra)

exp_file='/Users/mariacastillo/Desktop/HEATSTROKE/DATA CEL FILES/ALL/normalized_exprs.txt'
SDRF <- read_excel("/Users/mariacastillo/Desktop/HEATSTROKE/DATA CEL FILES/ALL/metadata.xlsx")
# Load the gene expression matrix 
data.exprs <- input_file(input='/Users/mariacastillo/Desktop/HEATSTROKE/DATA CEL FILES/ALL/normalized_exprs.txt')
exp_val=read.table(exp_file, header = TRUE, sep = "", dec = ".")
colnames(exp_val)=SDRF$`Sample number`

# ranked.exprs <- probe_ranking(input=exp_file,
#                               probe_number=600, 
#                               probe_num_selection="Fixed_Probe_Num",
#                               data.exp=data.exprs, 
#                               method="SD_Rank")
# 
# hclust_analysis <- cluster_analysis(sel.exp=exp_val,
#                                     cluster_type="HClust",
#                                     distance="euclidean", linkage_type="ward.D2", 
#                                     gene_distance="correlation",
#                                     num_clusters=3, data_name="Heat Stroke", 
#                                     probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
#                                     cluster_num_selection="Fixed_Clust_Num")
# kmeans_analysis <- cluster_analysis(sel.exp=exp_val,
#                                     cluster_type="Kmeans",
#                                     distance=NULL, linkage_type=NULL, 
#                                     gene_distance=NULL, num_clusters=3,
#                                     data_name="Heat Stroke", probe_rank="SD_Rank",
#                                     probe_num_selection="Fixed_Probe_Num",
#                                     cluster_num_selection="Fixed_Clust_Num")

rownames(Shared_prot_trans_up)=Shared_prot_trans_up$shared_upregulated_p
d1_up <- dist(Shared_prot_trans_up,method = "euclidean")
d3_up <- hclust(d1_up, method = "average")
dend_up <- as.dendrogram(d3_up)
labels_up=cutree(d3_up, k=5)
dend_up %>% labels
labels_up_f=data.frame(labels_up)
labels_up_f$genes=rownames(labels_up_f)
#labels(dend_up) <- Shared_prot_trans_up$shared_upregulated_p[dend_up %>% labels]
dend_up %>% plot
dend_up %>% set("branches_k_color", k = 5) %>% plot(main = "Dendogram shared upregulated genes")

dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Dendogram_upregulated.jpg");
dev.off ()

rownames(Shared_prot_trans_down)=Shared_prot_trans_down$shared_downregulated_p
d1_down <- dist(Shared_prot_trans_down,method = "euclidean")
d3_down <- hclust(d1_down, method = "average")
dend_down <- as.dendrogram(d3_down)
labels_down=cutree(d3_down, k=5)
dend_down %>% labels
labels_down_f=data.frame(labels_down)
labels_down_f$genes=rownames(labels_down_f)
dend_down %>% plot
dend_down %>% set("branches_k_color", k = 5) %>% plot(main = "Dendogram shared downregulated genes")

dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Dendogram_downregulated.jpg");
dev.off ()

# Load the required library
library(dplyr)

# Define the file path to your CSV file
csv_file <- "/Users/mariacastillo/Desktop/HEATSTROKE/Shared genes/upregulated_GO.csv"
# Load the CSV file into a data frame
data_up <- read.csv(csv_file, stringsAsFactors = FALSE)
# Define a list of genes you want to search for
genes_to_search <- Shared_prot_trans_up$shared_upregulated_p
# Function to find related term_ids for a list of genes
find_related_terms <- function(gene_list,data) {
  related_terms <- data %>%
    filter(intersection_size > 0) %>%
    filter(strsplit(intersections, ",") %>% sapply(function(x) any(x %in% gene_list)))
  return(related_terms$term_id)
}

finale_df_up <- data.frame(Gene = character(0), Term_id = numeric(0))

for (gene1 in genes_to_search) {
# Find related term_ids for the list of genes
related_term_ids <- find_related_terms(gene1,data_up)
# Print the related term_ids
related_term_chr=paste(related_term_ids, collapse=", ")
length(related_term_chr)
new_row <- data.frame(Gene = gene1, Term_id = related_term_chr)
# Append the new row to the data frame
finale_df_up <- rbind(finale_df_up, new_row)
}
# Define the output Excel file path
output_excel_file <- "/Users/mariacastillo/Desktop/HEATSTROKE/Shared genes/gene_functions_up.xlsx"
write.xlsx(finale_df_up, output_excel_file)


# Define the file path to your CSV file
csv_file <- "/Users/mariacastillo/Desktop/HEATSTROKE/Shared genes/downregulated_GO.csv"
# Load the CSV file into a data frame
data_down <- read.csv(csv_file, stringsAsFactors = FALSE)
finale_df_down <- data.frame(Gene = character(0), Term_id = numeric(0))
genes_to_search <- Shared_prot_trans_down$shared_downregulated_p

for (gene1 in genes_to_search) {
  # Find related term_ids for the list of genes
  related_term_ids <- find_related_terms(gene1,data_down)
  # Print the related term_ids
  related_term_chr=paste(related_term_ids, collapse=", ")
  new_row <- data.frame(Gene = gene1, Term_id = related_term_chr)
  # Append the new row to the data frame
  finale_df_down <- rbind(finale_df_down, new_row)
}

# Define the output Excel file path
output_excel_file <- "/Users/mariacastillo/Desktop/HEATSTROKE/Shared genes/gene_functions_down.xlsx"
write.xlsx(finale_df_down, output_excel_file)

# Load required libraries
library(dendextend)

All_functions_up=unique(unlist(strsplit(finale_df_up$Term_id, ",")))
# Create a binary matrix
binary_matrix_up <- matrix(0, nrow = nrow(finale_df_up), ncol = length(All_functions_up))
colnames(binary_matrix_up) <- All_functions_up
rownames(binary_matrix_up) <- finale_df_up$Gene

# Fill in the binary matrix
for (i in 1:nrow(finale_df_up)) {
  functions_gene_i=unique(unlist(strsplit(finale_df_up$Term_id[i], ",")))
  for (j in functions_gene_i){
    binary_matrix_up[finale_df_up$Gene[i], j] <- 1
  }
}

d1_up <- dist(binary_matrix_up,method = "binary")
d3_up <- hclust(d1_up, method = "average")
dend_up <- as.dendrogram(d3_up)
labels_up_GO=cutree(d3_up, k=5)
labels_up_GO_f=data.frame(labels_up_GO)
labels_up_GO_f$genes=rownames(labels_up_GO_f)
dend_up %>% labels
dend_up %>% plot
dend_up %>% set("branches_k_color", k = 3) %>% plot(main = "Dendogram shared upregulated genes")

dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Dendogram_upregulated_GO.jpg");
dev.off ()

All_functions_down=unique(unlist(strsplit(finale_df_down$Term_id, ",")))
# Create a binary matrix
binary_matrix_down <- matrix(0, nrow = nrow(finale_df_down), ncol = length(All_functions_down))
colnames(binary_matrix_down) <- All_functions_down
rownames(binary_matrix_down) <- finale_df_down$Gene

# Fill in the binary matrix
for (i in 1:nrow(finale_df_down)) {
  functions_gene_i=unique(unlist(strsplit(finale_df_down$Term_id[i], ",")))
  for (j in functions_gene_i){
    binary_matrix_down[finale_df_down$Gene[i], j] <- 1
  }
}

d1_down <- dist(binary_matrix_down,method = "binary")
d3_down <- hclust(d1_down, method = "average")
dend_down <- as.dendrogram(d3_down)
labels_down_GO=cutree(d3_down, k=5)
labels_down_GO_f=data.frame(labels_down_GO)
labels_down_GO_f$genes=rownames(labels_down_GO_f)
dend_down %>% labels
dend_down %>% plot
dend_down %>% set("branches_k_color", k = 3) %>% plot(main = "Dendogram shared downregulated genes")

dev.copy(jpeg,filename="/Users/mariacastillo/Desktop/HEATSTROKE/Transcriptomics/Dendogram_downregulated_GO.jpg");
dev.off ()

output_file <- '/Users/mariacastillo/Desktop/HEATSTROKE/Shared genes/Dendogram_levels.xlsx'
wb <- createWorkbook()
addWorksheet(wb, sheetName = 'Upregulated expression')
writeData(wb, sheet = 'Upregulated expression', x = labels_up_f)
addWorksheet(wb, sheetName = 'Downregulated expression')
writeData(wb, sheet = 'Downregulated expression', x = labels_down_f)
addWorksheet(wb, sheetName = 'Upregulated GO')
writeData(wb, sheet = 'Upregulated GO', x = labels_up_GO_f)
addWorksheet(wb, sheetName = 'Downregulated GO')
writeData(wb, sheet = 'Downregulated GO', x = labels_down_GO_f)

saveWorkbook(wb, output_file)








# library(umap)
# library(dplyr)      # For data manipulation
# library(ggplot2)    # For visualization
# 
# umap_result <- umap(exp_val, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
# # Create a ggplot visualization
# ggplot(umap_result$layout, aes(x = V1, y = V2)) +
#   geom_point() +
#   ggtitle("UMAP Visualization of Gene Expression Data")










# Loading package
library(ClusterR)
library(cluster)

n_clust<-fviz_nbclust(data.frame(Shared_prot_trans_up[,2]), kmeans, method = "silhouette",k.max = 30)
kmeans_up =kmeans(data.frame(Shared_prot_trans_up[,2]), centers=5, iter.max = 10, nstart = 1)
plot(Shared_prot_trans_up[, 2], col = kmeans_up$cluster)

n_clust<-fviz_nbclust(data.frame(Shared_prot_trans_down[,2]), kmeans, method = "silhouette",k.max = 30)
kmeans_down =kmeans(data.frame(Shared_prot_trans_down[,2]), centers=5, iter.max = 10, nstart = 1)
plot(Shared_prot_trans_down[, 2], col = kmeans_down$cluster)

