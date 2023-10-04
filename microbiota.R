setwd("/Users/zhoucheng/Desktop/肠道孟德尔/菌群文件汇总")
#* ---install or library-------------------------------------
install.packages("R.utils")
install.packages("data.table")
install.packages("devtools")
install.packages("bidirectional MR")
install.packages("MRDirection")
install.packages("MRPRESSO")
install.packages("ivpack")


if (!require("devtools")) { install.packages("devtools") } else {}
if (!require("data.table")) { install.packages("data.table") } else {}
if (!require("TwoSampleMR")) { devtools::install_github("MRCIEU/TwoSampleMR") } else {}
if (!require("MRInstruments")) { devtools::install_github("MRCIEU/MRInstruments") } else {}
library(MRInstruments)
library(plyr)
library(dplyr)  
library(data.table)
library(R.utils)
library(devtools)
library(data.table)
library(TwoSampleMR)
library(stringr)
library(readxl)
library(MRPRESSO)
library(MendelianRandomization)


FileNames<-list.files(paste0(getwd()),pattern=".txt.gz")

# ---本地菌群数据处理-------------------------------------需要大的数据包
`%+%` <- function(x,y) paste0(x,y)
dir.create(getwd()%+%"/data")
path_data <- getwd()%+%"/data"
dir.create(getwd()%+%"/result")
path_res <- getwd()%+%"/result"

for (i in 1:length(FileNames)) {
  d1 <- try(fread(paste0(getwd(), "/", FileNames[i]), sep = "\t"), silent = TRUE)
  if (inherits(d1, "data.frame")) {  # 检查是否返回的是数据框对象
    d2 <- subset(d1, P.weightedSumZ < 1e-5)  # 使用简化形式
    d3 <- d2[, c(1, 4, 5, 6, 7, 8, 10, 11)]
    write.table(d3, paste0(path_data, "/", FileNames[i]), sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
  } else {
  }
}
# ---本地菌群数据处理-------------------------------------需要大的数据包



########规范化数据并进行连锁不平衡分析##################
########规范化数据并进行连锁不平衡分析##################
exp_dat <- list()
ex_pore <- c()
for(i in c(1:length(FileNames))){
  IV <- fread(paste0(path_data,"/",FileNames[i]))
  IV$PHENO <- FileNames[i]
  IV$eaf <- ifelse(IV$ref.allele == IV$eff.allele, IV$beta, 1 - IV$beta)
  IV1<-format_data(IV,
                   type="exposure",
                   phenotype_col = "PHENO",
                   snp_col = "rsID",
                   beta_col = "beta",
                   se_col = "SE",
                   pval_col = "P.weightedSumZ",
                   samplesize_col = "N",
                   effect_allele_col = "eff.allele",
                   other_allele_col = "ref.allele",
                   eaf_col = "eaf") 
  
  IV1 <- clump_data(IV1,clump_kb = 500,clump_r2 = 0.1)
  exp_dat[[i]] <- IV1
  ex_pore<-c(ex_pore,FileNames[i])
}


##############集合写出菌群所有文件###########————————————————————————————
exp_dat_combined <- do.call(rbind, exp_dat)
write.csv(exp_dat_combined, "exp_dat.csv", row.names = FALSE)




#save.image("ex_dat.Rdata")
#ex_pore <- c()
#for(i in c(1:length(FileNames))){
#  ex_pore<-c(ex_pore,FileNames[i])
#}
#save.image("ex_pore.Rdata")
#allSNP  <- do.call(rbind, exp_dat)
#View(allSNP)
#暴露做完可以保存


#############选择结局的数据#######################
#############选择结局的数据#######################
##############本地结局111111################
dataname1="/Users/zhoucheng/Desktop/肠道孟德尔/菌群文件汇总/GCST90018804_buildGRCh37.tsv"
GWAS_1 <- fread(dataname1)
allSNP <- do.call(rbind, exp_dat)
GWAS_2 <- subset(GWAS_1, GWAS_1$variant_id %in% allSNP$SNP & !is.na(GWAS_1$beta))
rm(GWAS_1)
GWAS_2$PHENO<-"HBV"
out_data <- format_data(GWAS_2,
                        type="outcome",
                        phenotype_col = "PHENO",
                        snp_col = "variant_id",
                        beta_col = "beta",
                        se_col = "standard_error",
                        eaf_col = "effect_allele_frequency",
                        pval_col = "p_value",
                        effect_allele_col = "effect_allele",
                        other_allele_col = "other_allele",)

out_dat <- list()
out_dat[[1]] <- out_data


##############本地结局2222222222################
dataname1="/Users/zhoucheng/Desktop/肠道孟德尔/菌群文件及代码/summary_stats_finngen_R9_CIRRHOSIS_BROAD.gz"
GWAS_1 <- fread(dataname1)
allSNP <- do.call(rbind, exp_dat)
GWAS_2 <- subset(GWAS_1, GWAS_1$rsids %in% allSNP$SNP & !is.na(GWAS_1$beta))
rm(GWAS_1)
GWAS_2$PHENO<-"HBV"
out_data <- format_data(GWAS_2,
                        type="outcome",
                        phenotype_col = "PHENO",
                        snp_col = "rsids",
                        beta_col = "beta",
                        se_col = "sebeta",
                        eaf_col = "af_alt",
                        pval_col = "pval",
                        samplesize_col = "N",
                        effect_allele_col = "ref",
                        other_allele_col = "alt",)


##############在线的结局########################
out_data<-extract_outcome_data(allSNP$SNP,"bbj-a-105")
length(unique(allSNP$SNP))
View(out_data)
out_dat <- list()
out_dat[[1]] <- out_data



################开始跑组合循环代码################
################开始跑组合循环代码################
out_come<-c("HBV")
results <- list()
for (i in c(1:length(ex_pore))){
  for (j in c(1:length(out_come))){
    dat <- harmonise_data(
      exposure_dat = exp_dat[[i]],
      outcome_dat = out_dat[[j]],
      action = 2)
    if (length(dat$SNP) == 0) next
    res <- mr(dat)
    res$exposure=ex_pore[i]
    res$outcome=out_come[j]
    print(paste0("------", ex_pore[i], " & ",out_come[j],"------"))
    
    print(generate_odds_ratios(res))
    
    #primary results
    results[[length(out_come)*(i-1)+j]]    <- generate_odds_ratios(res)
  }
}

####################规范化结局并导出######################
####################规范化结局并导出######################
results_allIV  <- do.call(rbind, results)
#format(round(x, 2), nsmall = 2)
results_allIV$estimate <- paste0(format(round(results_allIV$or, 2), nsmall = 2), " (", 
                                 format(round(results_allIV$or_lci95, 2), nsmall = 2), "-",
                                 format(round(results_allIV$or_uci95, 2), nsmall = 2), ")")

row_x <- rownames (results_allIV[which(results_allIV$pval > 0.05), ])

results_allIV$pvalue          <- format(results_allIV$pval, scientific = TRUE, digits = 2)
results_allIV[row_x, ]$pvalue <- format(round(results_allIV[row_x, ]$pval, 2), nsmall = 2)

outname1="/results.csv"
write.table(results_allIV[,c(3:ncol(results_allIV))],path_res%+%outname1,sep = ",",row.names =FALSE,col.names =TRUE,quote =TRUE)








####其他分析####-------------------------------------
#####此部分为单一数据-----------------
for (i in c(1:length(ex_pore))){
  for (j in c(1:length(out_come))){
    dat <- harmonise_data(
      exposure_dat = exp_dat[[66]],
      outcome_dat = out_dat[[j]],
      action = 2)
  }
}

#####此部分出结果了-----------------
######Cochran's Q p值
het <- mr_heterogeneity(dat)
######MR-Egger intercept p值
pleio <- mr_pleiotropy_test(dat)
######MR-Presso p值
mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
          SignifThreshold = 0.05)

######离群图
single <-mr_leaveoneout(dat)
mr_leaveoneout_plot(single)
######散点图
res <- mr(dat)
mr_scatter_plot(res, dat)
#######漏斗图
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
mr_funnel_plot(res_single)
#######计算eaf.exposure的平均值，将缺失值替换为平均值，计算F值
mean_eaf <- mean(dat$eaf.exposure, na.rm = TRUE)
dat$eaf.exposure[is.na(dat$eaf.exposure)] <- mean_eaf
write.csv(dat, "dateaf.csv", row.names = FALSE)

###########进行PhenoScannerSNP混淆因素检测####--------------
confounders <- phenoscanner(dat$SNP)
write.csv(confounders[[1]], "confounders_with.csv", row.names = FALSE)
write.csv(confounders[[2]], "confounders_without.csv", row.names = FALSE)



#############反向孟德尔随机化，先计算出相关性系数
dat$r_exposure <- dat$beta.exposure / (dat$se.exposure * sqrt(dat$samplesize.exposure))
dat$r_outcome <- dat$beta.outcome / (dat$se.outcome * sqrt(dat$samplesize.outcome))
fx1 <- mr_steiger2(r_exp=dat$r_exposure,
                   r_out=dat$r_outcome,
                   n_exp=dat$samplesize.exposure,
                   n_out=dat$samplesize.outcome,
                   r_xxo=1,
                   r_yyo=1)



