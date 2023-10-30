library("loo")

setwd("C:/Users/nelso/OneDrive/BNelson/Consulting/Oceans initiative/Working/srkw_modeling")

# Read in model objects
mod_1<- readRDS("mod_1.rds")
mod_2<- readRDS("mod_2.rds")
mod_3<- readRDS("mod_3.rds")
mod_4<- readRDS("mod_4.rds")
mod_5<- readRDS("mod_5.rds")
mod_6<- readRDS("mod_6.rds")
mod_7<- readRDS("mod_7.rds")
mod_8<- readRDS("mod_8.rds")
mod_9<- readRDS("mod_9.rds")
mod_10<- readRDS("mod_10.rds")
mod_11<- readRDS("mod_11.rds")
mod_12<- readRDS("mod_12.rds")
mod_13<- readRDS("mod_13.rds")
mod_14<- readRDS("mod_14.rds")
mod_15<- readRDS("mod_15.rds")
mod_16<- readRDS("mod_16.rds")

mod_fits<- list(mod_1, mod_2, mod_3, mod_4, mod_5, mod_6, mod_7, mod_8, mod_9, mod_10, 
                mod_11, mod_12, mod_13, mod_14, mod_15, mod_16)
n_models<- length(mod_fits)

# Calculate LOO
loo_ic<- vector("list", n_models)

## extract log densities from JAGS objects
for(i in 1:n_models){
  ## convert mcmc.list to matrix
  tmp_lp <- as.matrix(mod_fits[[i]]$model$BUGSoutput$sims.matrix)
  ## extract pointwise likelihoods
  tmp_lp <- tmp_lp[,grepl("lp_", colnames(tmp_lp))]
  #tmp_lp <- tmp_lp[,grepl("lp_n", colnames(tmp_lp))]
  ## if numerical underflows, convert -Inf to 5% less than min(likelihood)
  if(any(is.infinite(tmp_lp))) {
    tmp_lp[is.infinite(tmp_lp)] <- NA
    tmp_min <- min(tmp_lp, na.rm = TRUE)
    tmp_lp[is.na(tmp_lp)] <- tmp_min * 1.05
  }
  ## calculate LOOIC
  loo_ic[[i]] <- loo(tmp_lp, cores = 8)
}

## LOOIC for all data
looic_table <- round(compare(x = loo_ic), 2)

rownames(looic_table) <- sub("model", "", rownames(looic_table))

looic_table <- looic_table[order(as.numeric(rownames(looic_table))), ]

looic_table <- cbind(model = rep(c("Null","SRKW","SRKW+NRKW", "SRKW+NRKW (High K)"), 4), as.data.frame(looic_table))

looic_table[order(looic_table[,"looic"]), ]

#write.csv(looic_table, "looic_table.csv", row.names = TRUE)