###This function provides comprehensive longitudinal statistics on subsets of a data

## Inputs:

## data=data set
## sub_vars=multiple columns used for stratification 
## tm_var=field indicating the time variable
## dep_vars=multiple columns of dependent variables for which statistics will be performed
## rand_var=single column specifying the random effect variable, which is usually a subject identifier 
## normality=indication of whether the dependent variable is normally distributed (in which 
## case a parametric test will be performed) or not normality distributed (in which case a non
## parametric test will be performed)
## adj= The method in which the p values are adjusted 
## pw=pairwise statistic input ("same", "separate")
##     "same" indicates the same longitudinal test (parametric linear mixed effects models or nonparametric Friedman test) will also be applied to pairwise comparisons
##     "separate" indicates a separate pairwise test (parametric paired t-test or nonparametric paired Wilcoxon Signed rank test)  will be applied to pairwise comparisons
## Outputs:

## A table that includes each subset of data with pairwise comparisons between groups with the following fields
## Measure=the dependent variable (determined by dep_vars)
## Additional subset variables (depends on sub_vars)
## group1=one of the time points (determined by tm_var)
## group2=another one of the time points (determined by tm_var)
## p-value-overall=p value obtained from linear mixed effects models (normally distributed data) or Friedman tests (not normally distributed data) relating each dependent variable to time
## p-value = pairwise p value from post HOC multiple comparison test identifying statistical difference between 2 time points
## padj = adjusted p value from the pairwise comparison between 2 time points
## p.signif=Asterisks indicating significance

long_stat_sum <- function(data, sub_vars, tm_var, dep_vars, rand_var, normality, adj="bonferroni", pw="same") {
  
  #First we pull the length of the subset variables
  sub_var_num <- length(sub_vars)
  #We generate a blank list that we subsequently add to
  sub_var_list <- list()
  #We go through each subset variables and add to create a list of variables with select values
  for (j in 1:sub_var_num) {
    #We pull the individual subset variables
    sub_name_j <- sub_vars[j]
    #We then pull their values
    sub_vals_j <- unlist(as.vector(unique(data[,which(colnames(data)==sub_name_j)])))
    #We then add them to the list
    sub_var_list[[sub_name_j]] <- sub_vals_j
  }
  #Now we need to pull the timepoint values to generate pairwise comparisons
  tm_vals <- as.vector(unique(data[,which(colnames(data)==tm_var)]))
  #We generate pairwise combinations of 2
  tm_combs <- combn(unlist(tm_vals), 2, simplify = FALSE)
  #We collapse them into one string 
  tm_strings <- sapply(tm_combs, function(x) paste(x, collapse = "-xxx-"))
  #We now add them to the sub_var_list to allow us to generate a table
  sub_var_list[["combo"]] <- tm_strings
  #Now we generate a table with all the subset variables with corresponding values as well as a combo field with the pairwise combination
  subset_table <- expand.grid(sub_var_list)
  #Now we broaden to include that with every dependent variable 
  #First we pull the number of dep_vars
  dep_var_num <- length(dep_vars)
  #Now we generate a blank table that we add to with each dependent variable
  stats_sum_table <- c()
  #We go through each dependent variable to add it to the stats_sum_table
  for (d in 1:dep_var_num) {
    #For each variable we adjoin a column in the front with the dependent variable name
    stats_sum_table_d <- cbind(rep(dep_vars[d],nrow(subset_table)),subset_table)
    #We rename that first adjoining column Measure
    colnames(stats_sum_table_d)[1] <- "Measure"
    #We now adjoin it to the stats_sum_table
    stats_sum_table <- rbind(stats_sum_table,stats_sum_table_d)
  }
  #We now add new fields for additional specified fields
  #First we split the combination into 2 groups
  stats_sum_table <- separate(stats_sum_table, combo, into = c("group1", "group2"), sep = "-xxx-")
  #Now we add blank fields for the stats
  #One for the overall p-value 
  stats_sum_table["p-value-overall"] <- NA
  #One for the pairwise multiple comparison p-value
  stats_sum_table["p-value"] <- NA
  #One for the pairwise multiple comparison adjusted p-value
  stats_sum_table["padj"] <- NA
  #One to indicate the significance
  stats_sum_table["p.signif"] <- ""
  
  #Before the statistics are performed we need to extract the columns for the time and random variable
  time_col <- which(colnames(data)==tm_var)
  rand_col <- which(colnames(data)==rand_var)
  #Now we go through th stats summary table and subset the data and then perform the specified analysis 
  for (r in 1:nrow(stats_sum_table)) {
    #First we pull the column for the Measure
    #Below includes an edit around NA coluumn name values
    #If the column name is not NA we just pull which column is names the same as the measure
    if (!is.na(stats_sum_table$Measure[r])) {
      meas_col <- which(colnames(data)==stats_sum_table$Measure[r])  
    #If the column name is NA we use the na function
    } else {
      meas_col <- which(is.na(colnames(data)))  
    }
    #Now we subset the data by repeatedly filtering based on the subset filters 
    #We start with the full data set filtered by the time point
    data_r <- data
    #First we filter by the timepoints to only include the specified groups
    #First we change the time column name for simplicity
    colnames(data_r)[time_col] <- "Time"
    #We then go through the columns and subsequently filter out to only include rows with the specified the values
    for (c in 2:(ncol(stats_sum_table)-6)) {
      #We pull the column name
      colname_c <- as.character(colnames(stats_sum_table)[c])
      #We then identify that column in the data set
      col_c <- which(colnames(data) == colname_c)
      colnames(data_r)[col_c] <- "Filter"
      #We pull the column value
      sub_val_c <- as.character(stats_sum_table[r,c])
      #We now iteratively filter the data by that each column being each specified value
      data_r <- data_r[which(data_r$Filter==sub_val_c),]
      colnames(data_r)[col_c] <- colname_c
    }
    #Now we focus on the columns of the measures and the groups for the statistics ignoring the remaining columns
    data_r_filt_pre <- data_r[,c(meas_col,time_col,rand_col)]
    #We will also be sure to filter so that each random variable appears at least more than once to avoid singularities 
    #To do that, we need to specify the column names
    colnames(data_r_filt_pre) <- c("Amount", "Time", "ID")
    #Before the stats we want to make sure our Amount field is structured numerically 
    data_r_filt_pre$Amount <- as.numeric(data_r_filt_pre$Amount)
    #Now we specify those that appear more than once and have no NA values 
    data_r_filt <- data_r_filt_pre %>%
      filter(!is.na(Amount)) %>%
      group_by(ID) %>%
      filter(n() > 1) %>%
      ungroup() 
    #To avoid singularities we ensure there are at least two random variables present to continue with the statistics
    if (nrow(data_r_filt)>4) {
      #If the dependent variable data is normally distributed we use a linear mixed effect model      
      if (normality) {
        ##Overall comparison##
        
        #First we do a linear mixed effect with the data relating the dependent variable to time 
        model_r <- lmer(Amount ~ Time + (1 | ID), data = data_r_filt)
        #The mixed function is used to estimate mixed models with lme4 and calculate p-values for time fixed effects
        mixed_r <- mixed(model_r, data=data_r_filt)
        #We then use the summary function to generate a coefficient matrix 
        sum_r <- summary(mixed_r)
        #We then pull the p-value and put it in the stat summary table
        stats_sum_table$`p-value-overall`[r] <- sum_r$coefficients[2,5]
        
        ###Pairwise comparison##
        
        #Now we do additional tests for pairwise comparisons
        #First we need to pull the two groups for the comparisons
        time_1r <- stats_sum_table$group1[r]
        time_2r <- stats_sum_table$group2[r]
        #We then generate a pairwise dataset for comparisons between the two groups and do all of those stats as well 
        data_r_filt_pw <- data_r_filt %>% 
                          filter(Time %in% c(time_1r,time_2r))
        #Now if the pairwise input is "same" we do the same test on the pairwise subset
        if (pw=="same") {
          #Now we do a linear mixed effect model with the data relating the dependent variable to time with a pairwise comparison
          model_r_pw <- lmer(Amount ~ Time + (1 | ID), data = data_r_filt_pw)
          #The mixed function is used to estimate mixed models with lme4 and calculate p-values for time fixed effects
          mixed_r_pw <- mixed(model_r_pw, data=data_r_filt_pw)
          #We then use the summary function to generate a coefficient matrix 
          sum_r_pw <- summary(mixed_r_pw)
          #We now pull the pairwise p-value and put it in the stat summary table
          stats_sum_table$`p-value`[r] <- sum_r_pw$coefficients[2,5]
        } 
        #On the other hand if pairwise input is "different" we do a separate test (for normal data this will be a paired t test) 
        if (pw=="separate") {
          #The data needs to be restructed to be in the wide format
          data_r_filt_pw_wide <- pivot_wider(data_r_filt_pw,
                                             id_cols="ID",
                                             names_from="Time",
                                             values_from="Amount")
          #For simplicity in the stats the output needs to be in a data frame
          data_r_filt_pw_wide <- as.data.frame(data_r_filt_pw_wide)
          #We perform a t-test
          t_test_r <- t.test(as.numeric(data_r_filt_pw_wide[,2]),as.numeric(data_r_filt_pw_wide[,3]), paired=TRUE)
          stats_sum_table$`p-value`[r] <- t_test_r$p.value
        }

      } 
      #If the dependent variable data is not normally distributed we use a Friedman Test      
      if (!normality) {
        ##Overall comparison##
        
        #We only include values that have all three timepoints

        overall_data_r_filt <- data_r_filt %>%
                                group_by(ID) %>%
                                filter(n() == length(unlist(tm_vals))) %>%
                                ungroup()
        
        stat_test_r <- friedman.test(Amount ~ Time | ID, data = overall_data_r_filt)
        
        #From that we pull out the  p-value
        stats_sum_table$`p-value-overall`[r] <- stat_test_r$p.value
        
        ###Pairwise comparison##
        
        #Now we do additional tests for pairwise comparisons
        #First we need to pull the two groups for the comparisons
        time_1r <- stats_sum_table$group1[r]
        time_2r <- stats_sum_table$group2[r]
        #We then generate a pairwise data set for comparisons between the two groups and do all of those stats as well 
        data_r_filt_pw <- data_r_filt %>% 
          filter(Time %in% c(time_1r,time_2r))  %>%    
          group_by(ID) %>%
          filter(n() > 1) %>%
          ungroup()
        #We have to re-level the times to perform the necessary statistics
        data_r_filt_pw$Time <- factor(data_r_filt_pw$Time, levels=c(time_1r,time_2r))
        #Now if the pairwise input is "same" we do the same test on the pairwise subset
        if (pw=="same") {
          #We now perform the Friedman test on the pairwise data
          stat_test_r_pw <- friedman.test(Amount ~ Time | ID, data = data_r_filt_pw)
          #From that we pull out the  p-value
          stats_sum_table$`p-value`[r] <- stat_test_r_pw$p.value
        } 
        #On the other hand if pairwise input is "different" we do a separate test (for data not normally distributed this will be a paired Wilcoxon) 
        if (pw=="separate") {
          #The data needs to be restructured to be in the wide format
          data_r_filt_pw_wide <- pivot_wider(data_r_filt_pw,
                                             id_cols="ID",
                                             names_from="Time",
                                             values_from="Amount")
          #For simplicity in the stats the output needs to be in a data frame
          data_r_filt_pw_wide <- as.data.frame(data_r_filt_pw_wide)
          #We perform a paired wilcoxon signed rank
          wilcox_r <- wilcox.test(as.numeric(data_r_filt_pw_wide[,2]),as.numeric(data_r_filt_pw_wide[,3]), paired=TRUE)
          stats_sum_table$`p-value`[r] <- wilcox_r$p.value
        }
      }
    } else {
      #If there is not enough data the other p values do not apply
      stats_sum_table$`p-value-overall`[r] <- NA
      stats_sum_table$`p-value`[r] <- NA
      
    }
  }

  #Now for each subset we perform a p-value adjustment (this is only relevant when there are more than two timepoints)
  #P-values are adjusted groupwise meaning we only correct p values for one measure and one subset (defined by all the sub_vars)
  #To do this we first need to generate a subs table that just includes each measure and subset variables 
  subs <- unique(stats_sum_table[,1:(ncol(stats_sum_table)-6)])
  
  #Now we go through those subsets and do p-value adjustment
  for (l in 1:nrow(subs)) {
    meas_l <- subs$Measure[l]
    #Now we pull the rows from the stat summary table
    stat_tab_rows <- which(stats_sum_table$Measure==meas_l)
    #We go to pull the rows by sequentially filtering based on the other subset variables 
    for (k in 2:ncol(subs)) {
      rows_k <- which(stats_sum_table[,k]==subs[l,k])
      stat_tab_rows <- intersect(stat_tab_rows,rows_k)
    }
    #Now that we have the rows in the stats table that correspond to the measure and subset, we adjust
    stats_sum_table$padj[stat_tab_rows] <- p.adjust(stats_sum_table$`p-value`[stat_tab_rows], method=adj)
    #Now we add significance asterisks
    for (m in stat_tab_rows) {
      if (stats_sum_table$padj[m]<0.001) {
        stats_sum_table$p.signif[m] <- "***"
        # If p-value is below not under 0.001 but is under 0.01 two stars are provided
      } else if (stats_sum_table$padj[m] < 0.01) {
        stats_sum_table$p.signif[m] <- "**"
        # If p-value is below not under 0.01 but is under 0.05 one star is provided
      } else if (stats_sum_table$padj[m] < 0.05) {
        stats_sum_table$p.signif[m] <- "*"
        # If p-value is below not under 0.05, "ns" provided to indicate not significant
      } else {
        stats_sum_table$p.signif[m] <- "ns"
      }
    }
  }
  
  
#Now we return the final stats table  
return(stats_sum_table)
} 