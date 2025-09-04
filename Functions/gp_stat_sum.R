###This function provides comprehensive group statistics on subsets of a data

## Inputs:

## data=data set
## sub_vars=multiple columns used for stratification 
## dep_vars=multiple columns of dependent variables for which statistics will be performed
## gp_var=single column specifying the group variable in which statistics are performed across 
## normality=indication of whether the dependent variable is normally distributed (in which 
## case a parametric test will be performed) or not normality distributed (in which case a non
## parametric test will be performed)
## adj= The method in which the Dunn's test adjust for p values (only applies to non-parametric analysis)

## Outputs:

## A table that includes each subset of data with pairwise comparisons between groups with the following fields
## Measure=the dependent variable (determined by dep_vars)
## Additional subset variables (depends on sub_vars)
## group1=one of the groups (determined by gp_var)
## group2=another one of the groups (determined by gp_var)
## p-value-overall=p value obtained from Kruskal Wallace or ANOVA indicating statistical difference between at least two of the groups
## p-value = pairwise p value from post HOC multiple comparison test identifying statistical difference between group 1 and group 2 (not available for the parametric Tukey's test)
## padj = adjusted p value from the pairwise comparison between group1 and group 2
## p.signif=Asterisks indicating significance


gp_stat_sum <- function(data, sub_vars, dep_vars, gp_var, normality, adj="bonferroni") {
  
  #First we pull the length of the subset variables
  sub_var_num <- length(sub_vars)
  #We generate a blank list that we subsequently add to
  sub_var_list <- list()
  #We go through each subset variables and add to create a list of variables with select values
  for (j in 1:sub_var_num) {
    #We pull the individual subset variables
    sub_name_j <- sub_vars[j]
    #We then pull their values
    sub_vals_j <- as.vector(unique(unlist(data[,which(colnames(data)==sub_name_j)])))
    #We then add them to the list
    sub_var_list[[sub_name_j]] <- sub_vals_j
  }
  #Now we need to pull the group values to generate pairwise comparisons
  gp_vals <- as.vector(unique(unlist(data[,which(colnames(data)==gp_var)])))
  #We generate pairwise combinations of 2
  gp_combs <- combn(gp_vals, 2, simplify = FALSE)
  #We collapse them into one string 
  gp_strings <- sapply(gp_combs, function(x) paste(x, collapse = "-xxx-"))
  #We now add them to the sub_var_list to allow us to generate a table
  sub_var_list[["combo"]] <- gp_strings
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
  #One to indicate thsesignificance
  stats_sum_table["p.signif"] <- ""
  #Now we go through ever combination of subsets of data and measure and perform the statistical comparisons
  #Because the 5th to last columns though the last column specify the pairwise groups and the p values we  can only focus on the measure and soubset variables in the subs dataframe
  subs <- unique(stats_sum_table[,1:(ncol(stats_sum_table)-6)])
  #Now we go thorugh each measure (first column in subs) and data subset filters (the remaining columns) and pull that 
  #data do statistics and add statistics to the stats table
  for (r in 1:nrow(subs)) {
    #We pull the measure
    meas_r <- as.character(subs$Measure[r])
    #We pull the measure columns from the data set
    meas_col <- which(colnames(data)==meas_r)
    #Now we subset the data by repeatedly filtering based on the subset filters 
    #We start with the full data set
    data_r <- data
    #We then go through the columns and subsequently filter out to only include rows with the specified the values
    for (c in 2:ncol(subs)) {
      #We pull the column name
      colname_c <- as.character(colnames(subs)[c])
      #We then identify that column in the data set
      col_c <- which(colnames(data) == colname_c)
      colnames(data_r)[col_c] <- "Filter"
      #We pull the column value
      sub_val_c <- as.character(subs[r,c])
      #We now iteratively filter the data by that each column being each specified value
      data_r <- data_r[which(data_r$Filter==sub_val_c),]
      colnames(data_r)[col_c] <- colname_c
    }
    #Now we focus on the columns of the measures and the groups for the statistics ignoring the remaining columns
    data_r_filt <- data_r[,which(colnames(data_r) %in% c(meas_r, gp_var))]
    
    #We change the column names to apply to all measures and group names for the stats as we go through the stats table
    colnames(data_r_filt)[which(colnames(data_r_filt)==meas_r)] <- "Amount"
    colnames(data_r_filt)[which(colnames(data_r_filt)==gp_var)] <- "Group"
    #We ensure that the amount column is numeric for statistics
    data_r_filt$Amount <- as.numeric(data_r_filt$Amount)
    #If the data is normal we do an ANOVA with a Post Hoc Tukey's test for multiple comparisons
    if (normality) {
      #The ANOVA model is put in
      model_r <- aov(data=data_r_filt, Amount ~ Group)
      #The overall p value is pulled
      overall_p_r <- summary(model_r)[[1]][["Pr(>F)"]][1]
      #Now we do a Tukey's Multiple Comparisons
      tuk_r <- TukeyHSD(model_r)
      #Now we pull out the groups for matching
      tuk_terms_r <- tuk_r$Group
      #Now we pull the rows from the stat summary table
      stat_tab_rows <- which(stats_sum_table$Measure==meas_r)
      #We go to pull the rows by sequentially filtering based on the other subset variables 
      for (l in 2:ncol(subs)) {
        rows_l <- which(stats_sum_table[,l]==subs[r,l])
        stat_tab_rows <- intersect(stat_tab_rows,rows_l)
      }
      #Now we go through those rows in the stats table and add in corresponding values
      for (k in stat_tab_rows) {
        #We add in the overall p value
        stats_sum_table$`p-value-overall`[k] <- overall_p_r
        #We pull the pairwise groups
        group_1k <- stats_sum_table$group1[k]
        group_2k <- stats_sum_table$group2[k]
        #Then we pull the Tukey's test results for those two groups
        tuk_row <- which(rownames(tuk_terms_r)==paste0(group_1k,"-",group_2k) | rownames(tuk_terms_r)==paste0(group_2k,"-",group_1k))
        stats_sum_table$padj[k] <- tuk_terms_r[tuk_row,"p adj"]
        if (tuk_terms_r[tuk_row,"p adj"]< 0.001) {
          stats_sum_table$p.signif[k] <- "***"
          # If p-value is below not under 0.001 but is under 0.01 two stars are provided
        } else if (tuk_terms_r[tuk_row,"p adj"] < 0.01) {
          stats_sum_table$p.signif[k] <- "**"
          # If p-value is below not under 0.01 but is under 0.05 one star is provided
        } else if (tuk_terms_r[tuk_row,"p adj"] < 0.05) {
          stats_sum_table$p.signif[k] <- "*"
          # If p-value is below not under 0.05, "ns" provided to indicate not significant
        } else {
          stats_sum_table$p.signif[k] <- "ns"
        }
      }
    }
    #If the data is NOT normal we do a Kruskall Wallace with a Post Hoc Dunn's test for multiple comparisons
    if (!normality) {
      #The Kruskal Wallace Test is put in
      model_r <- kruskal_test(data=data_r_filt, Amount ~ Group)
      #The overall p value is pulled
      overall_p_r <- model_r$p
      #Now we do a Dunn's test for multiple comparisons
      dun_r <- dunn_test(Amount ~ Group,
               data=data_r_filt,
               p.adjust.method=adj)
      #Now we pull the rows from the stat summary table
      stat_tab_rows <- which(stats_sum_table$Measure==meas_r)
      #We go to pull the rows by sequentially filtering these to narrow in on rows corresponding rows in the stats table
      for (l in 2:ncol(subs)) {
        rows_l <- which(stats_sum_table[,l]==subs[r,l])
        stat_tab_rows <- intersect(stat_tab_rows,rows_l)
      }
      #Now we go through those rows in the stats table and add in corresponding values
      for (k in stat_tab_rows) {
        #We add in the overall p value
        stats_sum_table$`p-value-overall`[k] <- overall_p_r
        #We pull the pairwise groups
        group_1k <- stats_sum_table$group1[k]
        group_2k <- stats_sum_table$group2[k]
        #Then we pull the Dunn's test results 
        dunn_row <- which((dun_r$group1==group_1k & dun_r$group2==group_2k) | (dun_r$group1==group_2k & dun_r$group2==group_1k))
        #We then enter the p values and p adjusted values
        stats_sum_table$`p-value`[k] <- dun_r$p[dunn_row]
        stats_sum_table$padj[k] <- dun_r$p.adj[dunn_row]
        #Based on the pairwise adjusted p value we add significance
        # If p-value is below 0.001 three stars are provided
        if (dun_r$p.adj[dunn_row]< 0.001) {
          stats_sum_table$p.signif[k] <- "***"
          # If p-value is below not under 0.001 but is under 0.01 two stars are provided
        } else if (dun_r$p.adj[dunn_row] < 0.01) {
          stats_sum_table$p.signif[k] <- "**"
          # If p-value is below not under 0.01 but is under 0.05 one star is provided
        } else if (dun_r$p.adj[dunn_row] < 0.05) {
          stats_sum_table$p.signif[k] <- "*"
          # If p-value is below not under 0.05, "ns" provided to indicate not significant
        } else {
          stats_sum_table$p.signif[k] <- "ns"
        }
      }
    }
  }
#Now we return the final stats table
return(stats_sum_table)
} 