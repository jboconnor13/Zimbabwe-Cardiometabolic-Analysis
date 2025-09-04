###This function takes a dataset and generates a volcano plot data set that includes fold changes differences and p-values

## Inputs:

## data=data set
## var_names=fields in the data set to be used for calculations (essentially the dependent variables are seeing difference in)
## sep_field=the field defining two variables for which differences will be calculated across (separating field)
## paired=indication of whether the data is paired, meaning for each measurement in one sep_field value there is a corresponding measurement in the other sep field measure (i.e Subject before and after an intervention)
## normality=indication of whether the data is normally distributed (in which a parametric test will be performed) or not normally distributed (in which a nonparametric test will be performed)
## id_var=a defining feature for paired data indicating the corresponding values
## upfield=optional indication of the value that we want to be treated as a positive change (i.e if you want to see how much the variables change from not treated to treated, treated would be the upfield)
## adj=adjustment technique for the p-value changes

## Outputs:

## A table that includes rows corresponding to each variable (defined by var_names) with the following variables
## var_names=the correpsoning variable from var_names
## p-value=p value obtained from differential abundance test (parametric t test for normal data and nonparametric wilcoxon for data not normally distributed)
## mean_diff=mean difference in two groups (difference in averages for unpaired data and average of difference for paired data) 
## mean_fc=mean fold change in two groups (FC of averages for unpaired data and average of FC values for paired data) 
## log2_fc=log2 of the mean fold change
## log2_fc=log2 of the mean difference
## padj=adjusted p value

volc_data_generator <- function(data, var_names, sep_field, paired, normality, id_var="StudyID", upfield="XXXXX", adj="bonferroni") {
  
  #We set up the data file to include rows that are the variables we are assessing
  volc_data <- as.data.frame(var_names)
  #We add columns for the difference statistics
  #One column for the p-values
  volc_data["p-value"] <- NA
  #One column for the mean difference (good for transformed data)
  volc_data["mean_diff"] <- NA
  #One column for the fold change
  volc_data["mean_fc"] <- NA
  
  #Now we pull the separating field column number and values 
  #Column for the separating variable 
  sep_col <- which(colnames(data)==sep_field)
  #We not pull the values for the separating fields (only should be two)
  sep_vals <- unique(data[which(!is.na(data[,sep_col])),sep_col])
  
  #Now we go through each variable (rows of the volc data set) and fill in the statistics 
  for (j in 1:nrow(volc_data)) {
    #We pull the variable name 
    var_j <- volc_data[j,1]
    #Depending on if the data is paired or not we will structure the temporary data set (data_j) differently
    #We select the separating field and the the variable
    select_cols <- c(sep_field, var_j, id_var)
    #Now we filter the data to only focus on those columns
    data_filt_j <- data[,select_cols]
    #Now we pivot to a wider format for the statistics
    data_wide_j <- pivot_wider(data=data_filt_j, #We take the filtered data
                               names_from=sep_field, #Names for each column will come from the seperating field
                               values_from=var_j) #We pull the values from the variable of interest
    #If there is a defined upfield, we place it to the right in the data frame
    if (upfield !="XXXXX") {
      #The upper field is the one specified
      up_val <- upfield
      #The lower field is the other value not indicated 
      dn_val <- setdiff(sep_vals, up_val)
      #The other columns are pulled as well
      other_cols <- setdiff(colnames(data_wide_j),sep_vals)
      #Now we place the upfield to the right (last column)
      data_wide_j <- data_wide_j[,c(other_cols,dn_val,up_val)]
    }
    #The data is changed to a dataframe 
    data_wide_j <- as.data.frame(data_wide_j)
    #Column names will be adjusted as well for simplicity
    colnames(data_wide_j)[(ncol(data_wide_j)-1):ncol(data_wide_j)] <- c("Down","Up")
    #We also ensure there are numeric
    data_wide_j$Down <- as.numeric(data_wide_j$Down)
    data_wide_j$Up <- as.numeric(data_wide_j$Up)
    #If the data is not paired, the statistics are simple to calculate
    if (!paired) {
      #Because the data is not paired we just look at averages across the columns 
      #We calulate the fold change
      meanfc_j <- mean(data_wide_j$Up, na.rm=TRUE)/mean(data_wide_j$Down, na.rm=TRUE)
      #We calsulate the mean difference
      meandiff_j <- mean(data_wide_j$Up, na.rm=TRUE)-mean(data_wide_j$Down, na.rm=TRUE)
    }
    #If the data is paired calculations are slightly adjusted (we generate pairwise fc and differences and then average)
    if (paired) {
      #We generate a new column for the FC
      data_wide_j$FC <- data_wide_j$Up/data_wide_j$Down
      #We generate a new column for the differences
      data_wide_j$Difference <- data_wide_j$Up-data_wide_j$Down
      #Now we calculate the average FC
      meanfc_j <- mean(data_wide_j$FC, na.rm=TRUE)
      #And we calculate the average difference
      meandiff_j <- mean(data_wide_j$Difference, na.rm=TRUE)
    }
    #Now we generate the stats based on the normality
    #First we need to specify there are enough observations
    
    if (nrow(data_wide_j)>1) {
      if (normality) {
        #If the data is normally distributed, we do a parametric t-test
        test_j <- t.test(data_wide_j$Up, data_wide_j$Down, paired=paired)
        #We pull the p0value
        p_j <- test_j$p.value
      }
      if (!normality) {
        #If the data is not normally distributed, we do a nonparametric wilcoxon sign-rank test 
        test_j <- wilcox.test(data_wide_j$Up, data_wide_j$Down, paired=paired)
        p_j <- test_j$p.value
      }  
    } else {
      print(paste0("Warning: Insufficient observations for p-value calclation for ", var_j))
      p_j <- NA
    }
    
    #We enter the fc in the volc data set
    volc_data$mean_fc[j] <- meanfc_j
    #We enter the difference in the volc data set
    volc_data$mean_diff[j] <- meandiff_j
    #We enter the p-value in the volc data set
    volc_data$`p-value`[j] <- p_j
  }
  #We now calculate the log transformed FC
  volc_data$log2_fc <- log2(volc_data$mean_fc)
  #We calculate the log transformed difference
  volc_data$log2_diff <- log2(volc_data$mean_diff)
  #We now calultate the adjusted p values
  volc_data$padj <- p.adjust(volc_data$`p-value`, method=adj)
  
  #Now we return the volc data set
  return(volc_data)
}



