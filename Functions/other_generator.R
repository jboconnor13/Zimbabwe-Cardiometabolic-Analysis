###This function generates a list of variables with "Other" for those that are either infrequent or low in abundance 
other_generator <- function(data, name_field, type, threshold, abundance_field="Abundance", id_var="Sample") {
library(dplyr)
#data=mdf_add
#name_field="Genus"
#type="Percentage"
#threshold=0.2
#abundance_field="Count"
#id_var="Sample"
nm_col <- which(colnames(data)==name_field)
ab_col <- which(colnames(data)==abundance_field)
id_col <- which(colnames(data)==id_var)
data_temp <- data[,c(id_col,nm_col,ab_col)]
colnames(data_temp) <- c("id","nm","ab")
  if (type=="Max Abundance") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(max = max(as.numeric(get(abundance_field))))
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$max)>threshold)]
  }
  if (type=="Max Percent Abundance") {
    max_data_temp <- data_temp
    max_data_temp["perc"] <- NA
    for (l in 1:nrow(max_data_temp)) {
      stud_l <- max_data_temp[l,1]
      sum_l <- sum(as.numeric(max_data_temp[which(as.vector(max_data_temp$id)==stud_l),3]))
      max_data_temp$perc <- as.numeric(max_data_temp$ab)/as.numeric(sum_l)
    }
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$perc)>threshold)]
  }
  if (type=="Top Average") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(avg = mean(as.numeric(ab))) %>%
      arrange(desc(avg))
    names_inc <- max_data_temp$nm[1:threshold]
  }
  if (type=="Top Median") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(med = median(as.numeric(ab))) %>%
      arrange(desc(med))
    names_inc <- max_data_temp$nm[1:threshold]
  }
  if (type=="Top Max") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(max = max(as.numeric(ab))) %>%
      arrange(desc(max))
    names_inc <- max_data_temp$nm[1:threshold]
  }
  if (type=="Average") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(avg = mean(as.numeric(ab))) %>%
      arrange(desc(avg))
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$avg)>threshold)]
  }
  if (type=="Median") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(med = median(as.numeric(ab))) %>%
      arrange(desc(med))
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$med)>threshold)]
  }
  if (type=="Max") {
    max_data_temp <- data_temp %>%
      group_by(nm) %>%
      summarise(max = max(as.numeric(ab))) %>%
      arrange(desc(max))
    names_inc <- max_data_temp$nm[which(as.numeric(max_data_temp$max)>threshold)]
  }
  if (type=="Frequency") {
    data_temp <- data_temp[which(data_temp[,3]>0),]
    id_num <- length(unique(data[,which(colnames(data)==id_var)]))
    freq_data_temp <-as.data.frame(table(data_temp[which(colnames(data_temp)=="nm")]))
    names_inc <- freq_data_temp[which(as.numeric(freq_data_temp$Freq)>threshold),1]
  }
  if (type=="Percentage") {
    data_temp <- data_temp[which(data_temp[,3]>0),]
    id_num <- length(unique(data[,which(colnames(data)==id_var)]))
    freq_data_temp <-as.data.frame(table(data_temp[which(colnames(data_temp)=="nm")]))
    names_inc <- freq_data_temp[which(as.numeric(freq_data_temp$Freq/id_num)>threshold),1]
  }
  new_varnames <- as.data.frame(matrix("Other",nrow(data),1))
  data_var <- data
  colnames(data_var)[which(colnames(data_var)==name_field)] <- "Test_Group_Other"
  non_other_ind <- which(data_var$Test_Group_Other %in% names_inc)
  other_ind <- setdiff(1:nrow(data), non_other_ind)
  new_varnames$V1[non_other_ind] <- data_var$Test_Group_Other[non_other_ind]

  return(new_varnames)
} 
