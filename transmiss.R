library(tidyverse)
# library(tidyr)
library(ggplot2)

characteristic <- "transmiss"
range <- c("none", "1pt25", "1pt5", "1pt75")
vacc <- c("novacc","vacc")
matched <- c("matched", "unmatched")
crossimm <- c("crossimm0pt75", "crossimm1")

combos <- crossing(var1 = characteristic, var2 = range, var3 = vacc, var4 = matched, var5 = crossimm)
combos <- combos %>%
    filter(var3 == "novacc") %>%
    unite("file", var1:var5) 
combos$infousfile <- sprintf("%s_newinfous.txt", combos$file)
combos$cumfile <- sprintf("%s_cum.txt", combos$file)


infous_data <- data.frame(time = numeric(), new_infous1 = numeric(), new_infous2 = numeric(), new_infous = numeric(), label = numeric())

for (i in combos$file) {
    table <- read.table(combos$infousfile[combos$file==i], header=T)
    table$label <- i
    
    infous_data <- rbind(infous_data, table)
}

infous_data <- infous_data %>%
    rename(time=TIME, new_infous1 = new_infous1.1, new_infous2 = new_infous2.1, new_infous = new_infous.1) %>% 
    separate(label, c("characteristic", "range", "vacc", "match", "crossimm"), sep="_")

ggplot(infous_data, aes(x=time, y=new_infous, colour=range)) +
    geom_line() +
    facet_grid(rows=vars(crossimm), cols=vars(match))

cum_data <- data.frame(time = numeric(), strain1only = numeric(), both_strains = numeric(), strain2only = numeric(), allinf = numeric(), label = numeric())

for (i in combos$file) {
    table <- read.table(combos$cumfile[combos$file==i], header=T)
    table$label <- i
    
    cum_data <- rbind(cum_data, table)
}

cum_data <- cum_data %>%
    rename(time=TIME, strain1only = strain1only.1, both_strains = both_strains.1, strain2only = strain2only.1, allinf = allinf.1) %>% 
    filter(time==300) %>% 
    separate(label, c("characteristic", "range", "vacc", "match", "crossimm"), sep="_") %>%
    pivot_longer(c(strain1only, both_strains, strain2only, allinf), names_to = "strain", values_to = "count") %>% 
    filter(strain!="allinf") 
    

ggplot(cum_data, aes(x=range, y=count, fill=strain)) +
    geom_bar(position="stack", stat="identity") +
    facet_grid(rows=vars(crossimm), cols=vars(match))