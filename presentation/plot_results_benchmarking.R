library(jsonlite)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(rcartocolor)

#initiate empty df with correct # columns and types
result <- as.data.frame(matrix(nrow=0,ncol=5)) %>% mutate(across(everything(), as.character))
colnames(result) <- c("name","real_time","cpu_time","time_unit","OS")
result$real_time <- as.numeric(result$real_time)
result$cpu_time <- as.numeric(result$cpu_time)

#OS tested
dirs <- list.dirs("./results_benchmark")

# make dataframes out of benchmark results
for (i in 2:length(dirs)){
  
  list <- list.files(dirs[i])
  
  for (j in 1:length(list)){
    result2 <- fromJSON(paste0(dirs[i],"/",list[j]))$benchmarks
    result2 <- result2 %>% 
      select(name,real_time,cpu_time,time_unit) %>%
      mutate(OS=gsub("\\./results_benchmark/","",dirs[i]))
    result <- bind_rows(result,result2)
    
  }
}

# plot Benchmark matrix assembly #----
result2 <- result %>% filter(grepl('Bench1', name))
result2$name <- gsub("Bench.*","",result2$name)
result2$name <- gsub("assembleMatrix*","",result2$name)
result2 <- result2 %>% #make long format
  pivot_longer(cols=real_time:cpu_time,
               names_to="type",
               values_to="time")
result2$type <- factor(result2$type, levels=c("real_time","cpu_time"))
result2$name <- factor(result2$name, levels=c("Dense","SparseMan","SparseCsr","SparseEigen"))

plot <- result2 %>% 
  ggplot(aes(x=OS, 
             y=time,
             fill=name)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_carto_d(palette = "Earth") +
  facet_wrap(~type) +
  theme_bw()+
  labs(fill="",x="",y="Time [ms]")
  
ggsave(plot, file="assembly_bench_result.png", width=15, height=6.5, units="cm")

# plot time evolution #----
result2 <- result %>% filter(grepl('Bench2', name))
result2$name <- gsub("Bench.*","",result2$name)
result2$name <- gsub("timeEvolution*","",result2$name)
result2 <- result2 %>% #make long format
  pivot_longer(cols=real_time:cpu_time,
               names_to="type",
               values_to="time")
result2$type <- factor(result2$type, levels=c("real_time","cpu_time"))
result2$name <- factor(result2$name, levels=c("Dense","SparseMan","SparseCsr","SparseEigen"))

plot <- result2 %>% 
  ggplot(aes(x=OS, 
             y=time,
             fill=name)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_carto_d(palette = "Earth") +
  facet_wrap(~type) +
  theme_bw()+
  labs(fill="",x="",y="Time [ms]")

ggsave(plot, file="timeevolution_bench_result.png", width=15, height=6.5, units="cm")

# plot assembly+full time evolution #----
result2 <- result %>% 
  filter(grepl('Bench3', name)) %>%
  separate(name, sep="/", into = c("name", "dt_index", "iterations"), remove=TRUE)
timesteps <- data.frame(dt_index = as.character(seq(1:10)-1),
                        dt = paste0("dt=",seq(1:10)/0.00001))
result2 <- full_join(result2,timesteps)
result2$name <- gsub("Bench.*","",result2$name)
result2$name <- gsub("ode*","",result2$name)
result2 <- result2 %>% #make long format
  pivot_longer(cols=real_time:cpu_time,
               names_to="type",
               values_to="time")
result2$type <- factor(result2$type, levels=c("real_time","cpu_time"))
result2$name <- factor(result2$name, levels=c("Dense","SparseMan","SparseCsr","SparseEigen"))
result2$dt <-factor(result2$dt, levels=paste0("dt=",seq(1:10)/0.00001))

plot <- result2 %>% 
  ggplot(aes(x=OS, 
             y=time,
             fill=name)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_carto_d(palette = "Earth") +
  facet_wrap(~type + dt) +
  theme_bw()+
  labs(fill="",x="",y="Time [s]")

#ggsave(plot, file="ode_bench_result.png", width=15, height=6.5, units="cm")

plot <- result2 %>% 
  filter(type=="real_time") %>%
  droplevels() %>%
  ggplot(aes(x=OS, 
             y=time,
             fill=name)) +
  geom_bar(stat="identity", position="dodge")+
  scale_fill_carto_d(palette = "Earth") +
  facet_wrap(~type + dt,nrow=2) +
  theme_bw()+
  labs(fill="",x="",y="Time [s]")

ggsave(plot, file="ode_bench_result.png", width=18, height=8.5, units="cm")

