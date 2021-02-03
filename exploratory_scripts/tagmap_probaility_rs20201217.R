TTAA_sites = 19226699*3

chromsizes<-read.table("/home/r.schep/mydata/data/genomes/GRCh38/hg38_chromsizes2.txt",header=F, stringsAsFactors = F) %>% 
  filter(V1 != "chrY")

# make 2nd X
chrXX = chromsizes %>% filter(V1 == "chrX") %>% 
  mutate(V1 = "chrXX")

# add 2nd X 
chromsizesXX = chromsizes %>%
  bind_rows(chrXX)

# make alleles
chromsizesXXA = chromsizesXX %>% mutate(V1 = paste0(V1, "A"))
chromsizesXXB = chromsizesXX %>% mutate(V1 = paste0(V1, "B"))
chromsizesXXC = chromsizesXX %>% mutate(V1 = paste0(V1, "C"))

# bind alleles
chrom = chromsizesXXA %>% 
  bind_rows(chromsizesXXB, chromsizesXXC)


bag_of_chroms = chrom %>% 
  distinct(V1, V2) %>% # make sure there is no duplicates
  mutate(V2 = round(V2/1000000, 0),  # bin de genome
         reps = list(rep(V1, V2))) %>% #ok this is shady but it works
  pull(reps) # pull out the list chromosomes

bag_of_chroms = unlist(bag_of_chroms[1]) # just take one of the lists from above, this is correct

# check if they are the same: 
chrom %>% distinct(V1, V2) %>% 
  mutate(V2 = round(V2/1000000, 0)) 
table(bag_of_chroms)


# sampling                                      
dupli_test = replicate(10000, {
  dupli = all(table(sample(bag_of_chroms, 7)) == 1)
  # dupli =  all(base::duplicated(sample(bag_of_chroms, 7)) == F) # works also
})

# percent with duplication
table(dupli_test)[1]/10000



library(data.table)
sizes = fread('/DATA/scratch/usr/c.leemans/data/hg38/hg38.chrom.sizes',
              col.names=c('seq', 'length'))[grep('chr[0-9X]+$',seq), ]
draw <- function(seq_vec, length_vec, n, nploid=3, bin_size=10^6){
  seq = paste0(rep(seq_vec, each=nploid), '.', rep(1:nploid, length(seq_vec)))
  n_vec = rep(round(length_vec/bin_size), each=nploid)
  seq_list = lapply(1:length(seq), function(i){
    paste(seq[i], 1:i,sep='.')
  })
  draw = sample(unlist(seq_list), n, replace=T)
  return(sum(table(draw)[table(draw)>1]))
}
t_list = mclapply(1:1000, function(i){
  t = sum(unlist(sizes[, lapply(1:286, function(i){draw(seq, length, 7)})]))
  return(t)
}, mc.cores=20)
mean(unlist(t_list))