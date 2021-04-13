#!/usr/bin/env Rscript

library(argparse)
suppressMessages(suppressWarnings(require(data.table)))
suppressMessages(suppressWarnings(require(rtracklayer)))


parser <- ArgumentParser(description='Process counts into normalized means for a track')
parser$add_argument('--input', help='counts of individual experiments')
parser$add_argument('--chrom_sizes', help='chrom sizes')
parser$add_argument('--out', help='output')
parser$add_argument('--binsize', help='binsize used')
parser$add_argument('--threads', help='number of threads', type="integer")

argv <- parser$parse_args()

setDTthreads(threads=argv$threads)
save(argv, file='test.Rdata')

input_dt = fread(argv$input, sep='\t')

colnames(input_dt) = gsub("[#]*'(.*)'", '\\1', colnames(input_dt))

na_sum = rowSums(input_dt[,-c('chr', 'start', 'end')])

input_dt = input_dt[!is.na(na_sum) & na_sum > 0, ]


track_gr = makeGRangesFromDataFrame(input_dt)

count_dt = input_dt[,-c('chr', 'start', 'end')] * as.numeric(argv$binsize)

norm <- function(exp, ctrl){
    n_vec = c(sum(exp), sum(ctrl))
    m_min = min(n_vec)
    exp=exp/n_vec[1] * m_min + 1
    ctrl=ctrl/n_vec[2] * m_min + 1
    return(exp / ctrl)
}

norm_dt = count_dt[,list(r1=norm(exp1, ctrl1),
                         r2=norm(exp2, ctrl2))]

score(track_gr) = rowMeans(log2(norm_dt))

seql_df = read.table(argv$chrom_sizes)

seq_vec = seql_df[,2]
names(seq_vec) = seql_df[,1]

seqlengths(track_gr) = seq_vec[seqlevels(track_gr)]

start(track_gr) = start(track_gr) + 1
export.bw(trim(track_gr), argv$out)
