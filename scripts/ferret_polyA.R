library(ggplot2)
library(tidyverse)
library(dplyr)
library(stats)
#install.packages("ggdist")

poly_a_d <- read.csv("polyA_with_gene_assigned_all_type.csv", header=TRUE)

#unique(poly_a_d$stype)
poly_a_d$join <- paste(poly_a_d$sample, poly_a_d$stype, sep="-")

unique(poly_a_d$join)
x <- poly_a_d[poly_a_d$join=="vivo-24H","polyA_length"]
y <- poly_a_d[poly_a_d$join=="vivo-Mock","polyA_length"]
wilcox.test(x, y) 
wilcox.test(y, b)

z <- poly_a_d[poly_a_d$join=="vitro-Mock", "polyA_length"]
wilcox.test(y, z)

#vitro-mock vs vivo-6h
wilcox.test(z, b)
#vitro-mock vs vivo-24h
wilcox.test(z, x)

a <- poly_a_d[poly_a_d$join=="vitro-Treated", "polyA_length"]
b <- poly_a_d[poly_a_d$join=="vivo-6H", "polyA_length"]
wilcox.test(x, b)

wilcox.test(y, b)
#vitro-Treated vs vitro-Mock
wilcox.test(a, z)
#vitro-Treated vs vivo-Mock
wilcox.test(a, y)
#vitro-Treated vs vio-6H
wilcox.test(a, b)
#vitro-Treated vs vivo-24H
wilcox.test(a, x)

poly_a_d%>%
  group_by(join)%>% 
  summarise(Mean=mean(polyA_length), Max=max(polyA_length), Min=min(polyA_length), Median=median(polyA_length), Std=sd(polyA_length))
d_t <- poly_a_d %>% count(join)
d_sub <- poly_a_d[poly_a_d$polyA_length >= 400, ] %>% count(join)
d_sub$n/d_t$n*100
d_t

d_sub


sub_poly <- poly_a_d[sample(nrow(poly_a_d), size=500000),]
poly_a_d_mock <- poly_a_d[poly_a_d$stype=="Mock",]
mock_sub <- poly_a_d_mock[sample(nrow(poly_a_d_mock), size=100000),]

vivo_d <- poly_a_d[poly_a_d$sample=="vivo"&poly_a_d$gene_type=="Nuclear", ]
vivo_d_m <- poly_a_d[poly_a_d$sample=="vivo"&poly_a_d$gene_type!="Nuclear", ]

#random select a subset to make things quicker
vivo_sub <- vivo_d[sample(nrow(vivo_d), size=100000),]

vitro_d <- poly_a_d[poly_a_d$sample=="vitro"&poly_a_d$gene_type=="Nuclear",]
vitro_d_sub <- vitro_d[sample(nrow(vitro_d), size=100000),]

vitro_d_m <-poly_a_d[poly_a_d$sample=="vitro"&poly_a_d$gene_type!="Nuclear",]


ggplot(vivo_sub, aes(x = polyA_length, y = stype)) +
  ggdist::stat_halfeye(
    #adjust = .5, 
    #width = .6, 
    #.width = 0, 
    justification = -.2, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  )

ggplot(vivo_sub, aes(x = polyA_length, y = stype))+ 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA,
    color="orange"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  geom_point(
    ## draw horizontal lines instead of points
    shape = 95,
    size = 15,
    alpha = .2,
    color = "#1D785A"
  )

#ggplot(vivo_sub, aes(x = stype, y = polyA_length)) + 
ggplot(sub_poly, aes(x = gene_type, y = polyA_length)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.2, 
    point_colour = NA,
    fill="pink"
  ) + 
  geom_boxplot(
    width = .15, 
    outlier.shape = NA
  ) +
  ## add justified jitter from the {gghalves} package
  gghalves::geom_half_point(
    ## draw jitter on the left
    side = "l", 
    ## control range of jitter
    range_scale = .4, 
    ## add some transparency
    alpha = .2,
    color = "#1D785A"
  )
