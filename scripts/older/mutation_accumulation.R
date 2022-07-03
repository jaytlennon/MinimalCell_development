########ROY: delete below if code is "bad"


# Make table for contingency analyses
#ma.snm.mat <- matrix(ma.snm$count, ncol = 2, byrow = F)#############This code is bad. Produces an incorrect matrix
ma.snm.mat <- matrix(ma.snm$count, ncol = 2, byrow = T)
colnames(ma.snm.mat) <- c("non-minimal","minimal")
rownames(ma.snm.mat) <- c("A:T to C:G","A:T to G:C","A:T to T:A", "C:G to G:C",
                          "C:G to T:A", "C:G to A:T")
ma.snm.tab <- as.table(ma.snm.mat)
ma.snm.tab.margins <- addmargins(ma.snm.tab)

# X-squared = 69.92, df = NA, p-value = 9.999e-05
ma.snm.chi <- chisq.test(ma.snm.tab, simulate.p.value = TRUE, B = 10000)

# Posthoc analysis
ma.snm.z <- as.data.frame(ma.snm.chi$stdres)
ma.snm.x2 <- ma.snm.z$Freq^2
ma.snm.p <- pchisq(ma.snm.x2, df = 1, lower.tail = FALSE)
ma.snm.p.adj <- p.adjust(ma.snm.p, method="BH")
ma.snm.post.hoc <- data.frame(ma.snm.z, ma.snm.x2, ma.snm.p, ma.snm.p.adj)
colnames(ma.snm.post.hoc) <- c("type", "strain", "z", "chi2", "p", "p.adj")

#         type   strain           z     chi2        p          p.adj
# 1  A:T to C:G  minimal  0.44851026  0.201161457 6.537850e-01 0.785
# 2  A:T to G:C  minimal -4.84321140 23.456696673 1.277572e-06 2.555145e-06
# 3  A:T to T:A  minimal  0.08639216  0.007463605 9.311547e-01 0.931
# 4  C:G to G:C  minimal -6.90461823 47.673752850 5.033862e-12 1.510159e-11
# 5  C:G to T:A  minimal  9.47528700 89.781063769 2.660284e-21 1.596170e-20
# 6  C:G to A:T  minimal -3.68512950 13.580179403 2.285864e-04 0.00034

# 7  A:T to C:G non-minimal -0.44851026  0.201161457 6.537850e-01 0.785
# 8  A:T to G:C non-minimal  4.84321140 23.456696673 1.277572e-06 2.555145e-06
# 9  A:T to T:A non-minimal -0.08639216  0.007463605 9.311547e-01 0.931
# 10 C:G to G:C non-minimal  6.90461823 47.673752850 5.033862e-12 1.510159e-11
# 11 C:G to T:A non-minimal -9.47528700 89.781063769 2.660284e-21 1.596170e-20
# 12 C:G to A:T non-minimal  3.68512950 13.580179403 2.285864e-04 0.00034
```