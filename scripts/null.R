null5.mean <- mean(null$plot5$mean)
null5.sd <- sd(null$plot5$mean)
null5.es <- (CWM$CWM[5] - null5.mean)/null5.sd
null5.p <- quantile(null$plot5$mean, probs = c(0.05, 0.95))
null5.p
null5.sd
CWM$CWM[5]


pnorm(CWM$CWM[5], mean = null5.mean, sd = null5.sd, lower.tail=FALSE)
pnorm(CWM$CWM[9], mean = null9.mean, sd = null9.sd, lower.tail=TRUE)
null5.es
null9.es


nullCWM <- c(mean(null$plot1$mean),
             mean(null$plot2$mean),
             mean(null$plot3$mean),
             mean(null$plot4$mean),
             mean(null$plot5$mean),
             mean(null$plot6$mean),
             mean(null$plot7$mean),
             mean(null$plot8$mean),
             mean(null$plot9$mean),
             mean(null$plot10$mean),
             mean(null$plot11$mean),
             mean(null$plot12$mean),
             mean(null$plot13$mean),
             mean(null$plot14$mean),
             mean(null$plot15$mean))

hydroCWM$nullCWM <- nullCWM

getStats(hydroCWM, hydroCWM$nullCWM, nullCWM)

nullRange <- c(mean(null$plot1$range),
               mean(null$plot2$range),
               mean(null$plot3$range),
               mean(null$plot4$range),
               mean(null$plot5$range),
               mean(null$plot6$range),
               mean(null$plot7$range),
               mean(null$plot8$range),
               mean(null$plot9$range),
               mean(null$plot10$range),
               mean(null$plot11$range),
               mean(null$plot12$range),
               mean(null$plot13$range),
               mean(null$plot14$range),
               mean(null$plot15$range))

hydroCWM$nullRange <- nullRange

getStats(hydroCWM, hydroCWM$nullRange, nullRange)


blah <- lm(CWM ~ C_MDFM, data = hydroCWM)
blah.null <- lm(nullCWM ~ C_MDFM, data = hydroCWM)
plot(hydroCWM$CWM, hydroCWM$nullCWM)             
hist(hydroCWM$CWM, breaks=5)


t.test(hydroCWM$CWM, hydroCWM$nullCWM)
