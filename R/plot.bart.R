plot.aftree = function(
   x,
   plquants=c(.05,.95), cols =c('blue','black'),
   ...
)
{
     ql <- apply(x$m.train,2,quantile,probs=plquants[1])
     qm <- apply(x$m.train,2,quantile,probs=.5)
     qu <- apply(x$m.train,2,quantile,probs=plquants[2])

     ii <- order(qm)
     n.obs <- length(qm)
     print(n.obs)
     plot(qm[ii], 1:n.obs,xlim=range(ql,qu), ylim=c(1,n.obs), xlab='E(log T|x)',main=
      '90 percent posterior intervals for E(log T|x)', ylab = 'Observation sorted by E(log T|x)', type='n',las=1, ...)
     for (i in 1:n.obs) {
       lines(c(ql[ii[i]], qu[ii[i]]), c(i,i),col=cols[1])
       points(qm[ii[i]], i, pch=16)
     }
}
