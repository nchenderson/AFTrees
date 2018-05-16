
ConcIndex <- function(times, pred.times, status) {
   indec <- which(status==1)
  
   IndexMat <- combn(indec, 2)
   T1 <- times[IndexMat[1,]]
   T2 <- times[IndexMat[2,]]
    
   That1 <- pred.times[IndexMat[1,]]
   That2 <- pred.times[IndexMat[2,]]
    
   ans <- mean((T1 > T2)*(That1 > That2) + (T2 <= T1)*(That2 <= That1))
   return(ans)
}
