extend.set.data=function(set.data,sep)

{
  set.id=strsplit(set.data$set.id,split=sep,fixed=T)
  set.id=lapply(set.id,unique)
  k=lapply(set.id,length)
  k=unlist(k)
  row.id=rep(set.data$row.id,k)
  mtx.id=rep(set.data$mtx.id,k)
  set.id=unlist(set.id)
  res=cbind.data.frame(set.id=set.id,
                       mtx.id=mtx.id,
                       row.id=row.id)
  return(res)

}
