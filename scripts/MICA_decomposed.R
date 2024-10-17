fye<-function(oa, ea, cf, cs, wf) {
  if(oa>ea) {
    return(-log(sum(dhyper(oa:cf, cf, (cs-cf), wf)), 10))
  } else {
    return(-log(sum(dhyper(0:oa, cf, (cs-cf), wf)), 10))
  }
}

llr<-function(oa, ob, oc, od, ea, eb, ec, ed) {
  s1<-ifelse(log((oa/ea), base=exp(1))*oa=="NaN", 0, log((oa/ea), base=exp(1))*oa)
  s2<-ifelse(log((ob/eb), base=exp(1))*ob=="NaN", 0, log((ob/eb), base=exp(1))*ob)
  s3<-ifelse(log((oc/ec), base=exp(1))*oc=="NaN", 0, log((oc/ec), base=exp(1))*oc)
  s4<-ifelse(log((od/ed), base=exp(1))*od=="NaN", 0, log((od/ed), base=exp(1))*od)
  return(2*sum(s1, s2, s3, s4))
}




#call to be investigated
construction.name <- "Bark"
#freq of the call of interest
construction.freq <- 120
#size of the corpus
corpus <- 1000

which.index <-3


which.sort <- 1

which.accuracy <- 3

data<-read.table(file.choose(), header=T, sep="\t", quote="", comment.char=""); cases<-length(data[,1]); cat("\n")

which.output<-1


# computation

words<-data[,1]; word.freq<-data[,2]; obs.freq<-data[,3]; exp.freq<-faith<-delta.p.constr.to.word<-delta.p.word.to.constr<-relation<-coll.strength<-c(rep(0, cases))

for (i in 1:cases) {
  #i <- 1
  obs.freq.a<-obs.freq[i]
  obs.freq.b<-construction.freq-obs.freq.a
  obs.freq.c<-word.freq[i]-obs.freq.a
  obs.freq.d<-corpus-(obs.freq.a+obs.freq.b+obs.freq.c)
  
  exp.freq.a<-construction.freq*word.freq[i]/corpus; exp.freq[i]<-round(exp.freq.a, which.accuracy)
  exp.freq.b<-construction.freq*(corpus-word.freq[i])/corpus
  exp.freq.c<-(corpus-construction.freq)*word.freq[i]/corpus
  exp.freq.d<-(corpus-construction.freq)*(corpus-word.freq[i])/corpus
  
  faith[i]<-round((obs.freq.a/word.freq[i]), which.accuracy)
  
  delta.p.constr.to.word[i]<-round((obs.freq.a/(obs.freq.a+obs.freq.b))-(obs.freq.c/(obs.freq.c+obs.freq.d)), which.accuracy)
  delta.p.word.to.constr[i]<-round((obs.freq.a/(obs.freq.a+obs.freq.c))-(obs.freq.b/(obs.freq.b+obs.freq.d)), which.accuracy)
  
  coll.strength[i]<-round(switch(which.index,
                                 fye(obs.freq.a, exp.freq.a, construction.freq, corpus, word.freq[i]),
                                 llr(obs.freq.a, obs.freq.b, obs.freq.c, obs.freq.d, exp.freq.a, exp.freq.b, exp.freq.c, exp.freq.d),
                                 log((obs.freq.a/exp.freq.a), 2),
                                 (corpus*(((obs.freq.a)*((corpus-construction.freq-word.freq[i]+obs.freq.a)))-((construction.freq-obs.freq.a)*(word.freq[i]-obs.freq.a)))^2)/(construction.freq*word.freq[i]*((construction.freq-obs.freq.a)+((corpus-construction.freq-word.freq[i]+obs.freq.a)))*((word.freq[i]-obs.freq.a)+((corpus-construction.freq-word.freq[i]+obs.freq.a)))),
                                 log(((obs.freq.a+0.5)/(obs.freq.b+0.5))/((obs.freq.c+0.5)/(obs.freq.d+0.5)), 10)), which.accuracy)
  if (obs.freq.a>exp.freq.a) {
    relation[i]<-"attraction"
  } else if (obs.freq.a<exp.freq.a) {
    relation[i]<-"repulsion"
  } else {
    relation[i]<-"chance"
  }
}

output.table<-data.frame(call_1 = words, call_1_freq = word.freq, obs.freq, exp.freq, relation, faith, delta.p.call_2.to.call_1 = delta.p.constr.to.word, delta.p.call_1.to.call_2 = delta.p.word.to.constr, coll.strength)

attach(output.table)

sort.index<-switch(which.sort, order(call_1), order(-obs.freq, call_1), order(-faith, call_1), order(relation, -coll.strength))
output.table<-output.table[sort.index,]

detach(output.table)

# hypothetical repulsion strength of unattested verbs
corp.size<-as.integer(log(corpus, 10))
absents.words<-absents.obs.freqs<-absents.exp.freqs<-absents.delta.p.constr.to.word<-absents.delta.p.word.to.constr<-absents.collstrengths<-c(rep(0, corp.size))
for (i in 1:corp.size) {
  absents.words[i]<-letters[i]
  absents.obs.freqs[i]<-10^i
  
  obs.freq.a<-0
  obs.freq.b<-construction.freq
  obs.freq.c<-10^i
  obs.freq.d<-corpus-(construction.freq+10^i)
  
  exp.freq.a<-construction.freq*10^i/corpus; absents.exp.freqs[i]<-round(exp.freq.a, which.accuracy)
  exp.freq.b<-construction.freq*(corpus-10^i)/corpus
  exp.freq.c<-(corpus-construction.freq)*10^i/corpus
  exp.freq.d<-(corpus-construction.freq)*(corpus-10^i)/corpus
  
  absents.delta.p.constr.to.word[i]<-round((obs.freq.a/(obs.freq.a+obs.freq.b))-(obs.freq.c/(obs.freq.c+obs.freq.d)), which.accuracy)
  absents.delta.p.word.to.constr[i]<-round((obs.freq.a/(obs.freq.a+obs.freq.c))-(obs.freq.b/(obs.freq.b+obs.freq.d)), which.accuracy)
  
  absents.collstrengths[i]<-round(switch(which.index,
                                         fye(obs.freq.a, exp.freq.a, construction.freq, corpus, word.freq[i]),
                                         llr(obs.freq.a, obs.freq.b, obs.freq.c, obs.freq.d, exp.freq.a, exp.freq.b, exp.freq.c, exp.freq.d),
                                         log((obs.freq.a/exp.freq.a), 2),
                                         (corpus*(((obs.freq.a)*((corpus-construction.freq-word.freq[i]+obs.freq.a)))-((construction.freq-obs.freq.a)*(word.freq[i]-obs.freq.a)))^2)/(construction.freq*word.freq[i]*((construction.freq-obs.freq.a)+((corpus-construction.freq-word.freq[i]+obs.freq.a)))*((word.freq[i]-obs.freq.a)+((corpus-construction.freq-word.freq[i]+obs.freq.a)))),
                                         log(((obs.freq.a+0.5)/(obs.freq.b+0.5))/((obs.freq.c+0.5)/(obs.freq.d+0.5)), 10)), which.accuracy)
}

output.table.hyp<-data.frame(absents.words, absents.obs.freqs, absents.exp.freqs, "repulsion", absents.delta.p.constr.to.word, absents.delta.p.word.to.constr, absents.collstrengths)
colnames(output.table.hyp)<-c("absents.call_1", "absents.obs.freqs", "absents.exp.freqs", "relation", "absents.delta.p.call_2.to.call_1", "absents.delta.p.call_1.to.call_2", "absents.collstrengths")
cat("\a") # progress beep
