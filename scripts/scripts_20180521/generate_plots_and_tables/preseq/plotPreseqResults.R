data <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/preseq/A10_Elut_CA_NIC_1_SN1.sea_otter_23May2016_bS9RH.deduped.99.preseq.lcExtrap.txt",header=T)
data
ggplot(data,aes(x=TOTAL_READS,y=EXPECTED_DISTINCT))+
  geom_line()
head(data)
data2 <- read.table("/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/preseq/A11_Elut_OR_T207T.Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.preseq.lcExtrap.txt",header=T)
ggplot(data2,aes(x=TOTAL_READS,y=EXPECTED_DISTINCT))+
  geom_line()
