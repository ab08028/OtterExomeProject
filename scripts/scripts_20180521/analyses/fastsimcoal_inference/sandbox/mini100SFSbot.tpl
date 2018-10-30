//Parameters for the coalescence simulation program : simcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NEURCUR
NASNCUR
//Samples sizes and samples age
100
100
//Growth rates: negative growth implies population expansion
EURGROWRATE
ASNGROWRATE
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
5 historical events
TASNGROW 1 1 0 1 0 0
TEURGROW 0 0 0 1 0 0
TDIV 0 1 1 REC_EUR_RESIZE 0 0
TBOT 1 1 BOT_RESIZE 0 0
TANC 1 1 NANC_RESIZE 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data
FREQ 1 0 2.5e-8 OUTEXP

