# `BCB420.2019.BioGRID`

#### (BioGRID data annotatation of human genes)

&nbsp;

## Installation

###### [Alison Wu] &lt;alison.wu@mail.utoronto.ca&gt;

## 1 About this package:


This package describe the workflow to download interaction datasets from [BioGRID](https://thebiogrid.org/),how to map the Entrez Gene IDs to HGNC symbols, and how to annotate the example gene set.


&nbsp;

#### In this project ...

```text
 --BCB420.2019.BioGrid/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.BioGrid.Rproj
   |__DESCRIPTION
   |__dev/
      |__rptTwee.R
      |__toBrowser.R              
   |__inst/
      |__extdata/
         |entrez2sym.RData         # Entrez ID to HGNC symbol mapping tool
         |__xSetInt.tsv          # annotated example edges
      |__img/
         |__[...]                  # image sources for .md document
      |__scripts/
         |__recoverIDs.R           # utility to use biomaRt for ID mapping
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R
   |__README.md                    # this file

```

&nbsp;

----

## 2 BioGRID Data

BioGRID is a curated biological database of interactions. BioGRID interactions include protein-protein interactions, genetic interactions, chemical interactions, and post-translational modifications.


&nbsp;

## 3 Data download and cleanup

To download the source data from BioGRID :

1. Navigate to the [**BioGRID** database](https://thebiogrid.org/) 
2. Click "downloads" in the left corner, click "Current Release" and choose "BIOGRID-3.5.169" or follow the link [here](https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-3.5.169/)
3. Download the following data file: BIOGRID-ORGANISM-3.5.169.tab2.zip (68.09 MB)
4. Unzip the file
5. Choose the file `BIOGRID-ORGANISM-Homo_sapiens-3.5.169.tab2.txt` (122 MB) and place it into the your local working directory which is called `data`.

&nbsp;

## 4 Mapping Entrez IDs to HGNC symbols

BioGrid genes are Entrez IDs. These IDs can be mapped to HGNC symbols.

#### Preparations: packages, functions, files
Install a few required packages to begin:

1. **`readr`** allows the program to read large datasets.
```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```

2. **`biomaRt`** is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```

3. **`igraph`** is to compute statistics graphs and plots.
&nbsp;

```R
if (! requireNamespace("igraph")) {
  install.packages("igraph")
}
```

&nbsp;

#### 4.1 Step one: which IDs do we have to map?
&nbsp;

```R
  # Read the interaction data
  tmp <-readr::read_delim(file.path("../data","BIOGRID-ORGANISM-Homo_sapiens-3.5.169.tab2.txt"),
                          delim = "\t",
                          skip = 1, col_names = c("BioGRID Interaction ID", "Entrez_ID1", "Entrez_ID2"))  #468,058 rows
                          
  
  head(tmp)
  #    A tibble: 6 x 3
  #     `BioGRID Interaction ID` Entrez_ID1 Entrez_ID2
  #     1             <dbl>      <dbl>      <dbl>
  #     2              103       6416       2318
  #     3              117      84665         88
  #     4              183         90       2339
  #     5              278       2624       5371
  #     6              418       6118       6774
  #     7              586        375      23163                        
                          

  # The number of unique IDs we have to map:
  uENTREZ <- unique(c(tmp$Entrez_ID1, tmp$Entrez_ID2)) #23108 IDs
  
  ```

&nbsp;

#### 4.2  Step two: mapping via biomaRt

We first map Entrez IDs to HGNC symbols in bulk using biomaRt. And then we can look for the remaining IDs that we could not map via UniProt IDs (from HGNC reference data).


&nbsp;

###### 4.2.1  Constructing an ID-mapping tool

```R
# Map ENSP to HGNC symbols: open a "Mart" object ..
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  tmp <- biomaRt::getBM(filters = "entrezgene",
                       attributes = c("entrezgene",
                                      "hgnc_symbol"),
                       values = uENTREZ,
                       mart = myMart)

  head(tmp)
  #      entrezgene hgnc_symbol
  # 1    10036      CHAF1A
  # 2    10048      RANBP9
  # 3    10059       DNM1L
  # 4    10084       PQBP1
  # 5    10128      LRPPRC
  # 6    10138        YAF2

  nrow(tmp)  #  17058 HGNC symbols have been retrieved for the 23108 Entrez IDs
  
```
&nbsp;

There are potential problems with the data that biomart returns

There might be duplicates (more than one value returned).
```R
sum(duplicated(tmp$entrezgene)) #96 duplicates
```

Or nothing returned for one Entrez ID
```R

  sum(! (uENTREZ) %in% tmp$entrezgene)  # 6146
```

&nbsp;

To fix the duplicates:

```R
  dupEntrez <- tmp$entrezgene[duplicated(tmp$entrezgene)]
  tmp[tmp$entrezgene %in% dupEntrez, ]
          entrezgene     hgnc_symbol
  # 321         5888           RAD51
  # 322         5888
  # 368         6606            SMN2
  # 369         6606            SMN1
  # 627        25788            FSBP
  # 628        25788          RAD54B
  # 936         8878          SQSTM1
  # 937         8878
  # 1566        2222           FDFT1
  # 1567        2222
  # 1717        5430          POLR2A
  # 1718        5430
  # 2049       23370        ARHGEF18
  # 2050       23370
  # 2919       11068        CYB561D2

  #Remove the symbols with space
  tmp[tmp$hgnc_symbol %in% c(""), ]
  tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("")), ]
  
  #Check the UniprotID for the duplicates one
  #SMN2 and SMN1 has the same UniprotID Q16637 - so assign arbitarily to either one
  tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("SMN2","FSBP", "HLA-DQA2", "LINC02210-CRHR1", "DEFB4B", "DEFB103A", "CCL3L1", "CRYAA2", "IQCJ-SCHIP1", "TBC1D3L", "SERF1B", "DDT", "STAG3L3", "ZNF468", "KIR2DS5", "FNTB", "CHURC1-FNTB", "GTF2H2", "LGALS7B", "PLEKHG7", "FAM187A", "CBWD3", "BOLA2B", "RPP21", "TRIM39-RPP21", "LGALS7", "GUSBP1", "LINC00680","ARL17B", "USP17L2", "GATD3A", "RNA5S3", "PRR4", "GAGE12F", "GAGE12G")), ]
  
  #Check duplicates
  any(duplicated(tmp$entrezgene)) #FALSE
  ```
  
  &nbsp;
  
  Define the mapping tool after this cleanup
  
  &nbsp;
  
  ```R
  entrez2sym <- tmp$hgnc_symbol
  names(entrez2sym) <- tmp$entrezgene

  head(entrez2sym)
  # 10036    10048    10059    10084    10128    10138
  #"CHAF1A" "RANBP9"  "DNM1L"  "PQBP1" "LRPPRC"   "YAF2"
  ```
  &nbsp;
  
  ###### 4.2.2  Cleanup and validation of `entrez2sym`
  
  ```R
  sel <- ! (uENTREZ %in% names(entrez2sym))
> x <- rep(NA, sum( sel))
> names(x) <- uENTREZ[ sel ]

> any(duplicated(c(names(x), names(entrez2sym))))  #FALSE
entrez2sym <- c(entrez2sym, x)
> all(uENTREZ %in% names(entrez2sym)) #TRUE

-check if there are empty strings
sel <- which(entrez2sym == "")
> length(sel) #0

```
&nbsp;

Add out dated symbols

```R

  sel <- ( ! (entrez2sym %in% HGNC$sym)) & ( ! (is.na(entrez2sym)))
  length(        entrez2sym[ sel ] ) #175
  length( unique(entrez2sym[ sel ])) #175
  unkSym <- data.frame(unk = entrez2sym[ sel ],
                       new = NA,
                       stringsAsFactors = FALSE)

  View(unkSym)
  
  #add the outdated symbols(which are usually stored in the prev column in HGNC to entrez2sym)
  for (i in seq_len(nrow(unkSym))) {
      iPrev <- grep(unkSym$unk[i], HGNC$prev)[1]
      if (length(iPrev) == 1) {
          unkSym$new[i] <- HGNC$sym[iPrev]
 }
 ```
 
 4.3 Final Validation
 
 Final validation for our mapping tool
 
 
```R
  all(uENTREZ %in% names(entrez2sym)) #TRUE, all the Entrez IDs are mapped
  
  sum(! is.na(entrez2sym)) #16916 symbols are found
  
  sum(! is.na(entrez2sym)) * 100 / length(entrez2sym) #73.20%
  
  #Save map
  save(entrez2sym, file = file.path("inst", "extdata", "entrez2sym.RData"))
  
  #Load file from RStudio project
  load(file = file.path("inst", "extdata", "entrez2sym.RData"))
  ```
  

&nbsp;

# 5 Annotating gene sets with STRING Data

We can now annotate the gene sets with BioGrid dta using our mapping tool. 


&nbsp;

```R
  tmp <-readr::read_delim(file.path("../data", "BIOGRID-ORGANISM-Homo_sapiens-3.5.169.tab2.txt"),
                          delim = "\t",
skip = 1, col_names = c("BioGRID Interaction ID", "Entrez_ID1", "Entrez_ID2", "ID1", "ID2", "sys_name1", "sys_name2", "sym1", "sym2", "alias1", "alias2", "exp_name", "interaction_type", "author", "pubID", "org_ID1", "org_ID2", "interaction_throughput", "quant_score", "mod", "pheno", "qual", "tags", "ext"))

 #the number of interaction that is genetic
 nrow(tmp[tmp$interaction_type == "genetic",]) #5299
 
 #the number of interaction that is physical
 nrow(tmp[tmp$interaction_type == "physical",]) #462759
 
 #Finally we map Entrez IDs to HGNC sysmbols
tmp$Entrez_ID1 <- entrez2sym[tmp$Entrez_ID1]
tmp$Entrez_ID2 <- entrez2sym[tmp$Entrez_ID2]

VALIDATE:
sum(is.na(tmp$Entrez_ID1)) #168627
sum(is.na(tmp$Entrez_ID2)) #169208
# we remove edges in which either one or the other node is NA to
# create our final data:
BioGrid <- tmp[( ! is.na(tmp$Entrez_ID1)) & ( ! is.na(tmp$Entrez_ID2)), ]
 BioGridInt <- tmp[( ! is.na(tmp$Entrez_ID1)) & ( ! is.na(tmp$Entrez_ID2)), ]
#DOne
save(BioGridInt, file = file.path(".", "data", "BioGridInt.RData"))
 
 
 #Network Graph
 t <- data.frame(Num_GeneticInteractions = num_genetic, Percentage=round(num_genetic / total * 100, 2), Num_PhysicalInteractions= num_phy,Percentage=round(num_phy / total * 100, 2),stringsAsFactors = FALSE)
  head(t)
  Num_GeneticInteractions Percentage Num_PhysicalInteractions
1                    5299       1.13                   462759
  Percentage.1
1        98.87


 #the distribution
 deg <-  table(c(BioGridInt$Entrez_ID1, BioGridInt$Entrez_ID2))
 
  hist(deg, breaks=50,
      xlim = c(0, 1400),
      col = "#3fafb388",
      main = "STRING nodes degree distribution",
      xlab = "degree (undirected graph)",
      ylab = "Counts")
      ```
![](./inst/img/BioGridDistribution.png?sanitize=true "BioGRID Interaction distribution") 

```


&nbsp;



####  6 Annotation of the example gene set

To conclude, we annotate the example gene set, validate the annotation, and store the data.

&nbsp;

```R

# The specification of the sample set is copy-paste from the 
# BCB420 resources project.

xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")


# Example genes are not among the known nodes?
x <- which( ! (xSet %in% c(BioGridInt$Entrez_ID1, BioGridInt$Entrez_ID2)))
> cat(sprintf("\t%s\t(%s)\n", HGNC[xSet[x], "sym"], HGNC[xSet[x], "name"]))


# ATP2A1	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 1)
# 	ATP2A2	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 2)
# BECN2	(beclin 2)
# BIRC6	(baculoviral IAP repeat containing 6)
# BLOC1S2	(biogenesis of lysosomal organelles complex 1 subunit 2)
# BORCS5	(BLOC-1 related complex subunit 5)
# BORCS7	(BLOC-1 related complex subunit 7)
# BORCS8	(BLOC-1 related complex subunit 8)
# CTTN	(cortactin)
# EPG5	(ectopic P-granules autophagy protein 5 homolog)
# GABARAPL2	(GABA type A receptor associated protein like 2)
# INPP5E	(inositol polyphosphate-5-phosphatase E)
# IRGM	(immunity related GTPase M)
# LAMP5	(lysosomal associated membrane protein family member 5)
# MAP1LC3B	(microtubule associated protein 1 light chain 3 beta)
# MGRN1	(mahogunin ring finger 1)
# OPTN	(optineurin)
# OSBPL1A	(oxysterol binding protein like 1A)
# PI4K2A	(phosphatidylinositol 4-kinase type 2 alpha)
# PIK3C3	(phosphatidylinositol 3-kinase catalytic subunit type 3)
# RAB20	(RAB20, member RAS oncogene family)
# RAB29	(RAB29, member RAS oncogene family)
# RAB34	(RAB34, member RAS oncogene family)
# RAB7A	(RAB7A, member RAS oncogene family)
# RPTOR	(regulatory associated protein of MTOR complex 1)
# RUBCNL	(rubicon like autophagy enhancer)
# SNAP29	(synaptosome associated protein 29)
# SNAP47	(synaptosome associated protein 47)
# SPG11	(SPG11, spatacsin vesicle trafficking associated)
# STX17	(syntaxin 17)
# STX6	(syntaxin 6)
# TFEB	(transcription factor EB)
# TGM2	(transglutaminase 2)
# TOM1	(target of myb1 membrane trafficking protein)
# TPCN1	(two pore segment channel 1)
# TPCN2	(two pore segment channel 2)
# TPPP	(tubulin polymerization promoting protein)
# VPS11	(VPS11, CORVET/HOPS core subunit)
# YKT6	(YKT6 v-SNARE homolog)

#select the gene  that are in the example set:
sel <- (BioGridInt$Entrez_ID1 %in% xSet) & (BioGridInt$Entrez_ID2 %in% xSet)
xSetInt <- BioGridInt[sel, c("Entrez_ID1", "Entrez_ID2")]

#Statitics:
nrow(xSetInt) #24

# Save the annotated set
writeLines(c("Entrez_ID1\tb",
                            sprintf("%s\t%s", xSetInt$Entrez_ID1, xSetInt$Entrez_ID2)),
                        con = "xSetInt.tsv")
# The data set can be read back
myXset <- read.delim(file.path("inst", "extdata", "xSetInt.tsv"),
                      stringsAsFactors = FALSE)


myXset <- read.delim(system.file("extdata",
                                  "xSetInt.tsv",
                                  package = "BCB420.2019.BioGrid"),
                      stringsAsFactors = FALSE)
 ```

&nbsp;
