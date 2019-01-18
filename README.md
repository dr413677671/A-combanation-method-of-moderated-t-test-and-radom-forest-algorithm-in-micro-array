# A combanation method of moderated t-test and radom forest algorithm in micro-array

%% The pipeline consists of 5 main functions: DE transcripts identification, random forest model, functional genes annotation, sub-group distinguish, and GO/KEGG annotation.

Introduction
-------------

RNASeq is a weel-developed sequencing technology. The pipeline is an reliable tool to identify DE transcripts or probeset in affy-metrix microarray data between two cohorts.
Distinguish the cohorts based on a combanation method of moderated t-test and random forest algorithm.
Generated the GO/KEGG annotation
Compared the results of two groupings



As of R >= 3.4.4, the script require some R standard library. For users who still need to support several packages are required: 
pd.hta.2.0, stringr , GEOquery, simpleaffy, GenomicRanges, oligoClasses, oligo, oligoData, RColorBrewer, affyPLM, ggplot2, limma, topGO, clusterProfiler, Rgraphviz and randomForest


Compatibility
-------------

The pipeline should work on R>=3 (it was tested on 3.4.4).


Configration
------------

The pipeline need to be modified in-codes.
                
Bugs
----

If you find a bug, please try to reproduce it with R 3.4.4 and latest bioconductor.

If it does not happen in 3.6.5, please file a bug in the python issue tracker.
If it happens there also, file a bug to 413677671@qq.com.


Ran Duan
Email: 413677671@qq.com
