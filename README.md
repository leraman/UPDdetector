## Arguments

`vcf.file = '/Users/leraman/PhD/BAF/batch1-gatk-haplotype-joint-annotated.vcf.gz'`

`singles = c('D1810466', 'D1810468', 'D1810470')` # use for separate sample testing; will detect LOHs

`trio = c()` # use for trio analysis: father, mother, child; will detect LOH and increased mendelian errors (MEs)

`output.dir = '/Users/leraman/PhD/BAF/output2'`

`cairo.bitmap = F` # some systems require cairo bitmap for plotting

`min.seq.depth = 15` # consider variants with at least a depth of min.seq.depth

`hist.bins = 100` # number of bins in beta value histograms

`bin.size = 5` # in MB, for LOH/ME calling

`min.support = 1` # number of variants required in bin for P-value calculation; for LOH/ME calling

`successive.bins.required = 3` # successive 'significant' bins required for actual significance of middle bin; for LOH/ME calling; odd number

`cytoband = '/Users/leraman/PhD/BAF/cyotband.hg38.txt'` # if given, filters out gvar and stalk