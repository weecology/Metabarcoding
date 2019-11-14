### JV26 & JV38 Data Summary

### Each amplicon summary folder contains:
	1) Mapping folder
		A) Basic information about each sample.
	2) OTU.Table folder
		A) OTU table in BIOM format
		B) OTU table in TAB-delimited format
		C) General OTU table summary
		D) If BLASTn was used, any OTU that did not meet taxonomy assignment criteria
			will be denoted with "Unclassified". If no BLASTn hit was returned at all, 
			then another set of files ('Unassigned') was written.
		E) If Closed-Reference OTU folder, if ran.
	3) OTU.Seqs folder
		A) The representative OTU sequences. Only, if `de novo` OTU clustering 
			was performed.

		
### OTU picking and taxonomy assignment by amplicon:
	1) ITS2:
		A) Used r1 reads from JV26 & JV38
			i) primer sequences were removed
			ii) trimmed to 150 bp
		B) Clustered OTUs at 99% similarity
			i) used de novo based chimera removal
			ii) Classified using BLASTn
		c) Taxonomy is based on the top BLASTn hit that has at least a 90% query cover
			and 85% identity.
	