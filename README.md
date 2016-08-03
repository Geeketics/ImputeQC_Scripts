# ImputeQC_Scripts
Scripts for various steps of the QC after whole-genome imputation

Summary of each files purpose:

  - 1_strand_align.sh: used by current workflow script to strand align
  - Current_workflow_pos_nofilter.sh: Murray's script, use with pos, doesn't filter on missingness (for sequence-based SNPs) - out of date
  - Current_workflow_pos.sh: Murray's script, use with pos - out of date
  - Current_workflow.sh: Murray's script - choose rsID or pos, choose seq or chip, choose phase or not
  - Current_workflow_wHWE.sh: script used to prep coreExome files for impute (does not HWE filter - already done previously), made from Current_workflow.sh
  - extractVCFinfo2.sh: script to get first 8 columns of VFC file - compensates for vcfs that specify if snp was 'typed'
  - extractVCFinfo.sh: script to get first 8 columns of VCF file
  - hwe_afterimpute_running.txt: summary of commands run
  - hwe_afterimpute.sh: calculates hwe, & freq per snp, outputs list of extreme hwe departures
  - hwe_filter.R: script to filter hwe p-value, specify multiple testing correction method
  - imputationinfo_bychr.R - manhattan plot of infoscores split by chromosome
  - imputationinfo.R - manhattan plot of infoscores
  - InfoScore_Filtering.txt: summary of commands run
  - merge_chromosome_summary.R: merge all per chromosome maf/hwe outputs & plot info score vs maf/hwe bins - makes a massive file, takes ages, defunct idea
  - merge_info_hwe2.R: merge info, hwe & freq files
  - merge_info_hwe.R: merge info & hwe files - out of date
  - new_v4_ComparePOS1000GWithBimFiles.r: used by current workflow script to strand align
  - new_v4_CompareRS100GWithBimFiles.r: used by current workflow script to strand align
  - plot_info.R: plots info score histograms and vs maf/hwe plots - out of date
  - plot_infoscore.R: plot info score vs maf/hwe bins - outputs summary file used in plotting
  - README_StrandFixing.txt: Murray's notes
  - renameSNP.sh: rename SNPs labelled "." in VCF file to be chr_pos_ref_alt - essential step
