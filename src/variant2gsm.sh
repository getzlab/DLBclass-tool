MAF=$1
SEG=$2
BP=$3
SETID="freeze_Aug14_all_248"
SETSAMPLE="../gsm/freeze_Aug14_all.248.samples.txt"
OUTDIR="/tmp/"
BP="/tmp/tsvcat_workflow_ea91634e-6af4-4349-ba3a-51cc0eb776c2_call-tsvcat_task_freeze_Aug14_all.dRanger_SvABA_forBP.txt"
MAF="/tmp/tsvcat_workflow_0f5336b2-ef1a-4a94-bade-c90b74daf713_call-tsvcat_task_freeze_Aug14_all.maf_pon_filter_all_mutations.maf"
SEG="/tmp/tsvcat_workflow_3096c8fd-7a39-4211-a945-f51b64dab821_call-tsvcat_task_freeze_Aug14_all.GATK4_cnv_postprocessing_tumor_acs_log2CR.seg"
CNVBLACKLIST="~/Projects/DLBCL-Classifier/data_tables/additional_gsm_inputs/CNV.hg19.bypos.111213.CR1_event_added.bed"
GISTICARMS="~/Projects/DLBCL-Classifier/data_tables/additional_gsm_inputs/DLBCL_broad_significance.18Aug2024.tsv"
GISTICFOCALS="~/Projects/DLBCL-Classifier/data_tables/additional_gsm_inputs/DLBCL_focal_peaks.18Aug2024.tsv"

echo "python3 sv2gsm.py -i $SETID -s $SETSAMPLE -v $BP -o $OUTDIR "
#python3 sv2gsm.py --i $SETID -s $SETSAMPLE -v $BP -o $OUTDIR

echo "python3 maf2gsm.py --i $SETID -s $SETSAMPLE -m $MAF -o $OUTDIR"
#python3 maf2gsm.py -i $SETID -s $SETSAMPLE -m $MAF -o $OUTDIR

echo "python3 seg2gsm.py -i $SETID -s $SETSAMPLE -v $SEG -x $CNVBLACKLIST -a $GISTICARMS -f $GISTICFOCALS -o $OUTDIR"
#python3 seg2gsm.py -i $SETID -s $SETSAMPLE -v $SEG -x $CNVBLACKLIST -a $GISTICARMS -f $GISTICFOCALS -o $OUTDIR


SVGSM="/tmp/freeze_Aug14_all_248.21Aug2024.SV.GSM.tsv"
MAFGSM="/tmp/freeze_Aug14_all_248.21Aug2024.MAF.GSM.tsv"
CNVGSM="/tmp/freeze_Aug14_all_248.21Aug2024.CNV.GSM.tsv"
FEATUREORDER="../gsm/feature_order.19Aug2024.txt"

python3 combine2gsm.py -i $SETID -v $SEG -v $SVGSM -m $MAFGSM -c $CNVGSM -f $FEATUREORDER -o $OUTDIR
