vep_annotation () {
docker run --rm -v /path:/path:Z -v /path/:/path/ -w /data -p port:port   lab/vep104_root \
/opt/vep/src/ensembl-vep/vep --cache \
--refseq \
-i /path/esomi_project_backedup/joint_call/vqsr/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier.vcf.gz \
-o /path/esomi_project_backedup/annotation/vep/annotations_2022/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier_annotated_canonic.vcf.gz \
--force_overwrite \
--use_transcript_ref \
--dir_cache /opt/vep/.vep \
--dir_plugins /opt/vep/.vep/Plugins \
--species homo_sapiens \
--use_transcript_ref \
--assembly GRCh38       \
--no_stats \
--fork 18 \
--offline \
--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--buffer_size 500 \
--vcf \
--compress_output bgzip \
--biotype \
--hgvs \
--symbol \
--canonical \
--regulatory \
--gene_phenotype \
--variant_class \
--plugin DisGeNET,file=/path/share/vep_data/DisGeNET/all_variant_disease_associations_final.tsv.gz \
--plugin CADD,/path/share/vep_data/CADD/whole_genome_SNVs.tsv.gz \
--plugin Mastermind,/path/share/vep_data/Mastermind/mastermind_cited_variants_reference-2021.08.03-grch38.vcf.gz,0,0,1 \
--plugin dbNSFP,/path/share/vep_data/dbNSFP/dbNSFP4.2a_grch38.gz,REVEL_score,REVEL_rankscore,SIFT_pred,FATHMM_pred,PROVEAN_pred,MutationAssessor_pred,VEST4_score,LRT_score,LRT_pred,MutationTaster_pred,Polyphen2_HVAR_pred \
--plugin dbscSNV,/path/share/vep_data/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
--plugin SpliceAI,snv=/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/home/asselta/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
--plugin MaxEntScan,/path/share/vep_data/MaxEntScan/fordownload,SWA
}
#vep_annotation

vep_annotation () {
docker run --rm -v /path:/path:Z -v /path/:/path/ -w /data -p port:port   lab/vep104_root \
/opt/vep/src/ensembl-vep/vep --cache \
--refseq \
-i /path/esomi_project_backedup/joint_call/vqsr/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier.vcf.gz \
-o /path/esomi_project_backedup/annotation/vep/annotations_2022/allChr_95_95recalibrated_demultiplex_max_mc_id_no_outlier_annotated_canonic_pop.vcf.gz \
--force_overwrite \
--use_transcript_ref \
--dir_cache /opt/vep/.vep \
--dir_plugins /opt/vep/.vep/Plugins \
--species homo_sapiens \
--use_transcript_ref \
--assembly GRCh38       \
--no_stats \
--fork 18 \
--offline \
--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
--buffer_size 500 \
--vcf \
--compress_output bgzip \
--biotype \
--af_1kg \
--af_gnomad \
--hgvs \
--symbol \
--canonical \
--regulatory \
--gene_phenotype \
--variant_class \
--plugin DisGeNET,file=/path/share/vep_data/DisGeNET/all_variant_disease_associations_final.tsv.gz \
--plugin CADD,/path/share/vep_data/CADD/whole_genome_SNVs.tsv.gz \
--plugin Mastermind,/path/share/vep_data/Mastermind/mastermind_cited_variants_reference-2021.08.03-grch38.vcf.gz,0,0,1 \
--plugin dbNSFP,/path/share/vep_data/dbNSFP/dbNSFP4.2a_grch38.gz,REVEL_score,REVEL_rankscore,SIFT_pred,FATHMM_pred,PROVEAN_pred,MutationAssessor_pred,VEST4_score,LRT_score,LRT_pred,MutationTaster_pred,Polyphen2_HVAR_pred \
--plugin dbscSNV,/path/share/vep_data/dbscSNV/dbscSNV1.1_GRCh38.txt.gz \
--plugin SpliceAI,snv=/path/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/home/asselta/share/vep_data/spliceAI/files/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz,cutoff=0.5 \
--plugin MaxEntScan,/path/share/vep_data/MaxEntScan/fordownload,SWA
}
vep_annotation

#### notes ####
#--fasta /opt/vep/.vep/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
#fork 4 is the advised value for multithreding here: https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#faster
#buffer size 5000 is the advised value

#-entrypoint /usr/bin/perl
#--fasta  /path/ccappadona/reference/vep_references_refseq/homo_sapiens_refseq/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \

#-o /path/esomi_project_backedup/annotation/vep/annotations/CCRL2_annotated.vcf \



exit
