# Neop
## Python Package for Predicting Neoantigens and Obtaining Amino Acid Context from MAFs

### Installation

You can install either the last stable version or current version of neop using pip.

#### Stable version:

`pip install neop`

#### Current version:

Download this directory. Navigate to it and run the following command:

`pip install .`

#### Regardless of which version you install:

Before beginning predictions, you must run the following two commands to retrieve the data necessary for the package dependencies to function properly.

```
pyensembl install --release 75 76
mhcflurry-downloads fetch
```

### Getting Amino Acid Context for Variants in a MAF

```python
from neop.main import annotate_maf
output_maf = annotate_maf(path_to_MAF_file, path_to_output_folder, PREFIX)
```

In Python, use the command above to annotate a MAF file with the amino acid context for every variant and its corresponding wild-type reference. Amino acid context is the sequence of amino acids between "Start_Position"-10 and "End_Position"+10 for in-frame mutations, and the sequence of amino acids between "Start_Position"-10 and the new stop codon for frame-shift mutations. Varcode is used to obtain amino acid context (https://github.com/openvax/varcode).

`path_to_MAF_file` is the full path of the MAF file to be annotated. Variants in the MAF must have values for the following fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short. Consider using GenomeNexus: https://github.com/genome-nexus/genome-nexus-annotation-pipeline to annotate an incomplete MAF file. Note: Tumor_Seq_Allele2 is considered to be the alt allele. The value of Tumor_Seq_Allele1 is ignored.

`path_to_output_folder` is the full path of the folder where results will be generated. If this path is invalid or has a value of: "", the output files (log file, annotated maf) will be generated in a folder called: "output_files" (which will automatically be created in whichever directory the Python script is located).

`PREFIX` is an optional string that will be used to prefix the output files (i.e. "progress.log" will become "PREFIX_progress.log"). If its value is: "", then the current time will be used as the prefix.

`output_maf` is the path to the newly generated, annotated MAF file.

The two files generated from the `annotate_maf` command will be PREFIX_progress.log and PREFIX_original_maf_name.maf (which is the file corresponding to `output_maf`). They are generated in `path_to_output_folder`. PREFIX_progress.log will have all warnings associated with the run. PREFIX_original_maf_name.maf is the original maf with three extra columns. They include (1) the transcript used for annotations, (2) variant amino acid context, (3) wild-type amino acid context. If a transcript is not specified in the MAF file (Transcript_ID column), the transcript used for annotations will be the one with the highest priority consequence (where the mutation has the most damaging effect). This is almost always the canonical transcript. In the rare case it isn't, it is important to note that variant amino acid context for different transcripts is identical nearly 100% of the time.

If a variant doesn't lead to a change in amino acid sequence (e.g. RNA, IGR, Intron, etc.) or cannot generate neoantigens (e.g. nonsense mutations), PREFIX_original_maf_name.maf will annotate blank fields for the transcript, variant amino acid context, and wild-type amino acid context. If there was an error obtaining amino acid context, it will be annotated as such (e.g. ERROR: NO VARIANT CONTEXT OBTAINED).

Results are generated very rapidly. For small MAFs (less than 1000 variants), it should take no more than a few seconds. For larger MAFs (100,000+ variants), it usually takes between a few minutes to a half-hour for results to be generated.

### Predicting Neoantigens for Variants in a MAF

```python
from neop.main import predict
output_json, output_csv = predict(path_to_patient_allele_reference, path_to_MAF_file, path_to_output_folder, PREFIX)
```

In Python, use the command above to predict neoantigens from an inputted MAF file. MHCflurry 2.0 (https://github.com/openvax/mhcflurry), which outperforms NetMHCpan 4.0 (1), is used for the predictions. When you run this command, ignore the tensorflow warnings.

`path_to_patient_allele_reference` is the full path of a patient allele reference. An example of this file is present in this repository https://github.com/as1000/neop/blob/main/tests/allele_ref.csv). Briefly, the file can have any number of columns. However, the first column must be sample IDs (corresponding to "Tumor_Sample_Barcode" in the MAF). When this file is being read, the first row (containing column names) is always skipped. However, starting from the second row, each sample ID is matched with its corresponding alleles (present in columns 2 through N), which are later used during neoantigen prediction. If the path for this file is invalid or the value of `path_to_patient_allele_reference` is set to: "", predictions are performed on 318 common HLA class I alleles covering ~98% of the population (https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1.full). Therefore, set the value of `path_to_patient_allele_reference` equal to: "", if the patient alleles are unknown. However, if WGS/WES .bam files are accessible, it is recommended that you identify patient HLA alleles using algorithms such as POLYSOLVER: https://software.broadinstitute.org/cancer/cga/polysolver_run. Alternatively expressed HLA alleles (suffixed with N, Q, etc.) are not supported. Thus, patients with these alleles will be skipped and a warning will be generated. Beyond HLA class I alleles, MHCflurry also supports HLA class II alleles, non-classical HLA alleles, and some HLA orthologs in primates, mice, cattle, sheep, pigs, salmon, and trout.

`path_to_MAF_file` is the full path of the MAF file to be annotated. Variants in the MAF must have values for the following fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short. Consider using GenomeNexus: https://github.com/genome-nexus/genome-nexus-annotation-pipeline to annotate an incomplete MAF file. Note: Tumor_Seq_Allele2 is considered to be the alt allele. The value of Tumor_Seq_Allele1 is ignored.

`path_to_output_folder` is the full path of the folder where results will be generated. If this path is invalid or has a value of: "", the output files (log file, JSON-formatted results, CSV-formatted results) will be generated in a folder called: "output_files" (which will automatically be created in whichever directory the Python script is located).

`PREFIX` is an optional string that will be used to prefix the output files (i.e. "progress.log" will become "PREFIX_progress.log"). If its value is: "", then the current time will be used as the prefix.

`output_json` is the path to the newly generated output JSON file.

`output_csv` is the path to the newly generated output CSV file.

The three files generated from the `annotate_maf` command will be PREFIX_progress.log, PREFIX_output.json (corresponding to `output_json`), and PREFIX_output.csv (corresponding to `output_csv`). They are generated in `path_to_output_folder`. PREFIX_progress.log will have the warnings associated with the run and the progress of the predictions. PREFIX_output.csv will have several columns. The first few are used to identify the variant and are identical to the original MAF file (Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2). The next fields include the transcript analyzed (Transcript), MHCflurry presentation score-based harmonic-mean best rank (Presentation_HBR), affinity-based harmonic-mean best rank (Affinity_HBR), wild-type amino acid context (WT_Amino_Acid_Context), variant amino acid context (Variant_Amino_Acid_Context), the HLA allele used in predictions (HLA_Allele), the number of variant-specific peptides binding with affinity of less than 500nM to this HLA allele (N_total_neoantigens_500nM_cutoff), the number of variant-specific peptides binding with affinity of less than 50nM to this HLA allele (N_strong_binders_50nM_cutoff), the variant-specific peptide with the best presentation score (Best_Presentation_Score_Peptide), the presentation score of that peptide (Best_Presentation_Score_Presentation_Score), the presentation score percentile of that peptide (Best_Presentation_Score_Presentation_Percentile), the affinity of that peptide (Best_Presentation_Score_Affinity), the affinity percentile of that peptide (Best_Presentation_Score_Affinity_Percentile), the variant-specific peptide with best affinity (Best_Affinity_Peptide), the presentation score of that peptide (Best_Affinity_Presentation_Score), the presentation score percentile of that peptide (Best_Affinity_Presentation_Percentile), the affinity of that peptide (Best_Affinity_Affinity), the affinity percentile of that peptide (Best_Affinity_Affinity_Percentile). PREFIX_output.json will have results from the CSV in JSON format.

The immunogenicity metric for a variant is its Harmonic-mean Best Rank (HBR). It is calculated as follows if N alleles are tested:

HBR = N/(1/R<sub>1</sub> + 1/R<sub>2</sub> ... + 1/R<sub>N</sub>)

R<sub>M</sub> is the best percentile rank for any peptide derived from the variant, on allele M.

MHCflurry presentation score-based HBR uses presentation percentile as the rank metric. Affinity-based HBR uses affinity percentile as the rank metric. The lower the HBR, the more likely a variant will give rise to a neoantigen.

This command takes about 0.5 seconds per variant per allele. For small MAFs (less than 1000 variants), results should be generated in less than one hour. For larger MAFs (100,000+ variants), it usually takes several hours or a few days for results to completely generate. However, output files will update every time predictions for five more variants have been completed.

(1) O’Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: Improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. Cell systems. 2020 Jul 22;11(1):42-8.
