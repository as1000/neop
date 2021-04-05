# Neop
## Python Package for Predicting Neoantigens and Obtaining Amino Acid Context from MAFs

You can install the package using pip. Download this repo, navigate to the root directory (containing setup.py) and enter the following command:

`pip install .`

Once this package is installed, you are ready to begin predicitions.

In Python, use the following command to obtain the amino acid context for each variant and its corresponding reference allele in a MAF file (i.e. the sequence of amino acids between "Start_Position"-10 and "End_Position"+10 for in-frame mutations, and the sequence of amino acids between "Start_Position"-10 and the new stop codon for frame-shift mutations):

```
from neop.main import annotate_maf
output_maf = annotate_maf(path_to_MAF_file, path_to_output_folder, PREFIX)
```python

`path_to_MAF_file` is the full path of the MAF file to be annotated. Note: Tumor_Seq_Allele2 is considered as the alt allele. Tumor_Seq_Allele1 is ignored. Variants in the MAF must have values for the following fields: Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, Variant_Classification, Tumor_Sample_Barcode, HGVSp_Short. Consider using GenomeNexus: https://github.com/genome-nexus/genome-nexus-annotation-pipeline to annotate an incomplete MAF file.

`path_to_output_folder` is the full path of the folder where results will be generated. If this path is invalid or has a value of: "", the output files (log file, annotated maf) will be generated in a folder called: "output_files" (which will automatically be created in whichever directory the Python script is located).
`PREFIX` is an optional string that will be used to prefix the output files (i.e. "progress.log" will become "PREFIX_progress.log"). If its value is: "", then the current time will be used as the prefix.

`output_maf` is the path to the annotated maf file.

The two files generated from `annotate_maf` will be PREFIX_progress.log and PREFIX_original_maf_name.maf (corresponding to `output_maf`). They are generated in `path_to_output_folder`. PREFIX_progress.log will have warnings associated with the run. PREFIX_original_maf_name.maf is the original maf with three extra columns. They include (1) variant amino acid context, (2) reference allele amino acid context, and (3) transcript used for annotations. Unless transcripts are specified in the MAF file, the transcript used for annotations will be the one with the highest priority consequence (where the mutation has the most damaging effect). This is almost always the canonical transcript. In the rare case it isn't, it is important to note that variant amino acid context for different transcripts is identical nearly 100% of the time.

If a variant doesn't lead to a change in amino acid sequence (e.g. RNA, IGR, Intron, etc.) or cannot generate neoantigens (e.g. nonsense mutations), PREFIX_original_maf_name.maf will annotate blank fields for the variant amino acid context, reference allele amino acid context, and transcript used for annotations. If there was an error obtaining amino acid context, it will be annotated as such (e.g. ERROR: NO VARIANT CONTEXT OBTAINED).

This command is generally very quick. For small MAFs (less than 1000 variants), results should be generated in less than a few seconds. For larger MAFs (100,000+ variants), it usually takes between a few minutes to a half-hour for results to be generated.

The second command associated with this package can be used to predict neoantigens if a MAF file is inputted:

`from neop.main import predict`

`output_json, output_csv = predict(path_to_patient_allele_reference, path_to_MAF_file, path_to_output_folder, PREFIX)`

When you run this command, ignore the tensorflow warnings.

`path_to_patient_allele_reference` is the full path of a patient allele reference. An example of this file is present in this repo (http://www.github.com/as1000/neop/tests/allele_ref.csv). Briefly, the file can have any number of columns. However, the first column must be sample IDs (corresponding to "Tumor_Sample_Barcode" in the MAF). The first row (column names) is skipped. However, starting from the second row, each sample ID is matched with its corresponding alleles (present in columns 2 through N). This matching is used later on during the neoantigen predictions. If the path is invalid or the value of the variable is set to: "", predictions are performed on 318 common HLA class I alleles covering ~98% of the population (https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1.full). Set the value of `path_to_patient_allele_reference` equal to: "", if the patient alleles are unknown. However, it is recommended that HLA alleles are identified through algorithms such as POLYSOLVER: https://software.broadinstitute.org/cancer/cga/polysolver_run if WGS/WES BAM files are accessible. Alternatively expressed HLA alleles (suffixed with N, Q, etc.) are not supported. Thus, patients with these alleles will be skipped and a warning will be generated. Beyond HLA class I alleles, there is also support for HLA class II alleles, non-classical HLA alleles, and some HLA orthologs in primates, mice, cattle, sheep, pigs, salmon, and trout.
`path_to_MAF_file` is the full path of the MAF file to be annotated. Note: Tumor_Seq_Allele2 is considered as the alt allele. Tumor_Seq_Allele1 is ignored.
`path_to_output_folder` is the full path of the folder where results will be generated. If this path is invalid or has a value of: "", the output files (log file, JSON-formatted results, CSV-formatted results) will be generated in a folder called: "output_files" (which will automatically be created in whichever directory the Python script is located).
`PREFIX` is an optional string that will be used to prefix the output files (i.e. "progress.log" will become "PREFIX_progress.log"). If its value is: "", then the current time will be used as the prefix.

`output_json` is the path to the output JSON file.
`output_csv` is the path to the output CSV file.

The three files generated from `annotate_maf` will be PREFIX_progress.log, PREFIX_output.json (corresponding to `output_json`), and PREFIX_output.csv (corresponding to `output_csv`). They are generated in `path_to_output_folder`. MHCflurry is used for the predictions. PREFIX_progress.log will have warnings associated with the run and progress on the predictions. PREFIX_output.json will have results from the CSV in JSON format. PREFIX_output.csv will have several columns. The first few are used to identify the variant and are identical to the original MAF file (Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2). The next fields include the transcript analyzed (Transcript), MHCflurry presentation score-based harmonic mean best rank (Presentation_HBR), affinity-based harmonic mean best rank (Affinity_HBR), reference amino acid context (WT_Amino_Acid_Context), variant amino acid context (Variant_Amino_Acid_Context), the HLA allele used in predictions (HLA_Allele), the number of peptides derived from the variant binding with affinity of less than 500nM to this HLA allele (N_total_neoantigens_500nM_cutoff), the number of peptides derived from the variant binding with affinity of less than 50nM to this HLA allele (N_strong_binders_50nM_cutoff), the peptide with best presentation score (Best_Presentation_Score_Peptide), the presentation score of that peptide (Best_Presentation_Score_Presentation_Score), the presentation score percentile of that peptide (Best_Presentation_Score_Presentation_Percentile), the affinity of that peptide (Best_Presentation_Score_Affinity), the affinity percentile of that peptide (Best_Presentation_Score_Affinity_Percentile), the peptide with best affinity (Best_Affinity_Peptide), the presentation score of that peptide (Best_Affinity_Presentation_Score), the presentation score percentile of that peptide (Best_Affinity_Presentation_Percentile), the affinity of that peptide (Best_Affinity_Affinity), the affinity percentile of that peptide (Best_Affinity_Affinity_Percentile).

The immunogenicity metric for a variant is the Harmonic-mean Best Rank (HBR). It is calculated as follows if N alleles are tested:

HBR = N/(1/R<sub>1</sub> + 1/R<sub>2</sub> ... + 1/R<sub>N</sub>)

R<sub>M</sub> is the best percentile rank for any peptide derived from the variant, on allele M.

The lower the HBR, the more likely a variant will give rise to a neoantigen.

This command is generally takes more time to run. For small MAFs (less than 1000 variants), results should be generated in less than one hour. For larger MAFs (100,000+ variants), it usually takes a few days for results to completely generate. A good rule of thumb is about 0.5 seconds per variant per allele.
