#Shebang line
#! /usr/bin/bash

# Load in necessary tool kits.
module load fastqc
module load cutadapt
module load bowtie2
module load samtools
module load bedtools

# Navigate to directory with reference file.
cd /home/FCAM/mcb5430/usr13/midterm_chip_M

# Create a script file for fastqc loop with nano. Save and exit when done.
# nano hg38_loop.sh

#################################################################################################################
	# BASH SCRIPT
#################################################################################################################

# Human Bowtie2 reference genome.
index="/home/FCAM/mcb5430/STUFF/genomes/hg38/hg38_base_bt2/hg38_base"

# Human chromosome size position reference.
genome="/home/FCAM/mcb5430/STUFF/genomes/hg38/hg38_chromInfo_base.txt"

for file in fastq/*fastq.gz
	do
		echo "Starting processing of $file"

		# Run fastqc on original reads.
		echo "Running fastqc on original reads for $file"
		fastqc "$file"

		# Quality trim reads and set minimum length.
		echo "Quality trimming reads at Q30 and minimum length of 30 for $file"
		cutadapt -q 30 -m 30 -o "${file%.fastq.gz}_qtrim.fastq" "$file"

		# Remove adapter sequence from quality-trimmed reads.
		echo "Removing adapter sequence from quality-trimmed reads for $file"
		cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 30 -o "${file%.fastq.gz}_qtrim_clipped.fastq" "${file%.fastq.gz}_qtrim.fastq"

		# Run fastqc on processed reads.
		echo "Running FastQC on quality-trimmed and adapter-clipped reads for $file"
		fastqc "${file%.fastq.gz}_qtrim_clipped.fastq"

		# Bowtie2 alignment: original
		echo "Mapping original reads to the human genome for $file"
		bowtie2 -p 16 -x "$index" -U "$file" -S "${file%.fastq.gz}_raw.sam"

		# Bowtie2 alignment: combined
		echo "Mapping quality-trimmed and adapter-clipped reads to the human genome for $file"
		bowtie2 -p 16 -x "$index" -U "${file%.fastq.gz}_qtrim_clipped.fastq" -S "${file%.fastq.gz}_qtrim_clipped.sam"




# Set sample prefix.
prefix="${file%.fastq.gz}"

# Assign colors to ChIP and Input data.
	if [[ "$(basename "$prefix")" == chip* ]]; then
		color="255,0,0"
	else
		color="0,0,255"
	fi

# Convert original and processed output files.
for stage in raw qtrim_clipped
	do

		# Convert and sort SAM to BAM files.
		echo "Converting and sorting SAM files to BAM files for ${prefix}_${stage}"
		samtools sort -o "${prefix}_${stage}_sorted.bam" "${prefix}_${stage}.sam"

		# Convert BAM to BED files.
		echo "Converting BAM files to BED files for ${prefix}_${stage}"
		bedtools bamtobed -i "${prefix}_${stage}_sorted.bam" > "${prefix}_${stage}.bed"

		# Sort BED files by chromosome and starting position.
		echo "Sorting BED file ${prefix}_${stage}"
		sort -k1,1 -k2,2n "${prefix}_${stage}.bed" > "${prefix}_${stage}_sorted.bed"

		# Create bedGraph.
		echo "Creating bedGraph for ${prefix}_${stage}"
		bedtools genomecov -bg -i "${prefix}_${stage}_sorted.bed" -g "$genome" > "${prefix}_${stage}_coverage.tmp.bedGraph"

		# Add genome browser track line and sample color.
		echo "Formatting genome browser track for ${prefix}_${stage}"
		{
			printf 'track type=bedGraph name="%s_%s" visibility=full color=%s\n' "$(basename "$prefix")" "$stage" "$color"
			cat "${prefix}_${stage}_coverage.tmp.bedGraph"
		} > "${prefix}_${stage}.bedGraph"

		# Remove temporary coverage file.
		rm "${prefix}_${stage}_coverage.tmp.bedGraph"

	done

		echo "Finished processing $file"

	done 2>&1 | tee log.txt

#################################################################################################################
	# PART 2: Calling peaks with MACS2
#################################################################################################################

# Make a directory to output peak_calls to.
mkdir -p peak_calls

# Load MACS.
module load macs2/2.2.9.1

# Call peaks using macs2.
macs2 callpeak -t fastq/chipM_qtrim_clipped_sorted.bam -c fastq/inputM_qtrim_clipped_sorted.bam -f BAM -g hs -n chipM --outdir peak_calls 2>&1 | tee -a peak_calls/chipM_peaks_log.txt

# Cut out the first five columns from chipM_peaks.narrowPeak.
cut -f 1-5 peak_calls/chipM_peaks.narrowPeak > peak_calls/chipM_peaks.bed

# Add track headers to the peaks file and summits file for display in the Genome Browser.
awk 'BEGIN {print "track name=chipM_peaks description=\"ChIP Peaks\" visibility=pack useScore=0"} {print}' \
peak_calls/chipM_peaks.bed > peak_calls/chipM_peaks_track.bed

awk 'BEGIN {print "track name=chipM_summits description=\"ChIP Summits\" visibility=pack useScore=0"} {print}' \
peak_calls/chipM_summits.bed > peak_calls/chipM_summits_track.bed

# Sort the narrowPeak file to display the largest and most significant peak.
sort -k9,9nr peak_calls/chipM_peaks.narrowPeak | head -1

#################################################################################################################
	# PART 3
#################################################################################################################

# Make a directory. 
mkdir -p motif

# Sort summit file by motifs.
sort -k5,5nr peaks/summits_M.bed | head -500 > motif/top500_summits.bed

# Create 100 bp windows surrounding each summit.
bedtools slop -i motif/top500_summits.bed -g "$genome" -l 50 -r 49 > motif/top500_windows.bed

# Set FASTA reference file.
FASTA="/home/FCAM/mcb5430/STUFF/genomes/hg38/fasta_hg38/hg38.fa"

# Extract the top 500 DNA sequences using the FASTA reference.
bedtools getfasta -fi "$FASTA" -bed motif/top500_windows.bed -fo motif/top500_summits.fasta -name

# Load MEME.
module load meme

# Analyze the output summits file for motifs and output 2 motifs.
meme motif/top500_summits.fasta -dna -revcomp -mod zoops -nmotifs 2 -minw 6 -maxw 12 -oc motif/meme_out

# Extract all peaks in preparation for MAST.
bedtools getfasta -fi "$FASTA" -bed peaks/peaks_M.bed -fo motif/all_peaks.fasta -name

# Run MAST on all peaks to find identified motif occurrence.
mast motif/meme_out/meme.txt motif/all_peaks.fasta -oc motif/mast_out

# Create a hit_list for locations that a motif-of-interest was found.
mast -hit_list motif/meme_out/meme.txt motif/all_peaks.fasta > motif/mast_hits.txt

# Count peaks with motifs and calculate percentage.
total=$(grep -c '^>' motif/all_peaks.fasta)
found=$(awk 'NF && $1 !~ /^#/ {print $1}' motif/mast_hits.txt | sort -u | wc -l)
awk -v found="$found" -v total="$total" 'BEGIN {printf "Peaks with motifs: %d / %d (%.2f%%)\n", found, total, 100*found/total}'

# Set chromosome 3=13 reference.
CHR13_FASTA="/home/FCAM/mcb5430/STUFF/genomes/hg38/fasta_hg38/ind_chr_hg38_fasta/chr13.fa"

# Set background model for FIMO.
background="/home/FCAM/mcb5430/STUFF/genomes/hg38/hg38_markov_bkgrnd.txt"

# Run FIMO
fimo --bgfile "$background" --oc motif/fimo_chr13 motif/meme_out/meme.txt "$CHR13_FASTA"


# Convert FIMO to BED.
awk 'BEGIN {FS=OFS="\t"} !/^#/ && $4 ~ /^[0-9]+$/ {print $3, $4-1, $5, $1}' motif/fimo_chr13/fimo.tsv > motif/fimo_chr13.bed


# Extract chromosome 13 peaks from BED file.
awk '$1=="chr13"' peaks/peaks_M.bed > motif/chipM_chr13_peaks.bed


# Count motif occurrences overlapping peaks.
inside=$(bedtools intersect -a motif/fimo_chr13.bed -b motif/chipM_chr13_peaks.bed -u | wc -l)


# Count motif occurrences outside peaks.
outside=$(bedtools intersect -a motif/fimo_chr13.bed -b motif/chipM_chr13_peaks.bed -v | wc -l)


# Display your results.
echo "Motifs inside peaks: $inside"
echo "Motifs outside peaks: $outside"
echo "Total motif occurrences: $((inside + outside))"

# Identify transcription factors binding motifs.
tomtom -oc motif/tomtom_out motif/meme_out/meme.txt /home/FCAM/mcb5430/STUFF/TFdb/JASPAR2018_CORE_vertebrates_redundant.meme
