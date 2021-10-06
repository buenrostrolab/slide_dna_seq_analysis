
## 0. Set up environment

home_dir='/broad_thechenlab/zack/slidedna_alignment'
data_dir='/broad_thechenlab/Nextseq/200930_NB501164_1146_AHC3LVAFX2/Data/Intensities/BaseCalls/201001_demultiplex/raw_fastq/'
#data_dir='/broad_thechenlab/Novaseq/201015_SL-NVK_0388_AHT3GLDRXX/Data/Intensities/BaseCalls/1_original_demultiplex/'
run='5_cerebellumNormal4x_puck16_S5'
#run='1_colonNormal4x_puck12_S1'
#run='6_cerebellumPEG4x_puck17_S2'
genome='/broad_thechenlab/atac_pipeline/genomes/bowtie2/mm10/mm10'

step=10
num_cores=8
#split_size=4000000 # 1M reads, use for NextSeq (need more files to fully parallelize)
split_size=20000000 # 5M reads, use for NovaSeq (less files to check)

bin='/broad_thechenlab/zack/slidedna_alignment/bin'
tmp=$home_dir/tmp/

cd $home_dir
mkdir -p $run

r1_path=`ls $data_dir/$run"_"*R1*.fastq.gz`
r2_path=`echo $r1_path | sed 's/R1/R2/g'`
r3_path=`echo $r1_path | sed 's/R1/R3/g'`

r1=$run"_R1"
r2=$run"_R2"
r3=$run"_R3"

## 1. Split FASTQ

mkdir -p $tmp/split/

if  [ $step -le 1 ]; then

	echo "$run - Spliting FASTQs (1)"
	bin/fastp -i $r1_path -o $tmp/split/$r1.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L  -Q 2>logs/split_R1.log & 
	bin/fastp -i $r2_path -o $tmp/split/$r2.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L  -Q 2>logs/split_R2.log &
	bin/fastp -i $r3_path -o $tmp/split/$r3.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L  -Q 2>logs/split_R3.log &
	wait

fi

ls $tmp/split/ | grep $run | grep "R1" | grep -P -o "[0-9]{4}" > $run/split_list.txt
num_split=`tail $run/split_list.txt`

## 2. Extract barcode (.ext)

mkdir -p $tmp/ext/

if  [ $step -le 2 ]; then

	echo "$run - Extracting barcodes (2)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		python $bin/extractBCfromR2.py \
		$tmp/split/{1}.$r1.fastq.gz $tmp/split/{1}.$r2.fastq.gz $tmp/split/{1}.$r3.fastq.gz \
		$tmp/ext/{1}.$r1.ext.fastq.gz $tmp/ext/{1}.$r2.ext.fastq.gz $tmp/ext/{1}.$r3.ext.fastq.gz \
		2>logs/ext.log :::: $run/split_list.txt

fi

## 3. Move R2 barcode to header (.umi)

mkdir -p $tmp/umi/

if  [ $step -le 3 ]; then
	
	echo "$run - Moving R2 barcode to header (3)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		python $bin/add_umis_from_R2.py \
		$tmp/ext/{1}.$r1.ext.fastq.gz $tmp/ext/{1}.$r2.ext.fastq.gz $tmp/ext/{1}.$r3.ext.fastq.gz \
		$tmp/umi/{1}.$r1.umi.fastq.gz $tmp/umi/{1}.$r2.umi.fastq.gz \
		2>logs/umi.log :::: $run/split_list.txt

fi

## 4. Trim adapters (.trim)

mkdir -p $tmp/trim/

if  [ $step -le 4 ]; then

	echo "$run - Trimming adapters (4)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		$bin/pyadapter_trim.py \
		-a $tmp/umi/{1}.$r1.umi.fastq.gz -b $tmp/umi/{1}.$r2.umi.fastq.gz \
		-c $tmp/trim/{1}.$r1.trim.fastq.gz -d $tmp/trim/{1}.$r2.trim.fastq.gz \
		1>logs/trim.log :::: $run/split_list.txt
fi

## 5. Alignment (.aln)

mkdir -p $tmp/aln/

if  [ $step -le 5 ]; then

	echo "$run - Aligning (5)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		bowtie2 -X2000 -p1 --rg-id $run \
		-x $genome \
		-1 $tmp/trim/{1}.$r1.trim.fastq.gz \
		-2 $tmp/trim/{1}.$r2.trim.fastq.gz '|' \
		samtools view -bS - -o $tmp/aln/{1}.$run.aln.bam \
		2>logs/align.log :::: $run/split_list.txt

fi

## 6. Sort (.sort)

mkdir -p $tmp/sort/

if  [ $step -le 6 ]; then


	echo "$run - Sorting (6)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools sort $tmp/aln/{1}.$run.aln.bam -o $tmp/sort/{1}.$run.sort.bam \
		2>logs/sort.log :::: $run/split_list.txt
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools index $tmp/sort/{1}.$run.sort.bam \
		2>logs/index.log :::: $run/split_list.txt

fi

## 7. Filter (.flt)

mkdir -p $tmp/flt/
chrs=`samtools view -H $tmp/sort/0001.$run.sort.bam | grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}'`

if  [ $step -le 7 ]; then

	echo "$run - Filtering by quality, chromosome, pairing (7)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools view -b -q 30 -f 0x2 $tmp/sort/{1}.$run.sort.bam -o $tmp/flt/{1}.$run.flt.bam `echo $chrs` \
		2>logs/filter.log :::: $run/split_list.txt
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools index $tmp/flt/{1}.$run.flt.bam \
		2>logs/index.log :::: $run/split_list.txt

fi

## 8. Add barcode to UM tag (.tag)

mkdir -p $tmp/tag/

if  [ $step -le 8 ]; then

	echo "$run - Adding barcode to UM tag (8)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
    		python $bin/move_umis_to_UM_tag.py \
		$tmp/flt/{1}.$run.flt.bam \
		$tmp/tag/{1}.$run.tag.bam \
		2>logs/tag.log :::: $run/split_list.txt
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools index $tmp/tag/{1}.$run.tag.bam \
		2>logs/index.log :::: $run/split_list.txt

fi

## 9. Merge all (.merge)

if  [ $step -le 9 ]; then
	
	echo "$run - Merging (9)"
	ls $tmp/tag/*bam > $run/merge_list.txt	
	samtools merge -f -b $run/merge_list.txt --threads $num_cores $run/$run.merge.bam
	samtools index $run/$run.merge.bam
fi

## 10. UMI-tools (.group)

if  [ $step -le 10 ]; then
	
	echo "$run - Running UMI-tools (10)"
	umi_tools group \
		-I $run/$run.merge.bam \
		--extract-umi-method=tag --umi-tag=UM:Z \
		--method=cluster --edit-distance-threshold=2 \
		--group-out=$run/$run.group.tsv \
		--log=$run/$run.group.log \
		--paired --output-bam \
		-S $run/$run.group.bam

    	cat $run/$run.group.tsv | grep -v "contig" | cut -f2,3,7,8,9 | sort | uniq | sort -k5,5 -n | grep -P -v "N" > $run/$run.reads.txt
    	cat $run/$run.reads.txt | cut -f3 | sort | uniq -c | sort -n -r > $run/$run.barcodes.txt
fi

