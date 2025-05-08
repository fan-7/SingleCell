########################################
# 1. 环境设置与变量定义模块（Shell）
#######################################
#!/bin/bash
# 定义环境变量，激活相关环境（假设atac环境已配置好相关软件）
export PATH=/software/biosoft/software/python/python2.7/bin:$PATH
source activate atac

# 定义关键路径和变量
rawdata_dir="your_rawdata_path"
ref_dir="your_ref_path"
software_dir="your_software_path"
result_dir="your_result_path"
type="fastq"
read="R1"
dataname="sgg-mESCWTKO19"
filename_list=("Input-2" "IP-2")
thres=6

########################################
# 2. 下机数据量和碱基质量质控模块（Shell）
#######################################
# 遍历样本组
for group in ${filename_list[@]}; do
    # 在结果目录下创建样本组文件夹
    cd ${result_dir}
    mkdir ${group}
    cd ${group}
    mkdir qc_raw
    # 复制linker.list文件到质控文件夹
    cp $ref_dir/linker.list ${result_dir}/${group}/qc_raw
    # 查找原始数据文件并记录文件名
    find $rawdata_dir -name ${dataname}-${group}*.${type}.gz >${result_dir}/${group}/${group}_rawdata.filename
    # 对每个原始数据文件进行fastqc质控
    less ${result_dir}/${group}/${group}_rawdata.filename | while read id; do 
        (fastqc $id -t $thres -a qc_raw/linker.list -o qc_raw/. );
    done
    # 汇总fastqc结果
    cd qc_raw
    multiqc *_fastqc.zip
done
########################################
#3. 单细胞连接效率质控模块（Shell）
#######################################
for group in ${filename_list[@]}; do
    # 创建tag文件夹
    cd $result_dir
    mkdir $group
    cd $group
    mkdir tag_raw
    cd tag_raw
    # 检查原始文件并记录
    find $rawdata_dir -name ${dataname}-${group}*.${type}.gz > ${result_dir}/${group}/${group}_rawdata.filename
    # 备份原始文件到结果文件夹
    find $rawdata_dir -name ${dataname}-${group}*.${type}.gz -type f|xargs -i scp {} ${result_dir}/${group}/tag_raw
    # 执行相关脚本处理数据
    cd $result_dir/$group/tag_raw
    find ${result_dir}/${group}/tag_raw -name \*${read}\* > input.txt
    cp ${software_dir}/tag-raw-fastq.sh ./
    cp ${ref_dir}/SGG-SP-barcode.txt ./
    sh tag-raw-fastq.sh input.txt SGG-SP-barcode.txt
done
########################################
#4. 单细胞标记成功的序列数量和长度分布质控模块（Shell）
########################################
for group in ${filename_list[@]}; do
    # 创建质控文件夹
    cd ${result_dir}/${group}
    mkdir qc_tag
    cd qc_tag
    # 合并R1、R2标记序列
    find ${result_dir}/${group} -name \*R1_N4uniq\* | xargs cat | gzip > R1_pool.fastq.gz
    find ${result_dir}/${group} -name \*R2_N4uniq\* | xargs cat | gzip > R2_pool.fastq.gz
    # 统计序列数量和长度
    gunzip -c R1_pool.fastq.gz | awk 'NR%4==2' | wc -l > readsCount
    gunzip -c R1_pool.fastq.gz | awk 'NR%4==2' | awk '{print length($1)}' | sort -n > length.raw
    # 绘制长度分布密度图（假设脚本已存在）
    cp ${software_dir}/plot_raw_reads_length_density.R ./
    Rscript plot_raw_reads_length_density.R
    # 去除接头和polyA等，保留至少15个碱基的reads
    cutadapt -m 15 -e 0.1 --trim-n -a GATCGGAAGA -a "polyG1=GG{5}" -a "polyC1=CC{5}" -a "polyT1=TT{5}" -a "polyA1=AA{5}" -o cut_R1_pool.fastq.gz  R1_pool.fastq.gz
    # 再次进行fastqc质控
    fastqc *.fastq.gz -t 6 -a ${result_dir}/${group}/qc_raw/linker.list -o ./
    # 统计处理后序列数量和长度并绘图
    gunzip -c cut_R1_pool.fastq.gz | awk 'NR%4==2' | wc -l >> readsCount
    cp ${software_dir}/plot_cut_reads_length_density.R ./
    gunzip -c cut_R1_pool.fastq.gz  | awk 'NR%4==2' | awk '{print length($1)}' | sort -n > length
    Rscript plot_cut_reads_length_density.R
done
########################################
#5. 核糖体扣除效率质控模块（Shell）
########################################
# 定义相关变量
ribo_fa="your_ribo_fa_path"
dataname="sgg-mESCWTKO19"
riboselect_dir="your_riboselect_script_path"
bamCoverage="/software/biosoft/software/python/python2.7_2018_12/bin/bamCoverage"
filename_list=("Input-2" "IP-2")
thres=6

for group in ${filename_list[@]}; do
    filename=${group}
    workdir=${result_dir}/${group}/qc_tag
    out="bwa"
    cd $workdir
    mkdir ribo
    cd  ${workdir}/ribo
    # 短序列比对到核糖体序列
    bwa aln  ${ribo_fa} ../cut_R1_pool.fastq.gz  > aln_fq.sai
    bwa samse ${ribo_fa} ../cut_R1_pool.fastq.gz | samtools view -Sb > ${out}.aln-ribo.bam
    # 常规比对
    bwa mem -t 6 -h 15 ${ribo_fa} ../cut_R1_pool.fastq.gz | samtools view -Sb > ${out}.mem-ribo.bam
    # 合并比对结果
    samtools merge -n -r -h ${out}.aln-ribo.bam --threads 6 ${out}.all-ribo.bam ${out}.aln-ribo.bam ${out}.mem-ribo.bam
    rm ${out}.aln-ribo.bam ${out}.mem-ribo.bam aln_fq.sai
    # 统计比对到核糖体的reads数量
    echo "ribo reads with no -q:" > readsCount
    samtools view ${out}.all-ribo.bam |grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
    # 生成bw文件用于IGV查看
    mkdir bw
    samtools sort ${out}.all-ribo.bam  -@ 6 > bw/sort_${filename}.bam
    samtools index bw/sort_${filename}.bam
    cd bw
    $bamCoverage -p 6 --normalizeUsing RPKM -bs 1 -b sort_${filename}.bam -o ${filename}-rpkm.bw
    $bamCoverage -p 6 --normalizeUsing None -bs 1 -b sort_${filename}.bam -o ${filename}-count.bw
    # 提取非核糖体reads
    cd  ${workdir}/ribo
    stranded="n"
    samtools sort -n --threads 6 ${out}.all-ribo.bam -O BAM -o ${out}.nsorted.all-ribo.bam
    rm ${out}.all-ribo.bam
    python ${riboselect_dir} ${out}.nsorted.all-ribo.bam $stranded ${out}
done
########################################
#6. 参考基因组比对模块（Shell）
########################################
# 定义相关变量
STAR_dir="/software/biosoft/software/STAR-2.7.7a/source/STAR"
ref_dir="your_ref_dir"
ref="your_ref_file"

for group in ${filename_list[@]}; do
    workdir=${result_dir}/${group}/qc_tag
    cd $workdir
    mkdir star
    cd star
    # 使用STAR比对到参考基因组
    ${STAR_dir} --runThreadN 6 --genomeDir ${ref_dir}/star_mm10 --readFilesIn $workdir/ribo/bwa.nonRibo.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 20 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes All  --outFileNamePrefix star_nonRibo
    cd ${workdir}/star
    mkdir count
    cd count
    # 统计比对的reads数量
    echo "reads with no -q:" > readsCount
    samtools view ${workdir}/star/star_nonRiboAligned.out.bam |grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
    # 过滤质量低于20的reads
    samtools view ${workdir}/star/star_nonRiboAligned.out.bam -q 20 -h > star2_20.sam
    echo "reads with -q:" >> readsCount
    cat star2_20.sam|grep -v "^@"|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
    # 提取唯一比对的reads
    cat star2_20.sam|grep -E "^@|NH:i:1$|NH:i:1[^0-9]" > uniqmap.sam
    echo "reads with uniqmap:" >> readsCount
    cat uniqmap.sam|cut -f1|awk '!x[$0]++'|wc -l >> readsCount
    # 转换文件格式为bam和bed
    samtools view -S uniqmap.sam -b -o uniqmap.bam
    bedtools bamtobed -i uniqmap.bam -split > uniqmap.bed
    rm uniqmap.sam
    # 注释转录本
    intersectBed -a uniqmap.bed -b ${ref} -wa -wb -f 0.51 > tmp.SG2
    echo "SG2 mapping reads:" >> readsCount
    awk '!x[$4]++' uniqmap.bed|wc >> readsCount
    # 统计RNA数量
    echo "SG2 ncRNA+mRNA reads:" >> readsCount
    awk '!a[$4]++' tmp.SG2|wc >> readsCount
    # 统计编码RNA区域分布
    echo "SG2 mRNA region reads:" >> readsCount
    awk '$13=="protein_coding"' tmp.SG2|cut -f14|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
    # 统计所有RNA区域分布
    echo "SG2 mRNA+ncRNA region reads:" >> readsCount
    cut -f14 tmp.SG2|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
    # 统计编码RNA类型分布
    echo "SG2 mRNA+ncRNA biotype reads:" >> readsCount
    cut -f13 tmp.SG2|sed -e 's/_/\t/'|cut -f1|awk -v OFS="\t" '{num[$1]++}END{for(i in num){print i,num[i]}}' >> readsCount
done
########################################
#7. 基因组分布质控模块（Shell）
########################################
for group in ${filename_list[@]}; do
    software_dir="your_software_path"
    bamCoverage="/software/biosoft/software/python/python2.7_2018_12/bin/bamCoverage"
    cd ${result_dir}/${group}/qc_tag/star/count/
    mkdir bw
    cd bw
    # 排序bam文件
    samtools sort ../uniqmap.bam -@ 6 > sort_uniqmap.bam
    samtools index sort_uniqmap.bam
    # 生成BW文件用于IGV查看
    $bamCoverage -p 6 --normalizeUsing RPKM -b ./sort_uniqmap.bam -bs 10 -o ${filename}-uniqmap-rpkm-bs10.bw
    $bamCoverage -p 6 --normalizeUsing None -b ./sort_uniqmap.bam -bs 10  -o ${filename}-uniqmap-count-bs10.bw
    $bamCoverage -p 6 --normalizeUsing RPKM -b ./sort_uniqmap.bam -bs 1 -o ${filename}-uniqmap-rpkm-bs1.bw
    $bamCoverage -p 6 --normalizeUsing None -b ./sort_uniqmap.bam -bs 1  -o ${filename}-uniqmap-count-bs1.bw
    # 计算基因体覆盖度
    export PATH=/software/biosoft/software/python/python2.7_2018_12/bin:$PATH
    export PYTHONPATH=/software/biosoft/software/python/python2.7_2018_12/lib/python2.7/site-packages
    export PATH=/software/biosoft/software/python/python2.7/bin:$PATH
    cd ${result_dir}/${group}/qc_tag/star/count/bw
    /software/biosoft/software/RSeQC-2.6.4/scripts/geneBody_coverage.py -r /p300s/yangyg_group/fanxiuA/project/brain-scm6A/1_ref/mm10/mm10_RefSeq.bed -i sort_uniqmap.bam -o uniqmap_genecoverage
    # 使用QoRTs进行二次质控
    QoRTs_dir="your_QoRTs_jar_path"
    chrSize="your_chrSize_path"
    ref_G="your_ref_G_path"
    cd ${result_dir}/${group}/qc_tag/star/count/bw
    mkdir qorts_alluniqmap
    java -Xmx60G -jar $QoRTs_dir QC \
        --stranded \
        --maxReadLength 150 \
        --singleEnded \
        --chromSizes $chrSize \
        sort_uniqmap.bam \
        $ref_G \
        ${result_dir}/${group}/qc_tag/star/count/bw/qorts_alluniqmap
    # 绘制count分布和coverage图（假设脚本已存在）
    cd ${result_dir}/${group}/qc_tag/star/count/
    cp $software_dir/plot_count_cover_order.R ./
    Rscript plot_count_cover_order.R
done
########################################
#8. m6A 分析质控模块（Shell）
########################################
# 定义相关变量
rawdata_dir="your_rawdata_path"
ref_dir="your_ref_dir"
software_dir="your_software_path"
result_dir="your_result_path"
macs2_dir="/software/biosoft/software/python/python2.7_2018_12/bin/macs2"
ref="your_m6A_ref_file"
thres=6
case="IP-2"
ctrl="Input-2"
cluster="wt"
filename="IP_Input_ext200"

# 调用peak
cd $result_dir/$case
mkdir peak
cd $result_dir/$case/peak
mkdir $cluster
cd $cluster
mkdir $filename
cd $filename
$macs2_dir callpeak -t  ${result_dir}/${case}/qc_tag/star/count/bw/sort_uniqmap.bam -c ${result_dir}/${ctrl}/qc_tag/star/count/bw/sort_uniqmap.bam -n $filename --nomodel --keep-dup all -q 0.05 -B --SPMR --extsize 200
wc -l IP_Input_ext200_peaks.narrowPeak
# peak分布分析
cd $result_dir/$case/peak/$cluster/${filename}
cp $software_dir/distribution_plot_Mus_m6A.pl ./
cut -f1,2,3,7 ${filename}_peaks.narrowPeak > ok.bed
perl distribution_plot_Mus_m6A.pl
# 注释peak
intersectBed -a ok.bed -b ${ref} -wa -wb -f 0.51 | awk '!seen[$1,$2,$3]++' > anno_ok.bed
file_anno_ok=anno_ok.bed
# 统计peak数量
echo "peak_with_anno_ok" > peaksCount
cat $file_anno_ok | wc -l >> peaksCount
echo "nonCodingRNA" >> peaksCount
grep -v -E 'intron|utr3|utr5|cds' $file_anno_ok | wc -l >> peaksCount
echo "utr3" >> peaksCount
grep "utr3" $file_anno_ok | wc -l >> peaksCount
echo "utr5" >> peaksCount
grep "utr5" $file_anno_ok | wc
