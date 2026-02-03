# Fungal genome annotation pipeline
### Two steps were applied for genome annotation. 
### Step 1: gene structure prediction and RNA prediction
### Step 2: coding gene annotation (Pfam, Interproscan, eggNOG, SignalP, TMHMM, EffectorP, VFDF)

## Usage

## Step 1: gene structure annotation

```bash
bash fungi_annotation_step1.sh -genome <genome files> -threads  -database <path to database> -predictsoftware <genemark/augustus> [-species <species name>] [-structurePrograms RNAmmer,tRNAscan,busco,quast]
```

## Step 2: functional gene annotation

```bash
bash fungi_annotation_step2.sh -faa <protein sequence> -threads  -database <path to database> [-annotationPrograms emapper,signalp,Deeptfactor,SwissProt,DFVF,FungalP450s,NR,CAZy,antiSMASH,EffectorP,Interproscan,KOG,Pfam]
```

```bash
#!/bin/bash

while true; do
    case "$1" in
        -genome) genome=$2; shift 2;;
        -threads) threads=$2; shift 2;;
        -database) database=$2; shift 2;;
        -predictsoftware) predictsoftware=$2; shift 2;;
        -species) species=$2; shift 2;;
        -structurePrograms) structurePrograms=$2; shift 2;;
        *) break;;
    esac
done

# default program
structurePrograms=${structurePrograms:-""}

# make output directory
mkdir -p tmp/pre-analysis
mkdir -p tmp/annotation/repeats/repeatmasker
mkdir -p tmp/annotation/repeats/repeatmodeler
mkdir -p tmp/annotation/predictgene
mkdir -p result/Structure

########## 1. extract sample name
samplename=$(perl get_name.pl ${genome} | head -n 1)
echo "Sample name: $samplename"

########## 2. upload genome file
genome_file="tmp/${samplename}.fasta"
cp ${genome} ${genome_file}

########## 3. Funannotate clean and sort
export PATH=/software/funannotate/bin:$PATH
funannotate clean -i ${genome_file} --minlen 1000 -o tmp/pre-analysis/${samplename}.clean.fasta
funannotate sort -i tmp/pre-analysis/${samplename}.clean.fasta -b scaffold -o tmp/pre-analysis/${samplename}.clean.sorted.fa
python funanno_name_cor.py -indir tmp/ -sample ${samplename}

sortfile="tmp/pre-analysis/${samplename}.clean.sorted.fa"

########## 4. repeat masking
# RepeatMasker
RepeatMasker -pa ${threads} -html -gff -a -e ncbi -species fungi \
    -dir tmp/annotation/repeats/repeatmasker \
    -libdir ${database}/Libraries \
    ${sortfile}

# RepeatModeler
BuildDatabase -name ${samplename} ${sortfile}
RepeatModeler -threads ${threads} -database ${samplename}

# RepeatModeler结果
RepeatMasker -pa ${threads} -html -gff -a -e ncbi \
    -lib RM*/consensi.fa.classified \
    -dir tmp/annotation/repeats/repeatmodeler \
    ${sortfile}

# merge results from repeatMaker and repeatmodeler
perl merge_repeatMasker_out.pl ${sortfile} \
    tmp/annotation/repeats/repeatmasker/${samplename}.clean.sorted.fa.out \
    tmp/annotation/repeats/repeatmodeler/${samplename}.clean.sorted.fa.out \
    > ${samplename}.repeat.stats

awk 'NR>3{print $0}' genome.repeat.out > ${samplename}.repeat.out
mv genome.repeat.gff3 ${samplename}.repeat.gff

# soft masking
cp ${sortfile} ${samplename}.sort.fa
perl /opt/RepeatMasker/util/maskFile.pl -fasta ${samplename}.sort.fa \
    -annotations ${samplename}.repeat.out -softmask
mv ${samplename}.sort.fa.masked ${samplename}.combined.masked.fasta

maskfa="${samplename}.combined.masked.fasta"
cp ${maskfa} tmp/annotation/repeats/
cp ${samplename}.repeat.gff tmp/annotation/repeats/
cp ${samplename}.repeat.out tmp/annotation/repeats/

# copy to output directory
cp tmp/pre-analysis/name_cor_list.txt result/Structure/
cp ${maskfa} result/Structure/
cp ${samplename}.repeat.gff result/Structure/repeats/
cp ${samplename}.repeat.out result/Structure/repeats/

if [[ -f "$maskfa" ]]; then
    echo "Repeat masking done!" >> log
else
    echo "Repeat masking Fail!" >> log
    exit 1
fi

########## 5. tRNA/sRNA prediction（optional）
IFS=',' read -ra PROGS <<< "$structurePrograms"
for strPro in "${PROGS[@]}"; do
    strPro=$(echo $strPro | xargs)  # remove space
    
    if [ "$strPro" == "RNAmmer" ]; then
        # RNAmmer (rRNA prediction)
        mkdir -p tmp/annotation/rnammer
        rnammer -S euk -multi -m ssu,lsu \
            -xml tmp/annotation/rnammer/${samplename}.rna.xml \
            -f tmp/annotation/rnammer/${samplename}.rna.fa \
            -h tmp/annotation/rnammer/${samplename}.rna.report \
            -gff tmp/annotation/rnammer/${samplename}.rna.gff \
            ${maskfa}
        
        if [[ -f "tmp/annotation/rnammer/${samplename}.rna.gff" ]]; then
            echo "RNAmmer done!" >> log
        else
            echo "RNAmmer Fail!" >> log
        fi
    fi
    
    if [ "$strPro" == "tRNAscan" ]; then
        # tRNAscan (tRNA prediction)
        mkdir -p tmp/annotation/trnascan
        tRNAscan-SE -o tmp/annotation/trnascan/${samplename}.trna.txt ${maskfa}
        
        if [[ -f "tmp/annotation/trnascan/${samplename}.trna.txt" ]]; then
            echo "tRNAscan done!" >> log
        else
            echo "tRNAscan Fail!" >> log
        fi
    fi
    
    if [ "$strPro" == "busco" ]; then
        # BUSCO completeness evaluation
        export PATH=/software/busco_v5.6.1/bin:$PATH
        busco -i ${genome_file} -o ${samplename}_busco \
            -l ${database}/busco_v5.6.0/fungi_odb10 \
            -m genome -f -c ${threads} --offline
        
        buscofile=$(ls ${samplename}_busco/short_summary.specific.*.busco.txt)
        perl busco.pl ${buscofile} ${samplename} > result/Structure/busco.txt
        
        if [[ -f "$buscofile" ]]; then
            echo "BUSCO done!" >> log
        else
            echo "BUSCO Fail!" >> log
        fi
    fi
    
    if [ "$strPro" == "quast" ]; then
        # QUAST评估
        quast.py -o ${samplename}_quast -e --fungus -t ${threads} --glimmer ${genome_file}
        cp ${samplename}_quast/report.txt result/Structure/quast.txt
        
        if [[ -f "${samplename}_quast/report.txt" ]]; then
            echo "QUAST done!" >> log
        else
            echo "QUAST Fail!" >> log
        fi
    fi
done


########## 6. gene prediction
mkdir -p tmp/annotation/predictgene/
cd tmp/annotation/predictgene/

if [ "$predictsoftware" == "genemark" ]; then
    # GeneMark-ES
    export PATH=/root/gmes_linux_64_4:$PATH
    cp /path/to/gm_key_64 ~/.gm_key
    gmes_petap.pl --ES --fungus --cores ${threads} --sequence ../../${maskfa}
    
    cut -f 1 genemark.gtf | cut -d ' ' -f 1 > part1
    cut -f 2-9 genemark.gtf > part2
    paste part1 part2 > genemark.st.gtf
    get_sequence_from_GTF.pl genemark.st.gtf ../../${maskfa}
    rm -rf part1 part2
    mv genemark.st.gtf ${samplename}.gff
    mv prot_seq.faa ${samplename}.faa
    mv nuc_seq.fna ${samplename}.fna
    
elif [ "$predictsoftware" == "augustus" ]; then
    # Augustus
    if [ -z "$species" ]; then
        echo "Error: Species name required for Augustus"
        exit 1
    fi
    
    perl get_spe.pl "${species}" > log
    spe=$(cat log | awk '{print $0}')
    augustus --species=$spe --strand=both --genemodel=partial \
        --codingseq=on --protein=on ../../${maskfa} > ${samplename}.gff
    
    python augustus_result.py -n ${samplename} -o ./
    mv ${samplename}.prot.fa ${samplename}.faa
fi


cd ../../

faafile="tmp/annotation/predictgene/${samplename}.faa"
gfffile="tmp/annotation/predictgene/${samplename}.gff"
genefile="tmp/annotation/predictgene/Gene_count.txt"

if [[ -f "$faafile" ]] && [[ -f "$gfffile" ]]; then
    echo "Gene prediction done!" >> log
    
    # copy to output directory
    mkdir -p result/Structure
    cp ${faafile} result/Structure/
    cp ${gfffile} result/Structure/
    cp ${genefile} result/Structure/
else
    echo "Gene prediction Fail!" >> log
    exit 1
fi

echo "Step 2 completed successfully!"
```
###Step 2 gene annotation
```bash
#!/bin/bash

# 参数解析
while true; do
    case "$1" in
        -faa) faa=$2; shift 2;;
        -threads) threads=$2; shift 2;;
        -database) database=$2; shift 2;;
        -annotationPrograms) annotationPrograms=$2; shift 2;;
        *) break;;
    esac
done

# 默认值
annotationPrograms=${annotationPrograms:-""}

# 创建输出目录
mkdir -p tmp/annotation
mkdir -p result/Annotation

########## 1. 获取样本名称
samplename=$(perl get_name.pl ${faa} | head -n 1)
echo "Sample name: $samplename"

########## 2. 复制FAA文件
mkdir -p tmp/
cp ${faa} tmp/${samplename}.faa
faafile="tmp/${samplename}.faa"

########## 3. 功能注释程序
IFS=',' read -ra ANNO_PROGS <<< "$annotationPrograms"
for anno in "${ANNO_PROGS[@]}"; do
    anno=$(echo $anno | xargs)  # 去除空格
    
    if [ "$anno" == "emapper" ]; then
        # eggNOG-mapper
        mkdir -p tmp/annotation/emapper
        emapper.py -i ${faafile} \
            --output tmp/annotation/emapper/${samplename} \
            -m diamond -d euk --cpu ${threads} \
            --data_dir ${database}/eggnog_data
        
        if [[ -f "tmp/annotation/emapper/${samplename}.emapper.annotations" ]]; then
            mkdir -p result/Annotation/emapper
            grep -v "^##" tmp/annotation/emapper/${samplename}.emapper.annotations \
                | sed "s/^\#//g" > result/Annotation/emapper/emapper.txt
            echo "emapper done!" >> log
        else
            echo "emapper Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "signalp" ]; then
        # TMHMM + SignalP
        # 先运行TMHMM
        mkdir -p tmp/annotation/tmhmm
        cp ${faafile} ${samplename}.faa
        sed -i "s#\*##g" ${samplename}.faa
        cat ${samplename}.faa | decodeanhmm.Linux_x86_64 \
            -f /tmhmm-2.0c/lib/TMHMM2.0.options \
            -modelfile /tmhmm-2.0c/lib/TMHMM2.0.model > tmhmm.result.txt
        
        # 提取有跨膜结构域的序列
        perl get_tmhmm.pl tmhmm.result.txt > ${samplename}.tmhmm.list
        seqkit grep -f ${samplename}.tmhmm.list ${faafile} > tmhmm.faa
        
        # 运行SignalP
        mkdir -p tmp/annotation/signalp
        export PERL5LIB=/software/signalp-4.1/lib$PERL5LIB
        signalp -t euk -T tmp/annotation/signalp -f long tmhmm.faa \
            > tmp/annotation/signalp/${samplename}.signalp_short_out
        
        if [[ -f "tmp/annotation/signalp/${samplename}.signalp_short_out" ]]; then
            mkdir -p result/Annotation/signalp
            cp tmp/annotation/signalp/${samplename}.signalp_short_out \
                result/Annotation/signalp/signalp.txt
            echo "SignalP done!" >> log
        else
            echo "SignalP Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "Deeptfactor" ]; then
        # DeepTFactor (转录因子预测)
        mkdir -p tmp/annotation/tf
        cp ${faafile} ${samplename}.faa
        sed -i "s#\*##g" ${samplename}.faa
        cp -rf /path/to/deeptfactor ./
        export PATH=/deeptfactor/bin:$PATH
        export LD_LIBRARY_PATH=/deeptfactor/lib:$LD_LIBRARY_PATH
        cd deeptfactor
        python tf_running.py -i ../${samplename}.faa -o ./result -g cpu
        cp ./result/prediction_result.txt ../../tmp/annotation/tf/
        mkdir -p ../../result/Annotation/transcription_factors
        cp ./result/prediction_result.txt \
            ../../result/Annotation/transcription_factors
        cd ../
        
        if [[ -f "tmp/annotation/tf/prediction_result.txt" ]]; then
            echo "DeepTFactor done!" >> log
        else
            echo "DeepTFactor Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "SwissProt" ]; then
        # UniProt SwissProt
        mkdir -p tmp/annotation/uniprot
        diamond blastp -p ${threads} \
            -d ${database}/uniprot_sprot/uniprot_sprot.dmnd \
            -q ${faafile} \
            -o tmp/annotation/uniprot/${samplename}.swiss-prot.tab \
            -e 1e-5 -k 1 --max-hsps 1 --id 70 --query-cover 50 --subject-cover 50 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
            sstart send evalue bitscore qlen slen qcovhsp scovhsp
        
        if [[ -f "tmp/annotation/uniprot/${samplename}.swiss-prot.tab" ]]; then
            mkdir -p result/Annotation/uniprot
            cp tmp/annotation/uniprot/${samplename}.swiss-prot.tab result/Annotation/uniprot
            echo "SwissProt done!" >> log
        else
            echo "SwissProt Fail!" >> log
        fi
    fi

    if [ "$anno" == "DFVF" ]; then
        # DFVF (真菌毒力因子数据库)
        mkdir -p tmp/annotation/DFVF
        diamond blastp -p ${threads} \
            -d ${database}/DFVF/DFVF.dmnd \
            -q ${faafile} \
            -o tmp/annotation/DFVF/${samplename}.DFVF.txt \
            -e 1e-5 -k 1 --max-hsps 1 --id 70 --query-cover 50 --subject-cover 50 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
            sstart send evalue bitscore qlen slen qcovhsp scovhsp
        
        if [[ -f "tmp/annotation/DFVF/${samplename}.DFVF.txt" ]]; then
            mkdir -p result/Annotation/DFVF
            cp tmp/annotation/DFVF/${samplename}.DFVF.txt result/Annotation/DFVF
            echo "DFVF done!" >> log
        else
            echo "DFVF Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "FungalP450s" ]; then
        # FungalP450s
        mkdir -p tmp/annotation/FungalP450s
        diamond blastp -p ${threads} \
            -d ${database}/FungalP450s/FungalP450s.dmnd \
            -q ${faafile} \
            -o tmp/annotation/FungalP450s/${samplename}.FungalP450s.txt \
            -e 1e-5 -k 1 --max-hsps 1 --id 70 --query-cover 50 --subject-cover 50 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
            sstart send evalue bitscore qlen slen qcovhsp scovhsp
        
        if [[ -f "tmp/annotation/FungalP450s/${samplename}.FungalP450s.txt" ]]; then
            mkdir -p result/Annotation/FungalP450s
            cp tmp/annotation/FungalP450s/${samplename}.FungalP450s.txt result/Annotation/FungalP450s/FungalP450s.txt
            echo "FungalP450s done!" >> log
        else
            echo "FungalP450s Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "CAZy" ]; then
        # CAZy (Carbohyr)
        mkdir -p tmp/annotation/CAZy
        run_dbcan ${faafile} protein \
            --db_dir ${database}/CAZy/ \
            --out_dir tmp/annotation/CAZy \
            --dbCANFile dbCAN-HMMdb-V12.txt
        
        if [[ -f "tmp/annotation/CAZy/diamond.out" ]]; then
            mkdir -p result/Annotation/CAZy
            cp tmp/annotation/CAZy/diamond.out result/Annotation/CAZy/CAZy.txt
            echo "CAZy done!" >> log
        else
            echo "CAZy Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "antiSMASH" ]; then
        # antiSMASH (biosynthetic gene cluster)
        # masked genome file and gff file from step 1 are required. 
        if [[ -f "tmp/annotation/repeats/${samplename}.combined.masked.fasta" ]] && \
           [[ -f "tmp/annotation/predictgene/${samplename}.gff" ]]; then
            mkdir -p tmp/annotation/antismash
            maskfa="tmp/annotation/repeats/${samplename}.combined.masked.fasta"
            gff="tmp/annotation/predictgene/${samplename}.gff"
            
            antismash ${maskfa} --asf --rre --cb-subclusters --cb-knownclusters \
                --cpus ${threads} --taxon fungi --pfam2go \
                --output-dir tmp/annotation/antismash \
                --genefinding-gff3 ${gff}
            
            if [[ -f "tmp/annotation/antismash/${samplename}.combined.masked.zip" ]]; then
                mkdir -p result/Annotation/Antismash
                cp tmp/annotation/antismash/${samplename}.combined.masked.zip \
                    result/Annotation/Antismash/
                echo "antiSMASH done!" >> log
            else
                echo "antiSMASH Fail!" >> log
            fi
        else
            echo "antiSMASH skipped: missing maskfa or gff file" >> log
        fi
    fi
    
    if [ "$anno" == "EffectorP" ]; then
        # EffectorP (effector prediction)
        cp ${faafile} ${samplename}.faa
        sed -i "s#\*##g" ${samplename}.faa
        python EffectorP.py -o ${samplename}.effector.txt \
            -E ${samplename}.effector.fasta -i ${samplename}.faa
        
        if [[ -f "${samplename}.effector.txt" ]]; then
            mkdir -p result/Annotation/Effector
            cp ${samplename}.effector.txt \
                result/Annotation/Effector/Effector.txt
            cp ${samplename}.effector.fasta \
                result/Annotation/Effector/Effector.fasta
            num=$(grep '>' ${samplename}.effector.fasta | wc -l)
            echo "EffectorP: $num" >> Annotation.stat.txt
            echo "EffectorP done!" >> log
        else
            echo "EffectorP Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "Interproscan" ]; then
        # InterProScan
        cp ${faafile} ${samplename}.faa
        sed -i "s#\*##g" ${samplename}.faa
        cp -rf /opt/interproscan ./
        ln -s ${database}/interproscan-5.75-106.0/data/* interproscan/data
        mkdir -p tmpdir output
        
        ./interproscan/interproscan.sh -i ${faafile} -f tsv -d output \
            -dp --pathways -goterms -t p -cpu ${threads} -T ./tmpdir
        
        if [[ -f "output/*.tsv" ]]; then
            mkdir -p result/Annotation/Interproscan
            cp output/*.tsv \
                result/Annotation/Interproscan/Interproscan.txt
            echo "Interproscan done!" >> log
        else
            echo "Interproscan Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "KOG" ]; then
        # KOG (A curated Eukaryotic orthologous database )
        mkdir -p tmp/annotation/KOG
        diamond blastp -p ${threads} \
            -d ${database}/KOG/kog.dmnd \
            -q ${faafile} \
            -o tmp/annotation/KOG/${samplename}.KOG.txt \
            -e 1e-5 -k 1 --max-hsps 1 --id 40 --query-cover 40 --subject-cover 40 \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend \
            sstart send evalue bitscore qlen slen qcovhsp scovhsp
        
        if [[ -f "tmp/annotation/KOG/${samplename}.KOG.txt" ]]; then
            mkdir -p result/Annotation/KOG
            cp tmp/annotation/KOG/${samplename}.KOG.txt result/Annotation/KOG
            echo "KOG done!" >> log
        else
            echo "KOG Fail!" >> log
        fi
    fi
    
    if [ "$anno" == "Pfam" ]; then
        # Pfam
        export PERLLIB=/opt/PfamScan:$PERLLIB
        perl /opt/PfamScan/pfam_scan.pl -fasta ${faafile} \
            -dir ${database}/Pfam \
            -outfile Pfam_diamond.txt -cpu ${threads}
        
        if [[ -f "Pfam_diamond.txt" ]]; then
            mkdir -p result/Annotation/Pfam
           cp Pfam_diamond.txt result/Annotation/Pfam/Pfam.txt
            echo "Pfam done!" >> log
        else
            echo "Pfam Fail!" >> log
        fi
    fi
done


echo "Step 3 completed successfully!"
```

