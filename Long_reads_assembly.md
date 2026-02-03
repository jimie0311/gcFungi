#############Long reads assembly################

```bash
bash fungi_assembly_long_reads.sh -inputdir -type <Pacbio/Nanopore/Pacbio_hifi>  -software <Flye/Hicanu/Hiasm>  -threads  -database <BUSCO database> [-size <genome size>] [-hicreads <path to HiC reads1, reads2; optional>]
```
```bash
#!/bin/bash

# data upload
while true; do
    case "$1" in
        -inputdir) inputdir=$2; shift 2;;
        -type) type=$2; shift 2;;
        -software) software=$2; shift 2;;
        -threads) threads=$2; shift 2;;
        -database) database=$2; shift 2;;
        -size) size=$2; shift 2;;
        -hicreads) hicreads=$2; shift 2;;
        -samplename) samplename=$2;shift 2;;
        *) break;;
    esac
done

# make tmp directory
mkdir -p tmp/inputfiles

########## 1. extract sample name
# check BAM file and format transformation
if [ -d "$inputdir" ]; then
    perl check_bam.pl ${inputdir}/ tmp/inputfiles
fi

########## 2. BAM to FASTA（if necessary）
if ls ${inputdir}/*.bam 1> /dev/null 2>&1; then
    ls ${inputdir}/*.bam > bam.list
    while read line; do
        ne=$(echo "$line" | awk -F '/' '{print $NF}' | sed "s#.bam##g")
        bam2fasta -o tmp/inputfiles/$ne.fasta $line
    done < bam.list
fi

readsdir="tmp/inputfiles"

########## 3. run assembly according to selected softwares
if [ "$software" == "Flye" ]; then
    # Flye组装
    export TMPDIR=/tmp
    ls ${readsdir}/* > list
    reads=$(cat list | perl -pe "s#\n# #g")
    
    if [[ "$type" =~ "Pacbio_hifi" ]]; then
        flye --pacbio-hifi $reads --out-dir ./out_flye --threads ${threads}
    elif [[ "$type" =~ "Pacbio" ]]; then
        flye --pacbio-raw $reads --out-dir ./out_flye --threads ${threads}
    else
        flye --nano-raw $reads --out-dir ./out_flye --threads ${threads}
    fi
    
    if [[ -f "out_flye/assembly.fasta" ]]; then
        cp out_flye/assembly.fasta tmp/${samplename}.fasta
        echo "Flye done!" >> log
    else
        echo "Flye Fail!" >> log
        exit 1
    fi
    
elif [ "$software" == "Hicanu" ]; then
    # Hicanu (Canu) assembly
    if [ -z "$size" ]; then
        echo "Error: Genome size required for Hicanu"
        exit 1
    fi
    
    ls ${readsdir}/* > list
    reads=$(cat list | perl -pe "s#\n# #g")
    
    if [[ "$type" =~ "Pacbio_hifi" ]]; then
        canu -p ${samplename} -d ./canu genomeSize=${size} -pacbio-hifi $reads maxThreads=${threads}
    elif [[ "$type" =~ "Pacbio" ]]; then
        canu -p ${samplename} -d ./canu genomeSize=${size} -pacbio $reads maxThreads=${threads}
    else
        canu -p ${samplename} -d ./canu genomeSize=${size} -nanopore $reads maxThreads=${threads}
    fi
    
    if [[ -f "canu/${samplename}.contigs.fasta" ]]; then
        cp canu/${samplename}.contigs.fasta tmp/${samplename}.fasta
        echo "Hicanu done!" >> log
    else
        echo "Hicanu Fail!" >> log
        exit 1
    fi
    
elif [ "$software" == "Hiasm" ]; then
    # Hifiasm assembly
    ls ${readsdir}/* > list
    reads=$(cat list | perl -pe "s#\n# #g")
    
    # Check if Hi-C data is available
    if [ -n "$hicreads" ] && [ $(echo $hicreads | tr ',' '\n' | wc -l) -ge 2 ]; then
        # run hifiasm with available Hi-C data
        hic1=$(echo $hicreads | cut -d',' -f1)
        hic2=$(echo $hicreads | cut -d',' -f2)
        
        if [[ "$type" =~ "Pacbio" ]]; then
            hifiasm -o ${samplename}.hifiasm -t ${threads} --h1 ${hic1} --h2 ${hic2} $reads
            gfa_file="${samplename}.hifiasm.hic.p_ctg.gfa"
        else
            hifiasm -o ${samplename}.hifiasm -t ${threads} --h1 ${hic1} --h2 ${hic2} --ont $reads
            gfa_file="${samplename}.hifiasm.hic.p_ctg.gfa"
        fi
    else
        # run Hifiasm without Hi-C
        if [[ "$type" =~ "Pacbio" ]]; then
            hifiasm -o ${samplename}.hifiasm -t ${threads} $reads
            gfa_file="${samplename}.hifiasm.bp.p_ctg.gfa"
        else
            hifiasm -o ${samplename}.hifiasm -t ${threads} --ont $reads
            gfa_file="${samplename}.hifiasm.bp.p_ctg.gfa"
        fi
    fi
    
    # GFA to FASTA
    if [[ -f "$gfa_file" ]]; then
        gfatools gfa2fa ${gfa_file} > ${samplename}.fasta
        cp ${samplename}.fasta tmp/${samplename}.fasta
        echo "Hifiasm done!" >> log
    else
        echo "Hifiasm Fail!" >> log
        exit 1
    fi
fi

########## 4. Assembly results selection
assemblyfasta="tmp/${samplename}.fasta"

if [[ ! -f "$assemblyfasta" ]]; then
    echo "Error: Assembly file not found!"
    exit 1
fi

########## 5. filter the short contig (<1000 bp)
mkdir -p spades_result
seqkit seq -m 1000 ${assemblyfasta} > spades_result/${samplename}.fasta

if [[ -f "spades_result/${samplename}.fasta" ]]; then
    echo "Seqkit done!" >> log
else
    echo "Seqkit Fail!" >> log
    exit 1
fi

########## 6. BUSCO completeness evaluation
busco -i spades_result/${samplename}.fasta -o ${samplename}_busco -m genome -c ${threads} -l ${database}

if [[ -f "${samplename}_busco/short_summary.specific.*.${samplename}_busco.txt" ]]; then
    buscofile=$(ls ${samplename}_busco/short_summary.specific.*.${samplename}_busco.txt)
    echo "BUSCO done!" >> log
else
    echo "BUSCO Fail!" >> log
    exit 1
fi

########## 7. QUAST evaluation
quast.py spades_result/${samplename}.fasta -o ${samplename}_quast -t ${threads}

if [[ -f "${samplename}_quast/report.txt" ]]; then
    quastfile="${samplename}_quast/report.txt"
    echo "QUAST done!" >> log
else
    echo "QUAST Fail!" >> log
    exit 1
fi

########## 8. Results integration

resultfile="${samplename}_result.txt"
echo "Sample: $samplename" > $resultfile
echo "Assembly: spades_result/${samplename}.fasta" >> $resultfile
echo "BUSCO: $buscofile" >> $resultfile
echo "QUAST: $quastfile" >> $resultfile

echo "All steps completed successfully!"
```