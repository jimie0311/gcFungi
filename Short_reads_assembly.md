###########short_reads_assembly############
```bash
bash fungi_assembly.sh -outputdir  -configfile  -threads
```

### data example

- `-outputdir`: path to output direction
- `-configfile`: config file in CSV inlcuding sample ID and the path to PE1 and PE2

## Script

```bash
while true; do
        case "$1" in            
                -outputdir) outputdir=$2; shift 2;;
                -configfile) configfile=$2; shift 2;;
                -threads) threads=$2; shift 2;;
                *) break;
        esac
done

##########load_config File
while read line
do
#       echo $line
        OLD_IFS="$IFS"
        IFS=","
        tmp=($line)
        IFS="$OLD_IFS"
        name=(${name[@]} ${tmp[0]})
#       echo $name
        pe1=(${pe1[@]} ${tmp[1]})
        pe2=(${pe2[@]} ${tmp[2]})
done < $configfile

for ((i=0;i<${#name[@]};i++))
do
        file2=${name[${i}]}
        echo $file2
        outputDir=${outputdir}/
		mkdir -p ${outputDir}/trim_galore
		trim_galore -a AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -a2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG --stringency 10 -length 75 -e 0.1 -o ${outputDir}/trim_galore/ --paired ${pe1[${i}]} ${pe2[${i}]} --gzip
		read1=$(ls ${outputDir}/trim_galore/${file2}_1_val_1.fq.gz)
		read2=$(ls ${outputDir}/trim_galore/${file2}_2_val_2.fq.gz)
		if [[ -f "$read1" ]]
		then
			echo "TrimGalore done!">>${outputDir}/log
		else
			echo "TrimGalore Fail!">>${outputDir}/log
			break
		fi
	
		####assembly
		mkdir -p ${outputDir}/spades
		spades.py --only-assembler -t ${threads}  -m 100 --careful -1 ${read1} -2 ${read2} -o ${outputDir}/spades/${file2}_spades_out
		if [[ -f "${outputDir}/spades/${file2}_spades_out/scaffolds.fasta" ]]
		then
			echo "Spades done!">>${outputDir}/log
		else
			echo "Spades Fail!">>${outputDir}/log
			break
		fi
			
		mv ${outputDir}/spades/${file2}_spades_out/scaffolds.fasta ${outputDir}/spades/${file2}_spades_out//${file2}.fasta
		mkdir -p ${outputDir}/spades/spades_result
		seqkit seq -m 1000   ${outputDir}/spades/${file2}_spades_out//${file2}.fasta > ${outputDir}/spades/spades_result/${file2}.fasta

done
```