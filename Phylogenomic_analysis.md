###########Phylogenomic analysis############
#extract busco_full_table.tsv from each genome's busco folder
Please run BUSCO first following the steps introduced by https://busco.ezlab.org/. All busco result folders are named with sample_name.

## Script

```bash
#!/bin/bash

base_dir="/data/CMTG_genome_project/tree/"

#creat temp directory
temp_dir=$(mktemp -d)
cd "$temp_dir"

first_file=1

# scan all directories
for dir in "$base_dir"/*/; do
  # extract sample name (directory name)
  sample_name=$(basename "$dir")
  
  target_file="$dir/run_fungi_odb10/full_table.tsv"
  
  # check if every target sample exists
  if [ ! -f "$target_file" ]; then
    echo "warning: sample $target_file doesn't exist，skip this sample"
    continue
  fi
  
  filtered_content=$(grep -v '^#' "$target_file")
  
   if [ $first_file -eq 1 ]; then
    # extract the first two columns from the first file, uniq and add heads
    echo -e "busco\t$sample_name" > "${sample_name}.col"
    echo "$filtered_content" | awk '{print $1 "\t" $2}' | sort -u -k1,1 | awk '{print $1 "\t" $2}' >> "${sample_name}.col"
    first_file=0
  else
    # 
    echo -e "$sample_name" > "${sample_name}.col"
    echo "$filtered_content" | awk '{print $1 "\t" $2}' | sort -u -k1,1 | awk '{print $2}' >> "${sample_name}.col"
  fi
done

# merge all columns
paste *.col > "$base_dir/combined_busco.tsv"

rm -rf "$temp_dir"

echo "done successfully! the result is saved in $base_dir/combined_busco.tsv"

```
# select single_busco ID for phylogenomic analysis.
# Calculate the occurrence proportion of each BUSCO gene across all samples. Typically, select genes with a missing rate below 95%. Be sure to exclude samples where the BUSCO score is significantly lower than closely related taxa. Store the final determined single BUSCO gene IDs into the file single_busco.ID. When calculating the missing rate, do not include outgroup samples.

```bash
#!/bin/bash

base_path="/data/CMTG_genome_project/tree"
busco_ids_file="single_busco.ID"
#check if single_busco.ID exists
if \[\[ ! -f "`$base_path/$`busco\_ids\_file" ]]; then
echo "Error: `$base_path/$`busco\_ids\_file does not exist."
exit 1
fi

while read -r busco_id; do 
# name output_name 
output_file="$busco_id.faa"
> "$output_file"

# scan all folders under base_path
for sample_path in "$base_path"/*; do
    [[ -d "$sample_path" ]] || continue

    # ****** Warning：only extract folder name ******
    sample_id=$(basename "$sample_path")
        
        target_file="$sample_id/run_fungi_odb10/busco_sequences/single_copy_busco_sequences/$busco_id.faa"

                if [[ -f "$target_file" ]]; then
               sed "s/^>.*/>$sample_id/" "$target_file" >> "$output_file" 
        else
            echo ">$sample_id" >> "$output_file"
        fi
 done

echo "Processed $busco_id. Output written to $output_file"
```
## blast with mafft
# cut the unalligned head and end sequences and generate buscoID.aligned.fasta buscoID.trimed.fasta buscoID.trimed.phylip

```bash
#!/bin/bash
while read line
do
if [[ {line}.faa
mafft --auto --thread 50 --quiet {line}.aligned.fasta
/data/Liangjunmin/opt/biosoft/trimal/source/trimal -in {line}.trimed.fasta -gappyout
#/data/Liangjunmin/opt/biosoft/trimal/source/trimal -in {line}.trimed.fasta -gt 0.8 -st 0.001 -cons 80
#python fasta_to_phylip.py {line}.trimed.auto.phylip
python fasta_to_phylip.py {line}.trimed.phylip
#seqkit seq {line}.aligned.fasta
done<single_busco.ID
```
## merge all phylip files
```bash
#!/bin/env python3
import os

def merge_phylip_files(base_path, sample_list, output_file):
    # --- 1. read samples ---
	    with open(sample_list) as f:
        seen_lines = f.readlines()
    all_samples = list(dict.fromkeys(ln.strip() for ln in seen_lines if ln.strip()))
    expected_cnt = len(all_samples)          

    # --- 2. collect and merge all phylip files ---
    files = [f for f in os.listdir(base_path) if f.endswith('.phylip')]
    if not files:
        print("Error: 目录下找不到任何 .phylip 文件")
        return

    # --- 3. Init Container ---
    seqs = {s: [] for s in all_samples}
    gene_lengths = []
    seen_samples = set()      # existed samples

    # --- 4. check every sample in order ---
    for fname in files:
        fpath = os.path.join(base_path, fname)
        with open(fpath) as f:
            head = f.readline().strip().split()
            nseq, L = int(head[0]), int(head[1])
            gene_lengths.append(L)

            present = {}
            for _ in range(nseq):
                line = f.readline().rstrip('\n')
                if len(line) < 30:
                    line += ' ' * (30 - len(line))
                sid = line[:30].rstrip()
                present[sid] = line[30:].strip()
                seen_samples.add(sid)

            for s in all_samples:
                seqs[s].append(present.get(s, '-' * L))

       # --- 5. double check all expected samples are merged ---
    set_expect = set(all_samples)          
    set_actual = seen_samples              

    missing   = [s for s in all_samples if s not in set_actual]   
    redundant = [s for s in set_actual if s not in set_expect]    

    if missing or redundant:
        if missing:
            print(f"Warning: included in list but missed phylip file {len(missing)}，name：")
            print('\n'.join(missing))
        if redundant:
            print(f"Warning: identified phylip file but not included in list {len(redundant)}，name：")
            print('\n'.join(redundant))
    else:
        print("sample list is fully consistent with detected phylip files.")
    # --- 6. print merged file ---
    total_length = sum(gene_lengths)
    with open(output_file, 'w') as out:
        out.write(f"{len(all_samples)} {total_length}\n")
        for s in all_samples:
            full_seq = ''.join(seqs[s])
            out.write(f"{s.ljust(30)}{full_seq}\n")

    print(f"merged sucessfully！output is {output_file}")
    print(f"sample count：{len(all_samples)}，total length is{total_length}")
``

## run phylogenomic tree

```bash
 iqtree -s merged_sequences.phylip -m MFP -bb 1000 -bnni --seqtype AA -pre fungi_genome -T AUTO --threads-max 160 --safe --quiet
```