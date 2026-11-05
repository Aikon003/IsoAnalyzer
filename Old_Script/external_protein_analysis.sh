#!/bin/bash


set -e

# --- 1. ARGUMENT VALIDATION ---
if [ "$#" -ne 4 ]; then
    echo "Error: Invalid number of arguments."
    echo "Usage: $0 <full_path_nt.fasta> <full_path_AA.fasta> <Number_CPU> <Output_Directory>"
    echo "Example: ./phase_6_external_analysis.sh /analysis/fasta/isoform_nt.fasta /analysis/fasta/isoform_AA.fasta 48 /analysis/output/"
    exit 1
fi

NT_FASTA="$1"
AA_FASTA="$2"
CPU="$3"
OUTPUT_DIR="$4"

CONDA_ACTIVATE_SCRIPT="/opt/conda/bin/activate" 

if [ ! -f "$NT_FASTA" ]; then echo "Error: NT file not found: $NT_FASTA"; exit 1; fi
if [ ! -f "$AA_FASTA" ]; then echo "Error: AA file not found: $AA_FASTA"; exit 1; fi

echo "--- Phase 6: Starting External Analysis ---"
echo "Cores to use: $CPU"
echo "Output Directory: $OUTPUT_DIR"
echo "------------------------------------------"

# --- 2. STEP 6.1: CPAT (Coding Potential) ---
echo "[1/4] Starting CPAT..."

(
  echo "Activating 'cpat' environment..."
  source "$CONDA_ACTIVATE_SCRIPT" cpat
  
  
  cpat.py -g "$NT_FASTA" \
          -d /usr/share/cpat_data/Human_logitModel.RData \
          -x /usr/share/cpat_data/Human_Hexamer.tsv \
          -o "$OUTPUT_DIR/cpat_output.txt"
)
echo "CPAT completed."


# --- 3. STEP 6.2: Pfam (Protein Domains) ---
echo "[2/4] Starting Pfam (Protein Domains)..."

(
  echo "Activating 'pfam' environment..."
  source "$CONDA_ACTIVATE_SCRIPT" pfam
  
  
  pfam_scan.pl -fasta "$AA_FASTA" \
               -dir /usr/share/pfam_data \
               -outfile "$OUTPUT_DIR/pfam_output.txt" \
               -cpu "$CPU"
)
echo "Pfam completed."


# --- 4. STEP 6.3: SignalP (Signal Peptides) ---
echo "[3/4] Starting SignalP (Signal Peptides)..."

OriginalD=$(pwd)
mkdir -p "$OUTPUT_DIR/signalp_temp"

cd "/usr/share/signalp-5.0b/bin"

./signalp -fasta "$AA_FASTA" \
        -org euk \
        -format short \
        -prefix "$OUTPUT_DIR/signalp_output" \
        -tmp "$OUTPUT_DIR/signalp_temp"
rm -rf "$OUTPUT_DIR/signalp_temp"
echo "SignalP completed."

cd "$OriginalD"



# --- 5. STEP 6.4: IUPred2A (Disordered Regions) ---
echo "[4/4] Starting IUPred2A (using Python 3.7)..."
SPLIT_DIR="$OUTPUT_DIR/iupred_split_fasta"
IUPRED_FINAL_OUT="$OUTPUT_DIR/iupred2a_output.txt"
mkdir -p "$SPLIT_DIR"
rm -f "$IUPRED_FINAL_OUT"

echo "  Step A: Splitting FASTA file (using awk)..."

awk -F'>' 'NF>1 {
    gsub(/[|: ]/,"_", $2); 
    out=sprintf("%s/%s.fa", "'$SPLIT_DIR'", $2); 
    print ">"$2 > out; 
    close(out)
} 
NF==1 {
    if(out) print $0 >> out
}' "$AA_FASTA"

echo "  Step B: Running IUPred2A in parallel (using bash loop)..."
cd "$SPLIT_DIR"

IUPRED_PYTHON="/opt/conda/bin/python"
IUPRED_SCRIPT="/usr/share/iupred2a/iupred2a.py"
IUPRED_DATA="/usr/share/iupred2a/"

if [ ! -f "$IUPRED_SCRIPT" ]; then
    echo "ERROR: iupred2a.py not found at path $IUPRED_SCRIPT."
    exit 1
fi

count=0
Nproc="$CPU" 

echo "  Starting $Nproc parallel processes..."
for f in *.fa; do
   
    "$IUPRED_PYTHON" "$IUPRED_SCRIPT" -a -d "$IUPRED_DATA" "${f}" long > "iupred2a_output_${f}.txt" &
    
    let count+=1
    [[ $((count%Nproc)) -eq 0 ]] && echo "  ...waiting for process batch..." && wait 
done
wait 

echo "  Step C: Merging IUPred2A results..."
for f in iupred2a_output_*.fa.txt; do
    id=$(echo "${f}" | sed 's/iupred2a_output_//g' | sed 's/.fa.txt//g')
    sed "s/# POS\t/>${id}\n# POS\t/g" "${f}" >> "$IUPRED_FINAL_OUT"
    echo -e "\n\n################" >> "$IUPRED_FINAL_OUT"
done

cd "$OUTPUT_DIR"
rm -rf "$SPLIT_DIR"
echo "IUPred2A completed."



echo "------------------------------------------"
echo "--- Phase 6 Completed ---"
echo "All output files have been created in: $OUTPUT_DIR"