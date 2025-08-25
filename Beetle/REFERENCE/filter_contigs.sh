#!/bin/bash

# extract_sequences.sh - Extract FASTA sequences based on contig IDs
# Usage: ./extract_sequences.sh input.fasta ids.txt > output.fasta

set -euo pipefail

# Function to display usage information
usage() {
    echo "Usage: $0 <input.fasta> <ids.txt>"
    echo ""
    echo "Extract FASTA sequences matching contig IDs from ids.txt"
    echo ""
    echo "Arguments:"
    echo "  input.fasta  - Input FASTA file containing sequences"
    echo "  ids.txt      - Text file with one contig ID per line"
    echo ""
    echo "Output:"
    echo "  Extracted sequences in FASTA format to stdout"
    echo ""
    echo "Example:"
    echo "  $0 genome.fasta contig_ids.txt > filtered_sequences.fasta"
    exit 1
}

# Function to log messages to stderr
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Main execution
main() {
    # Check for help flag
    if [[ "${1:-}" == "-h" ]] || [[ "${1:-}" == "--help" ]]; then
        usage
    fi
    
    # Get command line arguments
    local fasta_file="${1:-}"
    local ids_file="${2:-}"
    
    # Check if both arguments are provided
    if [[ -z "$fasta_file" ]] || [[ -z "$ids_file" ]]; then
        log "ERROR: Both FASTA file and IDs file must be provided"
        usage
    fi
    
    # Check if files exist
    if [[ ! -f "$fasta_file" ]]; then
        log "ERROR: FASTA file '$fasta_file' not found"
        exit 2
    fi
    
    if [[ ! -f "$ids_file" ]]; then
        log "ERROR: IDs file '$ids_file' not found"
        exit 2
    fi
    
    # Create clean IDs file
    local temp_ids=$(mktemp)
    trap "rm -f '$temp_ids'" EXIT
    
    log "Reading contig IDs from '$ids_file'..."
    
    # Extract clean IDs
    sed 's/^[[:space:]]*>//;s/[[:space:]].*//' "$ids_file" | grep -v '^$' > "$temp_ids"
    
    local id_count=$(wc -l < "$temp_ids")
    log "Loaded $id_count target contig IDs"
    
    if [[ $id_count -eq 0 ]]; then
        log "ERROR: No valid IDs found in '$ids_file'"
        exit 2
    fi
    
    log "Processing FASTA file '$fasta_file'..."
    
    # Use awk to extract sequences
    awk '
    BEGIN {
        # Read target IDs into array
        while ((getline id < "'$temp_ids'") > 0) {
            target[id] = 1
        }
        close("'$temp_ids'")
        found = 0
    }
    /^>/ {
        # Extract ID from header (first word after >)
        gsub(/^>[ \t]*/, "")
        split($0, parts, /[ \t]/)
        current_id = parts[1]
        
        # Check if this ID is a target
        if (current_id in target) {
            print_seq = 1
            print ">" $0
            found++
        } else {
            print_seq = 0
        }
        next
    }
    {
        # Print sequence lines if we should extract this sequence
        if (print_seq) {
            print
        }
    }
    END {
        print "Found " found " sequences" > "/dev/stderr"
    }
    ' "$fasta_file"
    
    log "Extraction complete"
}

# Run main function with all arguments
main "$@"
