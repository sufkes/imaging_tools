#!/bin/bash
# Run this command after running rsnAtlasLabels.sh to combine the results in a CSV file. 
rr=./report.csv; rm -f "$rr"; for ii in $(seq 1 25); do ff=ic${ii}.txt; echo -n "${ii}," >> "$rr"; cat "$ff" | cut -d " " -f 4- | tr '\n' ', ' >> "$rr"; echo >> "$rr"; done
