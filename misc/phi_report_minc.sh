cd "$1"
ls | while read ff; do mincheader "$ff" | grep -E 'patient\:|study\:'; done | sort | uniq > ../"$(basename "$1")"_phi_report.txt
