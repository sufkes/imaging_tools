#!/bin/bash

in_path="$1"

# Remove carriage returns from file.

# Backup file.
backup_path="$in_path".bak
cp -n "$in_path" "$backup_path"

# Remove carriage returns from file.
sed 's/\r//' "$backup_path" > "$in_path"
