#!/bin/bash

# Script to remove .simpo from filenames in the current directory

# Loop through all files containing .simpo in their name
for file in *.simpo.*; do
    # Check if the pattern matched any files
    if [ ! -e "$file" ]; then
        echo "No files with .simpo found in current directory"
        break
    fi
    
    # Generate new filename by removing .simpo
    newfile="${file/.simpo/}"
    
    # Check if target file already exists
    if [ -e "$newfile" ]; then
        echo "Warning: $newfile already exists, skipping $file"
        continue
    fi
    
    # Rename the file
    mv "$file" "$newfile"
    echo "Renamed: $file -> $newfile"
done

echo "Done!"
