import re

# Define input and output file paths
input_file = "cleaned_requirements.txt"
output_file = "final_cleaned_requirements.txt"

# Read and process the file
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        cleaned_line = re.sub(r"=.*", "", line).strip()  # Remove everything after and including '@'
        outfile.write(cleaned_line + "\n")

print(f"Cleaned requirements saved to {output_file}")
