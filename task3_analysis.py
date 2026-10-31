
# SWYNEX INTERNSHIP - TASK 3
# Visualization and Interpretation of Biological Data

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Entrez, SeqIO
from collections import Counter

# Download SARS-CoV-2 genome from NCBI
Entrez.email = "your_email@gmail.com"

handle = Entrez.efetch(
    db="nuccore",
    id="NC_045512.2",
    rettype="fasta",
    retmode="text"
)

record = SeqIO.read(handle, "fasta")
handle.close()

sequence = str(record.seq).upper()

# Calculate nucleotide counts
counts = Counter(sequence)

nucleotides = ["A", "T", "G", "C"]
values = [counts.get(n, 0) for n in nucleotides]

total = len(sequence)

gc_content = (
    (counts.get("G", 0) + counts.get("C", 0)) / total
) * 100

at_content = (
    (counts.get("A", 0) + counts.get("T", 0)) / total
) * 100

# Display results
print("Genome Length:", total)
print("Nucleotide Counts:", dict(zip(nucleotides, values)))
print("GC Content:", round(gc_content, 2), "%")
print("AT Content:", round(at_content, 2), "%")

# Chart 1: Nucleotide Composition
plt.figure(figsize=(8, 5))
plt.bar(nucleotides, values)
plt.title("Nucleotide Composition of SARS-CoV-2")
plt.xlabel("Nucleotides")
plt.ylabel("Number of Bases")
plt.savefig("Chart1_Nucleotide_Composition.png")
plt.show()

# Chart 2: GC and AT Content
plt.figure(figsize=(7, 5))
plt.pie(
    [gc_content, at_content],
    labels=["GC Content", "AT Content"],
    autopct="%1.2f%%"
)
plt.title("GC and AT Content")
plt.savefig("Chart2_GC_AT_Content.png")
plt.show()

# Chart 3: Nucleotide Frequency
frequencies = [(v / total) * 100 for v in values]

plt.figure(figsize=(8, 5))
plt.bar(nucleotides, frequencies)
plt.title("Nucleotide Frequency Distribution")
plt.xlabel("Nucleotides")
plt.ylabel("Frequency (%)")
plt.savefig("Chart3_Nucleotide_Frequency.png")
plt.show()

# Save results
results = pd.DataFrame({
    "Nucleotide": nucleotides,
    "Count": values,
    "Frequency_Percentage": frequencies
})

results.to_csv("Nucleotide_Analysis.csv", index=False)

# Save FASTA sequence
with open("Cleaned_SARS_CoV_2.fasta", "w") as file:
    file.write(">SARS-CoV-2_NC_045512.2\n")
    file.write(sequence + "\n")

print("Analysis completed successfully!")
