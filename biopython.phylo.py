from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import subprocess
import sys
import tempfile
import argparse

parser = argparse.ArgumentParser(description='Make a phylo SVG')
parser.add_argument('--distance_calculator', type = str, help = '', default = 'blastn')
parser.add_argument('--fasta_file', type = str, help = 'A fasta file that you want a phylogeny of', required = True)
parser.add_argument('--output_svg', type = str, help = 'output SVG file`')
parser.add_argument('--title', type = str, help = '')
args = parser.parse_args()

# made with help from Grok and Google Gemini

def execute(cmd):
	try:
		result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
		stdout = result.stdout
		stderr = result.stderr
		return [stdout, stderr]
	except subprocess.CalledProcessError as e:
		print(f"exit = {e.returncode}",	file=sys.stderr)
		print(f"STDOUT = {e.stdout}",	file=sys.stderr)
		print(f"STDERR = {e.stderr}", file=sys.stderr)
		raise RuntimeError(f"{cmd} failed") from e
		
def hide_inner_labels(node):
	if node.is_terminal():
		return node.name or ""  # Show terminal node names
	return ""  # Hide internal node labels
# Step 1: Read the FASTA file
fasta_file = args.fasta_file
sequences = list(SeqIO.parse(fasta_file, "fasta"))

# Step 2: Perform multiple sequence alignment with MUSCLE
tmp = tempfile.NamedTemporaryFile(suffix='.fa', dir='/tmp', delete=False)
tmp.close()
execute(f'muscle -align {fasta_file} -output {tmp.name}')
alignment = AlignIO.read(tmp.name, "fasta")

# Step 3: Calculate distance matrix and build phylogenetic tree
# https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html
calculator = DistanceCalculator('blastn')
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)  # Use Neighbor-Joining method

# Step 4: Draw and save the phylogenetic tree as SVG
fig = plt.figure(figsize=(12, 8))
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes, do_show = False, label_func=hide_inner_labels)
plt.title(args.title)
plt.savefig(args.output_svg, bbox_inches="tight")
