import pandas as pd
import sys

# Task 1: Read in nodes.tsv and edges.tsv
nodes = pd.read_csv('nodes_test.tsv', sep='\t')
edges = pd.read_csv('edges_test.tsv', sep='\t')

# Task 2: For each compound, compute the number of genes that are BIND (CbG) to it
compound_gene_bind = edges[edges['metaedge'] == 'CbG']
compound_gene_counts = compound_gene_bind['source'].value_counts()
top_5_compound_gene_counts = compound_gene_counts.sort_values(ascending=False).head(5)
print(top_5_compound_gene_counts)

# Task 3: For each DISEASE, compute the number of GENE(s) that are UPREGULATES (DuG)
disease_gene_upregulate = edges[edges['metaedge'] == 'DuG']
disease_gene_counts = disease_gene_upregulate['source'].value_counts()
top_5_disease_gene_counts = disease_gene_counts.sort_values(ascending=False).head(5)
print(top_5_disease_gene_counts)

# Define the hash functions
def mid_square_hash(key, r):
    # Convert the string to an integer using ASCII values
    key_int = sum(ord(c) for c in key)
    # Square the key
    square = key_int ** 2
    # Convert to string for easy slicing
    square_str = str(square)
    # Get the middle r digits
    hash_value = square_str[len(square_str)//2 - r//2 : len(square_str)//2 + r//2]
    return int(hash_value)

def folding_hash(key, digit_size):
    # Convert the string to an integer using ASCII values
    key_int = sum(ord(c) for c in key)
    # Convert to string for easy slicing
    key_str = str(key_int)
    # Divide into parts and sum
    hash_value = sum(int(key_str[i:i+digit_size]) for i in range(0, len(key_str), digit_size))
    return hash_value


# Task 5: 
# Compute the hash tables using mid-square method or Folding Method
hash_tables = [{} for _ in range(10)]
for _, row in edges.iterrows():
    key = row['source']
    for i in range(10):
        if i < 5:
            hash_value = mid_square_hash(key, r=3 if i < 2 else 4)
        else:
            hash_value = folding_hash(key, digit_size=2 if i < 7 else 3)
        hash_tables[i][hash_value] = key

# Compare the sizes of the tables
for i, table in enumerate(hash_tables):
    print(f"Size of table {i+1}: {sys.getsizeof(table)} bytes")
