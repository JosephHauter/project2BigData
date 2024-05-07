import pandas as pd
import sys
from functools import reduce

class DataProcessor:
    def __init__(self, nodes_file, edges_file):
        self.nodes = pd.read_csv(nodes_file, sep='\t')
        self.edges = pd.read_csv(edges_file, sep='\t')

    def compound_gene_counts(self):
        # Filter the edges DataFrame to select only the edges with metaedge 'CbG'
        compound_gene_bind = self.edges[self.edges['metaedge'] == 'CbG']
        
        # Define a mapper function to extract gene names and count them
        def mapper(compound_gene_pair):
            compound, gene = compound_gene_pair
            return [(compound, 1)]

        # Define a reducer function to add the total counts for each genes
        def reducer(counts, pair):
            compound, count = pair
            counts[compound] = counts.get(compound, 0) + count
            return counts

        # Map each compound-gene pair 
        mapped_data = [mapper((row['source'], row['target'])) for _, row in compound_gene_bind.iterrows()]
        mapped_data = reduce(lambda x, y: x + y, mapped_data)

        # Reduce the mapped data by adding the counts for each gene
        compound_gene_counts = reduce(reducer, mapped_data, {})

        # Sort the genes based on their counts in descending order
        sorted_compounds = sorted(compound_gene_counts.items(), key=lambda x: x[1], reverse=True)
        
        # Return the sorted list of compounds
        return sorted_compounds

    def disease_gene_upregulate_counts(self):
        # Filter the edges DataFrame to select only the edges with metaedge 'DuG'
        disease_gene_upregulate = self.edges[self.edges['metaedge'] == 'DuG']

        # Define a mapper function to extract gene names and count them
        def mapper(disease_gene_pair):
            disease, gene = disease_gene_pair
            return [(gene, 1)]

        # Define a reducer function to add the total counts for each genes
        def reducer(counts, pair):
            gene, count = pair
            counts[gene] = counts.get(gene, 0) + count
            return counts

        # Map each disease-gene pair 
        mapped_data = [mapper((row['source'], row['target'])) for _, row in disease_gene_upregulate.iterrows()]
        # Reduce the mapped data by adding the counts for each gene
        mapped_data = reduce(lambda x, y: x + y, mapped_data)
        disease_gene_counts = reduce(reducer, mapped_data, {})
        # Sort the genes based on their counts in descending order
        sorted_genes = sorted(disease_gene_counts.items(), key=lambda x: x[1], reverse=True)

        # Return the sorted list of genes
        return sorted_genes

    def mid_square_hash(self, key, r):
        # Convert the characters of the key to their ASCII values and add them
        key_int = sum(ord(c) for c in key)
        # Square the key integer
        square = key_int ** 2
        # Convert the squared integer to a string
        square_str = str(square)
        # Extract a substring from the middle of the squared string, determined by r which represents the hash value
        hash_value = square_str[len(square_str)//2 - r//2 : len(square_str)//2 + r//2]
         # Convert the hash value back to an integer and return it
        return int(hash_value)

    def compute_hash_tables(self, data, r_values):
        # Initialize a list of hash tables, one for each r value
        hash_tables = [[] for _ in range(len(r_values))]
        
        # Iterate over each compound and its count in the data
        for compound, count in data:
            # Iterate over each r value and its index
            for i, r in enumerate(r_values):
                # Compute the hash value for the compound using the mid-square method
                hash_value = self.mid_square_hash(compound, r=r)
                # Append the compound and its hash value to the appropriate hash table
                # The hash table index is determined by the index of the r value in r_values
                hash_tables[i % len(r_values)].append((hash_value, compound))
        # Return the list of hash tables
        return hash_tables

    def get_table_sizes(self, hash_tables):
        # Initialize a list to store the sizes of each hash table
        table_sizes = []

        # Iterate over each hash table in the list of hash tables
        for table in hash_tables:
            # Initialize the size of the current hash table with the size of the list structure itself
            table_size = sys.getsizeof(table)
            
            table_sizes.append(table_size)
        # Return the list of sizes for each hash table
        return table_sizes
    
def compare_storage_sizes(total_size_r3, total_size_r4):
    if total_size_r3 < total_size_r4:
        return "Method r = 3 requires less storage."
    elif total_size_r3 > total_size_r4:
        return "Method r = 4 requires less storage."
    else:
        return "Both methods require the same amount of storage."

if __name__ == "__main__":
    #read in the data
    processor = DataProcessor('nodes_test.tsv', 'edges_test.tsv')

    # Task 2: Compute compound gene counts
    top_5_compounds = processor.compound_gene_counts()[:5]
    print(top_5_compounds)
    print("Top 5 compounds:")
    for compound, count in top_5_compounds:
        print(f"{compound}: {count} genes")

    # Task 3: Compute disease gene upregulate counts
    top_5_diseases = processor.disease_gene_upregulate_counts()[:5]
    print("Top 5 diseases:")
    for gene, count in top_5_diseases:
        print(f"{gene}: {count} diseases")

    # Task 5: Compute hash tables with r=3 and r=4
    compound_data = processor.edges[processor.edges['metaedge'] == 'CbG']
    disease_data = processor.edges[processor.edges['metaedge'] == 'DuG']
    r_values_3 = [3, 3, 3, 3, 3]
    r_values_4 = [4, 4, 4, 4, 4]
    
    hash_tables_r3_compound = processor.compute_hash_tables(top_5_compounds, r_values_3)
    hash_tables_r3_disease = processor.compute_hash_tables(top_5_diseases, r_values_3)
    hash_tables_r4_compound = processor.compute_hash_tables(top_5_compounds, r_values_4)
    hash_tables_r4_disease = processor.compute_hash_tables(top_5_diseases, r_values_4)

    # Compute and print table sizes (compound)
    table_sizes_r3_compound = processor.get_table_sizes(hash_tables_r3_compound)
    table_sizes_r4_compound = processor.get_table_sizes(hash_tables_r4_compound)
    table_sizes_r3_compound_sum = sum(table_sizes_r3_compound)
    table_sizes_r4_compound_sum = sum(table_sizes_r4_compound)

    print("Table sizes for r=3 compound data:", sorted(table_sizes_r3_compound))
    print("Table sizes for r=4 compound data:", sorted(table_sizes_r4_compound))
    print("Total size for r=3 compound data:", table_sizes_r3_compound_sum)
    print("Total size for r=4 compound data:", table_sizes_r4_compound_sum)
    print(compare_storage_sizes(table_sizes_r3_compound_sum, table_sizes_r4_compound_sum))
    
    # Compute and print table sizes (disease)
    table_sizes_r3_disease = processor.get_table_sizes(hash_tables_r3_disease)
    table_sizes_r4_disease = processor.get_table_sizes(hash_tables_r4_disease)
    table_sizes_r3_disease_sum = sum(table_sizes_r3_disease)
    table_sizes_r4_disease_sum = sum(table_sizes_r4_disease)

    print("Table sizes for r=3 disease data:", sorted(table_sizes_r3_disease))
    print("Table sizes for r=4 disease data:", sorted(table_sizes_r4_disease))
    print("Total size for r=3 disease data:", table_sizes_r3_disease_sum)
    print("Total size for r=4 disease data:", table_sizes_r4_disease_sum)
    print(compare_storage_sizes(table_sizes_r3_disease_sum, table_sizes_r4_disease_sum))