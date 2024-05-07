import pandas as pd
import sys
from functools import reduce

class DataProcessor:
    def __init__(self, nodes_file, edges_file):
        self.nodes = pd.read_csv(nodes_file, sep='\t')
        self.edges = pd.read_csv(edges_file, sep='\t')

    def compound_gene_counts(self):
        compound_gene_bind = self.edges[self.edges['metaedge'] == 'CbG']
        
        def mapper(compound_gene_pair):
            compound, gene = compound_gene_pair
            return [(compound, 1)]

        def reducer(counts, pair):
            compound, count = pair
            counts[compound] = counts.get(compound, 0) + count
            return counts

        mapped_data = [mapper((row['source'], row['target'])) for _, row in compound_gene_bind.iterrows()]
        mapped_data = reduce(lambda x, y: x + y, mapped_data)

        compound_gene_counts = reduce(reducer, mapped_data, {})
        sorted_compounds = sorted(compound_gene_counts.items(), key=lambda x: x[1], reverse=True)
        return sorted_compounds

    def disease_gene_upregulate_counts(self):
        disease_gene_upregulate = self.edges[self.edges['metaedge'] == 'DuG']

        def mapper(disease_gene_pair):
            disease, gene = disease_gene_pair
            return [(gene, 1)]

        def reducer(counts, pair):
            gene, count = pair
            counts[gene] = counts.get(gene, 0) + count
            return counts

        mapped_data = [mapper((row['source'], row['target'])) for _, row in disease_gene_upregulate.iterrows()]
        mapped_data = reduce(lambda x, y: x + y, mapped_data)

        disease_gene_counts = reduce(reducer, mapped_data, {})
        sorted_genes = sorted(disease_gene_counts.items(), key=lambda x: x[1], reverse=True)
        return sorted_genes

    def mid_square_hash(self, key, r):
        key_int = sum(ord(c) for c in key)
        square = key_int ** 2
        square_str = str(square)
        hash_value = square_str[len(square_str)//2 - r//2 : len(square_str)//2 + r//2]
        return int(hash_value)

    def folding_hash(self, key, digit_size):
        key_int = sum(ord(c) for c in key)
        key_str = str(key_int)
        hash_value = sum(int(key_str[i:i+digit_size]) for i in range(0, len(key_str), digit_size))
        return hash_value

    def compute_hash_tables(self, data, r_values):
        hash_tables = [[] for _ in range(len(r_values))]
        
        for compound, count in data:
            for i, r in enumerate(r_values):
                hash_value = self.mid_square_hash(compound, r=r)
                hash_tables[i % len(r_values)].append((hash_value, compound))
        
        return hash_tables

    def get_table_sizes(self, hash_tables):
        table_sizes = []
        for table in hash_tables:
            table_size = sys.getsizeof(table)
            for hash_value, key in table:
                table_size += sys.getsizeof(hash_value)
                table_size += sys.getsizeof(key)
            table_sizes.append(table_size)
        return table_sizes

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
    table_sizes_r4_compound = processor.get_table_sizes(hash_tables_r3_compound)
    print("Table sizes for r=3 compound data:", sorted(table_sizes_r3_compound))
    print("Table sizes for r=4 compound data:", sorted(table_sizes_r4_compound))
    
    # Compute and print table sizes (disease)
    table_sizes_r3_disease = processor.get_table_sizes(hash_tables_r3_disease)
    table_sizes_r4_disease = processor.get_table_sizes(hash_tables_r4_disease)
    print("Table sizes for r=3 disease data:", sorted(table_sizes_r3_disease))
    print("Table sizes for r=4 disease data:", sorted(table_sizes_r4_disease))