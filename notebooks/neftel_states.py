import numpy as np
import loompy
from typing import List

class NeftelStates:
    def __init__(self, genesets : List[np.ndarray]) -> None:
        self.genesets = genesets
    
    def fit(self, ds: loompy.LoomConnection) -> np.ndarray:
        genes = ds.ra.Gene
        n_genes = len(genes)
        print("Loaded %s genes" % n_genes)
        expression = ds[:,:] # (gene, cell)
        total_umi = ds.ca.TotalUMIs if "TotalUMIs" in ds.ca else np.sum(expression, axis=0, dtype='uint32')
        print("Loaded expression matrix (%s x %s) and total UMIs" % expression.shape)
        gene_total_umi = ds.ra.GeneTotalUMIs if "GeneTotalUMIs" in ds.ra else np.sum(expression, axis=1, dtype='uint32')
        print("Normalizing to median total UMI")
        expression = (expression / total_umi) * np.median(total_umi)
        expression = expression.T # (cell, gene)
        print("Calculating relative expression")
        expression = expression - np.mean(expression, axis=0)
        print("Ordering genes by absolute expression")
        expression_order = np.argsort(gene_total_umi)
        genes = genes[expression_order]
        expression = expression[:, expression_order]
        print("Computing signature scores")
        bin_size = n_genes / 30

        def gene_signature_score(signature_genes):
            score = np.mean(expression[:, np.isin(genes, signature_genes)], axis=1)
            control_genes = np.zeros(n_genes, dtype=bool)
            for g in signature_genes:
                ix = np.where(genes == g)[0][0]
                random100 = np.random.choice(int(bin_size), size=100, replace=False)
                bin_ix = ix // bin_size
                #print (ix, bin_ix, bin_size)
                control_idx = np.ceil(random100 + bin_ix * bin_size).astype('int')
                control_genes[control_idx] = True
            control_score = np.mean(expression[:, control_genes], axis=1)
            return score - control_score
        
        scores = list(map(gene_signature_score, self.genesets))
        return np.hstack([x[:, None] for x in scores])

if __name__ == '__main__':
    import sys, numpy
    import loompy
    genes = ["ACTB","ACTG1","ACTN1","ACTN2","ACTN3","ACTN4","AFDN","AKT1","AKT2","AKT3","AMOTL1","ASH1L","CASK","CDC42","CDK4","CGN","CLDN1","CLDN10","CLDN11","CLDN14","CLDN15","CLDN16",
"CLDN18","CLDN19","CLDN2","CLDN20","CLDN22","CLDN23","CLDN3","CLDN4","CLDN5","CLDN6","CLDN7","CRB3","CSNK2A1","CSNK2A2","CSNK2B","CTNNA1","CTNNA2",
"CTNNA3","CTNNB1","CTTN","EPB41","EPB41L1","EPB41L2","EPB41L3","EXOC3","EXOC4","F11R","GNAI1","GNAI2","GNAI3","HCLS1","HRAS","IGSF5","JAM2","JAM3","KRAS","LLGL1","LLGL2","MAGI1",
"MAGI2","MAGI3","MAP3K20","MPDZ","MRAS","MYH1","MYH10","MYH11","MYH13","MYH14","MYH15","MYH2","MYH3","MYH6","MYH7","MYH7B","MYH8","MYH9","MYL10","MYL12A","MYL12B",
"MYL5","MYL7","MYL9","NRAS","OCLN","PARD3","PARD6A","PARD6B","PARD6G","PATJ","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R2A","PPP2R2B","PPP2R2C","PPP2R2D","PRKCA",
"PRKCB","PRKCD","PRKCE","PRKCG","PRKCH","PRKCI","PRKCQ","PRKCZ","PTEN","RAB13","RAB3B","RHOA","RRAS","RRAS2","SPTAN1","SRC","SYMPK","TJAP1","TJP1","TJP2","TJP3","VAPA","YBX3","YES1"]
    geneset1 = numpy.array(genes)
    genesets = [ geneset1 ]
    with loompy.connect(sys.argv[1], "r") as d:
        nf = NeftelStates(genesets)
        scores = nf.fit(d)

