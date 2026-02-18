import loompy
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import tomotopy as tt
import logging
import pickle

test_components = np.arange(0, 40, 5)[1:]

if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    print("Arguments: LOOMFILE [NUM_COMPONENTS]+")
    print("Default number of components to generate: %s" % test_components)
    sys.exit(0)

# Get arguments from command line
loomfile = sys.argv[1]
if len(sys.argv) > 2:
    test_components = sys.argv[2].split(',')
    test_components = np.array(test_components).astype('int')
prefix = os.path.splitext(loomfile)[0] + "_LDA"
print("Output will be to subdirs %s/{n}/" % prefix)

human_ncRNA_path = "human_ncRNA.txt"
print("Reading human non-coding RNA from %s" % human_ncRNA_path)
ncRNA_human = [ l.strip() for l in open(human_ncRNA_path).readlines() ]

#def train_model(corpus, k, alpha, eta, max_iter=500, verbose=True):
def train_model(corpus, accessions, data, k, alpha, eta, max_iter=500, verbose=True):
    model = tt.LDAModel(k=k, alpha=alpha, eta=eta)
    for i in range(data.shape[1]):
        try:
            words = list(np.repeat(accessions, data[:, i].astype('int').T, axis=0))
            model.add_doc(words, ignore_empty_words=False)
        except:
            print(len(words))
    last_perp = 0
    for i in np.arange(0, max_iter, 10):
        model.train(iter=10, workers=0, parallel=tt.ParallelScheme.PARTITION)
        perp = model.perplexity
        if verbose:
            logging.info('Iteration: {}\tPerplexity: {}'.format(i, perp))
        # Escape loop when perplexity is stagnant
        if np.abs(last_perp - perp) < 0.5:
            break
        last_perp = perp
    return model

def plot_top_words(transformed, n_top_words, gene_names, filebase):
    offset = 0
    while offset < transformed.shape[0]:
        fig, axes = plt.subplots(5, 5, figsize=(15, 15), sharex=True)
        axes = axes.flatten()
        for topic_idx, topic in enumerate(transformed[offset: offset + 25]):
            top_features_ind = topic.argsort()[: -n_top_words - 1 : -1]
            top_features = [gene_names[i] for i in top_features_ind]
            weights = topic[top_features_ind]

            ax = axes[topic_idx]
            ax.barh(top_features, weights, height=0.7)
            ax.set_title(f"Topic {topic_idx + offset}", fontdict={"fontsize": 15})
            ax.invert_yaxis()
            ax.tick_params(axis="both", which="major", labelsize=10)
            for i in "top right left".split():
                ax.spines[i].set_visible(False)
        plt.subplots_adjust(top=0.90, bottom=0.05, wspace=0.90, hspace=0.3)
        plt.savefig(f'{filebase}_{offset}.png', dpi=144)
        offset += 25
        plt.close()

def plot_topics(transformed, xy, topics, gene_names, filebase):
    offset = 0
    n_factors = transformed.shape[1]
    while offset < n_factors:
        fig = plt.figure(figsize=(16, 16))
        fig.subplots_adjust(hspace=0, wspace=0)
        for nnc in range(offset, offset + 16):
            if nnc >= n_factors:
                break
            ax = plt.subplot(4, 4, (nnc - offset) + 1)
            plt.xticks(())
            plt.yticks(())
            plt.axis("off")
            plt.scatter(x=xy[:, 0], y=xy[:, 1], c='lightgrey', marker='.', alpha=0.5, s=20, lw=0)
            cells = np.random.permutation(transformed.shape[0])
            plt.scatter(x=xy[cells, 0], y=xy[cells, 1], c=transformed[:, nnc][cells], marker='.', alpha=0.5, s=20, cmap='viridis', lw=0)
            ax.text(.01, .99, '\n'.join(gene_names[topics[nnc].argsort()[::-1]][:5]), horizontalalignment='left', verticalalignment="top", transform=ax.transAxes)
            ax.text(.99, .9, f"{nnc}", horizontalalignment='right', transform=ax.transAxes, fontsize=12)
        plt.savefig(f"{filebase}_{offset}.png", dpi=144)
        offset += 16
        plt.close()

def plot_top_gene_exp(transformed, n, xy, ds, embedding, genes_to_plot, filebase):
    plt.figure(figsize=(16, 16))
    gs = plt.GridSpec(5, 5)
    plt.subplot(gs[0])
    plt.scatter(xy[:, 0], xy[:, 1], c=transformed[:, n], s=0.1)
    plt.axis('off')
    plt.title('Topic')
    for i, g in enumerate(genes_to_plot):
        plt.subplot(gs[i + 1])
        exp = np.log(ds[np.where(ds.ra.Gene == g)[0][0], :].flatten() + 1)
        cells = exp > 0
        plt.scatter(embedding[:, 0], embedding[:, 1], color='grey', alpha=0.2, s=0.1)
        plt.scatter(embedding[cells, 0], embedding[cells, 1], c=exp[cells], s=0.1)
        plt.axis('off')
        plt.title(g)
        plt.margins(0.01, 0.01)
    plt.tight_layout()
    plt.savefig(f'{filebase}.png', dpi=144)
    plt.close()


# Logging
logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

# Load data
with loompy.connect(loomfile, 'r') as ds:
    # Select genes that are expressed in > 100 cells and < 60% of all cells, and not non-coding
    print(f'Counting cells per genes')
    data = ds.sparse().astype('int')
    is_ncRNA = np.isin(ds.ra.Accession, ncRNA_human)
    nnz = data.getnnz(axis=1)    
    selected = (nnz > 100) & (nnz < data.shape[1] * 0.6) & (~is_ncRNA)
    accessions = ds.ra.Accession[selected]
    gene_names = ds.ra.Gene[selected]
    logging.info(f'Selecting {selected.sum()} genes')
    data = data.toarray()[selected, :]
    cellsums = np.sum(data, axis=0)
    nonzerocells = cellsums.flatten() > 0
    nonzerocellidx = np.nonzero(nonzerocells)[0]
    # Selected a subset of cells
    n_cells = min(len(nonzerocellidx), 200_000)
    logging.info(f'Selecting {n_cells} cells')
    cells = np.random.choice(nonzerocellidx, n_cells, replace=False)
    cells = np.sort(cells)
    print("Loading data")
    data = data[:, cells]
    # Get embedding for selected cells
    embedding = ds.ca.UMAP if 'UMAP' in ds.ca else ds.ca.Embedding if 'Embedding' in ds.ca else ds.ca.TSNE
    xy = embedding[cells, :]

# Create corpus by adding each cell's gene content as a document of (integer) words, one per UMI
logging.info("Creating corpus")
corpus = tt.utils.Corpus()
for i in range(data.shape[1]):
    corpus.add_doc(list(np.repeat(accessions, data[:, i].astype('int').T, axis=0)))

# Fit models for each k in test components
logging.info(f"Training with {test_components} topics")
coh = []
perplexity = []
for k in test_components:
    print("%s components..." % k)
    folder = f'{prefix}/{k}'
    if not os.path.exists(folder):
        os.makedirs(folder)
    # save cells
    np.savetxt(f'{folder}/cells.txt', cells)

    # train model with alpha and eta set to 1/k, like sklearn
    # model = train_model(corpus, k=k, alpha=1/k, eta=1/k, verbose=True)
    model = train_model(corpus, accessions, data, k=k, alpha=50/k, eta=0.1, verbose=True)

    # calculate coherence with 5 genes
    perplexity.append(model.perplexity)
    coh.append(
        tt.coherence.Coherence(model, coherence='u_mass', top_n=5).get_score()
        )

    # save into file
    model.save(f'{folder}/model.bin')

    # Normalize the components to estimate P(t|w) instead of P(w|t)
    topics = np.vstack([model.get_topic_word_dist(i, normalize=False) for i in range(k)])
    total_per_gene = topics.sum(axis=0)
    normed_topics = topics / total_per_gene

    # Get top genes
    logging.info("Filtering genes and plotting top scorers")
    filtered_top_genes = {}
    ix = np.array([np.where(accessions == x)[0][0] for x in list(model.vocabs)])
    model_gene_names = gene_names[ix]
    for n in range(k):
        # only keep genes that are > 100 and > 0.7 normed
        selected = (normed_topics[n] > 0.6) & (topics[n] > 500) 
        filtered_normed_topics = normed_topics[n, selected]
        ix = np.array([np.where(accessions == x)[0][0] for x in list(model.vocabs)])
        filtered_genes = model_gene_names[selected]
        filtered_top_genes[n] = filtered_genes[filtered_normed_topics.argsort()[::-1]]
        plot_top_words(topics, 10, model_gene_names, f"{folder}/top_gene_scores")

    # Transform data 
    logging.info("Transforming and plotting factors")
    transformed = np.vstack([doc.get_topic_dist() for doc in model.docs])
    # Plot topics
    plot_topics(transformed, xy, topics, model_gene_names, f"{folder}/factors")

    # Plot top genes on the tSNE - does not work since model training removes some docs, unknown which
    with loompy.connect(loomfile, 'r') as ds:
        # for every topic
        for n in range(k):
            # get top genes from unnormalized topics
            genes_to_plot = model_gene_names[topics[n].argsort()[::-1]][:24]
            plot_top_gene_exp(transformed, n, xy, ds, embedding, genes_to_plot, f"{folder}/top_gene_exp_{n}")
            # get top genes from normalized topics
            genes_to_plot = filtered_top_genes[n][:24]
            plot_top_gene_exp(transformed, n, xy, ds, embedding, genes_to_plot, f"{folder}/filtered_genes{n}")


# Plot results
logging.info("Plotting search results")
plt.figure(None, (5, 5))
plt.plot(coh, marker='o')
plt.xticks(range(len(coh)), test_components)
plt.xlabel('Number of components')
plt.ylabel('u_mass coherence')
plt.tight_layout()
plt.savefig(f'{prefix}/coherence.png', dpi=144)

plt.figure(None, (5, 5))
plt.plot(perplexity, marker='o')
plt.xticks(range(len(coh)), test_components)
plt.xlabel('Number of components')
plt.ylabel('perplexity')
plt.tight_layout()
plt.savefig(f'{prefix}/perplexity.png', dpi=144)

