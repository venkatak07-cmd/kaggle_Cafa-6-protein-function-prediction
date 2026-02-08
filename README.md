# 🧬 CAFA 6 Protein Function Prediction - GO Term Propagation

[![Kaggle](https://img.shields.io/badge/Kaggle-Competition-blue?logo=kaggle)](https://www.kaggle.com/competitions/cafa-6-protein-function-prediction)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-GO_Analysis-green.svg)](https://geneontology.org/)

> **CAFA 6 Competition Solution**: Protein function prediction using Gene Ontology (GO) graph propagation with hard-max BFS algorithm and ontology-aware scoring.

## 🏆 Competition Overview

**CAFA 6** (Critical Assessment of protein Function Annotation) is an international challenge to predict protein function from amino acid sequences using the Gene Ontology (GO) framework.

| Metric | Approach |
|--------|----------|
| **Evaluation** | F-max, S-min, AUPR (per GO aspect) |
| **Ontology** | GO: BP (Biological Process), MF (Molecular Function), CC (Cellular Component) |
| **Method** | Graph-based score propagation |

## 📂 Repository Structure

```
├── go_propagation_pipeline.ipynb       # Main Kaggle notebook
├── data/
│   ├── Train/
│   │   ├── go-basic.obo               # GO ontology graph
│   │   ├── train_sequences.fasta
│   │   └── train_terms.tsv
│   ├── test_sequences.fasta
│   └── submission.tsv                 # Input predictions
├── submission_improved.tsv            # Output (propagated)
└── visualizations/                    # Generated plots
    ├── top_go_terms.png
    ├── score_distribution.png
    ├── new_terms_heatmap.png
    └── insights_dashboard.png
```

## 🚀 Quick Start

### 1. Setup (Kaggle Notebook)
```python
# Required packages
import pandas as pd
import numpy as np
import networkx as nx
from tqdm.auto import tqdm
from collections import deque
import matplotlib.pyplot as plt
import seaborn as sns
```

### 2. Run Pipeline
- **Runtime**: ~5-8 minutes (depending on submission size)
- **Memory**: ~4GB RAM
- **GPU**: Not required

## 🔬 Solution Approach

### Core Algorithm: Hard-Max BFS Propagation

The solution propagates predicted GO term scores through the ontology graph using **Breadth-First Search** with **hard-max** aggregation.

```python
# Propagation Logic
For each protein:
    Initialize queue with predicted GO terms and scores
    While queue not empty:
        Pop (term, score)
        For each parent in GO graph:
            new_score = score × edge_weight
            if new_score > current_parent_score:
                Update parent score
                Add parent to queue
```

### Graph Construction from OBO File

Parses `go-basic.obo` to build parent-child relationships:

| Relationship | Weight | Description |
|-------------|--------|-------------|
| **is_a** | 1.0 | Hierarchical inheritance (default) |
| **part_of** | 0.7 | Component relationship (softer weight) |

```python
# Example OBO parsing
id: GO:0008150
name: biological_process
is_a: GO:0005575 ! cellular_component
relationship: part_of GO:0005622 ! intracellular
```

### Key Parameters

```python
MIN_SCORE_THRESHOLD = 0.01     # Post-propagation cutoff
MAX_PROP_DEPTH = 8             # Prevents root flooding
PART_OF_WEIGHT = 0.7           # Softer than is_a
```

**Why These Values?**
- **Min Threshold**: Filters noise, keeps meaningful predictions
- **Max Depth**: Prevents over-propagation to root nodes (biological_process)
- **Part-of Weight**: Reflects weaker semantic relationship than is_a

## 🧮 Algorithm Details

### 1. Graph Building (`build_go_graph`)
```python
def build_go_graph(go_obo_path):
    parent_map = {}      # child -> [parents]
    edge_weight = {}     # (child, parent) -> weight

    # Parse OBO file
    # Extract is_a and part_of relationships
    # Assign weights: is_a=1.0, part_of=0.7

    return parent_map, edge_weight
```

**Output**: ~47,000 GO terms with hierarchical relationships

### 2. Hard-Max Propagation (`propagate_hard_max`)
```python
def propagate_hard_max(df, parent_map, edge_weight):
    for protein_id, group in df.groupby("protein_id"):
        scores = {}  # term -> max_score
        depth = {}   # term -> propagation_depth
        queue = deque(initial_predictions)

        while queue:
            term, score = queue.popleft()

            # Stop conditions
            if depth[term] >= MAX_PROP_DEPTH: continue
            if score < MIN_SCORE_THRESHOLD: continue

            # Propagate to parents
            for parent in parent_map[term]:
                new_score = score × edge_weight[(term, parent)]
                new_score = min(1.0, new_score)  # Clamp

                # Hard-max update
                if new_score > scores.get(parent, 0.0):
                    scores[parent] = new_score
                    depth[parent] = depth[term] + 1
                    queue.append((parent, new_score))
```

**Why Hard-Max?**
- ✅ Preserves strongest prediction path
- ✅ Prevents score dilution from multiple children
- ✅ Computationally efficient
- ✅ Biologically interpretable (most specific annotation)

### 3. Score Decay

```
Level 0 (Predicted): score = s
Level 1 (Parent):    score = s × 1.0 (is_a) or s × 0.7 (part_of)
Level 2 (Grandparent): score = s × 1.0² or s × 0.7²
...
Level N: score = s × weight^N
```

## 📊 Results & Visualizations

### Propagation Statistics

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| **Total Annotations** | X | Y | +Z% |
| **Avg GO Terms/Protein** | A | B | +C terms |
| **Unique GO Terms** | P | Q | +R new terms |

### Generated Visualizations

#### 1. Top GO Terms by Frequency
```python
top_terms = submission_improved['go_term'].value_counts().head(20)
# Bar plot of most propagated terms
```
**Insight**: Identifies high-level biological processes enriched across proteins

#### 2. Score Distribution
```python
sns.histplot(submission_improved['score'], bins=50, kde=True)
# Distribution of propagated confidence scores
```
**Insight**: Validates thresholding strategy, identifies score compression

#### 3. New Terms Added Per Protein
```python
new_terms = propagated_counts - original_counts
sns.histplot(new_terms_per_protein, bins=30, kde=True)
```
**Insight**: Shows propagation depth variability across proteins

#### 4. Protein-GO Term Heatmap
```python
heatmap_df = submission_improved.pivot_table(
    index='protein_id', 
    columns='go_term', 
    values='score', 
    fill_value=0
)
sns.heatmap(heatmap_df, cmap="cool")
```
**Insight**: Visualizes annotation density and protein functional similarity

#### 5. Insights Dashboard (4-Panel)
- New GO terms distribution
- Top proteins by new annotations
- GO term frequency (before/after)
- Comparative analysis

## 💻 Code Structure

```python
# 1. Configuration
OBO_PATH = "/kaggle/input/.../go-basic.obo"
MIN_SCORE_THRESHOLD = 0.01
MAX_PROP_DEPTH = 8
PART_OF_WEIGHT = 0.7

# 2. Build GO Graph
parent_map, edge_weight = build_go_graph(OBO_PATH)
# → ~47K terms with relationships

# 3. Load Predictions
submission = pd.read_csv(SUBMISSION_INPUT, sep="	")
# → protein_id, go_term, score

# 4. Propagate
submission_improved = propagate_hard_max(submission, parent_map, edge_weight)
# → Ancestor terms with propagated scores

# 5. Save & Visualize
submission_improved.to_csv(SUBMISSION_OUTPUT, sep="	")
# → Ready for CAFA evaluation
```

## 🎯 Key Features

### Ontology-Aware Propagation
- ✅ Respects GO DAG structure (Directed Acyclic Graph)
- ✅ Handles multiple inheritance (term → multiple parents)
- ✅ Differentiates relationship types (is_a vs part_of)
- ✅ Prevents cyclic propagation

### Computational Efficiency
- ✅ BFS with early stopping (depth limit)
- ✅ Score thresholding (prunes low-confidence paths)
- ✅ Protein-level parallelization (groupby)
- ✅ Progress tracking with tqdm

### Biological Soundness
- ✅ Hard-max preserves most specific annotations
- ✅ Score decay reflects semantic distance
- ✅ Part-of relationships weighted appropriately
- ✅ Root flooding prevention

## 🔧 Configuration Tuning

### Parameter Sensitivity

| Parameter | Low Value | High Value | Recommended |
|-----------|-----------|------------|-------------|
| **MIN_SCORE_THRESHOLD** | 0.001 (noisy) | 0.1 (strict) | **0.01** |
| **MAX_PROP_DEPTH** | 3 (shallow) | 15 (deep) | **8** |
| **PART_OF_WEIGHT** | 0.5 (weak) | 0.9 (strong) | **0.7** |

### Optimization Strategies

1. **Threshold Tuning**: Balance precision/recall trade-off
2. **Depth Limiting**: Prevent over-propagation to root
3. **Weight Adjustment**: Calibrate part_of vs is_a contribution
4. **Aspect-Specific**: Different parameters for BP/MF/CC

## 🚀 Improvements & Next Steps

### Potential Enhancements
- [ ] **Soft-Max Alternative**: Average scores from multiple children
- [ ] **Aspect-Specific Propagation**: Different weights for BP/MF/CC
- [ ] **Information Content Weighting**: Weight by GO term specificity
- [ ] **Neural Propagation**: Learnable propagation weights
- [ ] **Ensemble Propagation**: Combine multiple propagation strategies
- [ ] **Taxonomic Constraints**: Respect species-specific GO usage

### Advanced Features
- [ ] **IC-Based Filtering**: Filter by Information Content threshold
- [ ] **Semantic Similarity**: Resnik/Lin similarity for score calibration
- [ ] **Graph Embeddings**: Node2Vec/GraphSAGE for term representations
- [ ] **Probabilistic Propagation**: Bayesian belief propagation

### Evaluation Integration
- [ ] **F-max Optimization**: Tune for CAFA evaluation metric
- [ ] **Precision@K**: Top-K prediction analysis
- [ ] **Per-Aspect Analysis**: Separate BP/MF/CC optimization

## 📚 References

1. **CAFA 6 Competition**: https://www.kaggle.com/competitions/cafa-6-protein-function-prediction
2. **Gene Ontology**: http://geneontology.org/
3. **GO OBO Format**: http://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html
4. **CAFA Paper**: Radivojac et al. (2013). A large-scale evaluation of computational protein function prediction. Nature Methods.
5. **Graph Propagation**: Jiang et al. (2016). A new method for calculating information content of GO terms. Bioinformatics.

## 🛠️ Technical Details

### Dependencies
```
pandas >= 1.5.0
numpy >= 1.23.0
networkx >= 2.8.0
tqdm >= 4.65.0
matplotlib >= 3.6.0
seaborn >= 0.12.0
```

### Memory Optimization
- Streaming OBO parsing (line-by-line)
- Protein-level batch processing
- Score rounding (3 decimal places)
- Garbage collection between proteins

### Runtime Complexity
- **Graph Building**: O(N) where N = OBO file lines
- **Propagation**: O(P × T × D) where:
  - P = number of proteins
  - T = avg terms per protein
  - D = max propagation depth

## 🙋‍♂️ About 
**Competition**: CAFA 6 - Protein Function Prediction  
**License**: MIT

---

**How to Use**:
1. Upload base predictions to `submission.tsv`
2. Run propagation pipeline
3. Submit `submission_improved.tsv` to CAFA evaluation

**Questions?** Open an issue or discuss on Kaggle!
