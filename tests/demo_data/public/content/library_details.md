# Demo Library Details

This is a **demonstration version** of the Molecular Library Web Application, showcasing the platform's capabilities with a curated set of 10 cyanoarene molecules.

## What You Can Explore

### 🔍 Chemical Search
- **Substructure Search**: Find molecules containing specific chemical patterns
- **Similarity Search**: Discover molecules similar to your query structure  
- **Interactive Drawing**: Use the built-in molecule editor to draw structures

### 📊 Data Visualization  
- **PCA Projections**: Explore molecular relationships in reduced dimensions
- **UMAP Visualization**: Alternative dimensionality reduction for clustering
- **Interactive Plots**: Click and explore molecules directly from visualizations

### 🧪 Molecular Properties
The demo includes DFT-calculated properties for each molecule:
- **HOMO Energy**: Highest Occupied Molecular Orbital
- **LUMO Energy**: Lowest Unoccupied Molecular Orbital  
- **Dipole Moment**: Molecular polarity measure
- **Chemical Hardness (η)**: Electronic property descriptor

### 🌐 3D Visualization
- **NGL Viewer**: Interactive 3D molecular structures
- **Conformer Data**: Multiple geometric conformations per molecule
- **Real-time Rendering**: Smooth 3D manipulation and visualization

### 📥 Data Export
- **CSV Downloads**: Export molecular data and properties
- **JSON Format**: Structured data for programmatic access
- **Batch Export**: Download multiple molecules at once

## Sample Dataset: Cyanoarenes

The demo features 10 carefully selected cyanoarene molecules that demonstrate:
- **Structural Diversity**: Various aromatic substitution patterns
- **Property Range**: Different electronic and steric characteristics  
- **Search Examples**: Good targets for substructure and similarity queries

## Technology Showcase

This demo illustrates the full technology stack:
- **FastAPI Backend**: High-performance Python web framework
- **React Frontend**: Modern JavaScript user interface
- **PostgreSQL + RDKit**: Chemical database with structure search
- **Docker Containers**: Scalable, reproducible deployment
- **Multi-tenant Architecture**: Namespace-based library separation

## Ready for Production

This same codebase powers production molecular libraries with thousands of compounds. The namespace-based architecture allows multiple independent libraries to run from a single deployment.

**Explore the demo to see how this platform can serve your molecular data!**