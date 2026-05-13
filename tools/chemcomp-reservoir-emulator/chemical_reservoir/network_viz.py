"""
Bipartite graph visualization for chemical reaction networks.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import networkx as nx
from networkx.algorithms import bipartite
from typing import Tuple, Dict, Optional

from .config import get_output_config, OutputConfig


def save_figure(fig, filename: str, config: OutputConfig = None):
    """Save figure to configured directory."""
    if config is None:
        config = get_output_config()
    
    if config.save_figures:
        filepath = config.get_path(filename)
        fig.savefig(filepath, format=config.figure_format, 
                   dpi=config.figure_dpi, bbox_inches='tight')
        print(f"  Saved: {filepath}")
    
    if config.show_figures:
        plt.show()
    else:
        plt.close(fig)


class ReactionNetworkVisualizer:
    """
    Visualize chemical reaction networks as bipartite graphs.
    """
    
    def __init__(self, reservoir):
        """Initialize with a reservoir instance."""
        self.reservoir = reservoir
        self.species_names = reservoir.species_names
        self.n_species = reservoir.n_species
        self.n_reactions = reservoir.n_reactions
        self.substrate_matrix = reservoir.substrate_matrix
        self.stoich_matrix = reservoir.stoich_matrix
        self.rate_constants = reservoir.rate_constants
    
    def build_directed_bipartite_graph(self) -> nx.DiGraph:
        """Build directed bipartite graph."""
        G = nx.DiGraph()
        
        # Add species nodes
        for sp in self.species_names:
            is_input = sp in self.reservoir.input_channels
            is_output = sp in self.reservoir.output_channels
            G.add_node(sp, bipartite=0, node_type='species',
                      is_input=is_input, is_output=is_output)
        
        # Add reaction nodes
        for rxn_idx in range(self.n_reactions):
            rxn_name = f"R{rxn_idx + 1}"
            G.add_node(rxn_name, bipartite=1, node_type='reaction',
                      rate=self.rate_constants[rxn_idx])
        
        # Add directed edges
        for rxn_idx in range(self.n_reactions):
            rxn_name = f"R{rxn_idx + 1}"
            
            for sp_idx, sp in enumerate(self.species_names):
                substrate_stoich = self.substrate_matrix[rxn_idx, sp_idx]
                product_stoich = substrate_stoich + self.stoich_matrix[rxn_idx, sp_idx]
                
                if substrate_stoich > 0:
                    G.add_edge(sp, rxn_name, stoich=substrate_stoich, edge_type='substrate')
                
                if product_stoich > 0:
                    G.add_edge(rxn_name, sp, stoich=product_stoich, edge_type='product')
        
        return G
    
    def get_bipartite_layout(self, G: nx.Graph) -> Dict:
        """Get layout positions for bipartite visualization."""
        species_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'species']
        reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'reaction']
        
        pos = {}
        
        for i, sp in enumerate(sorted(species_nodes)):
            pos[sp] = (0, -i)
        
        for i, rxn in enumerate(sorted(reaction_nodes)):
            pos[rxn] = (2, -i * (len(species_nodes) / max(len(reaction_nodes), 1)))
        
        return pos
    
    def plot_bipartite_network(self,
                               title: str = "Chemical Reaction Network (Bipartite)",
                               filename: str = "reaction_network_bipartite",
                               figsize: Tuple[float, float] = (14, 10),
                               config: OutputConfig = None) -> plt.Figure:
        """Plot the reaction network as a bipartite graph."""
        if config is None:
            config = get_output_config()
        
        G = self.build_directed_bipartite_graph()
        G_undirected = G.to_undirected()
        
        species_nodes = [n for n, d in G.nodes(data=True) if d.get('bipartite') == 0]
        is_bipartite_graph = bipartite.is_bipartite_node_set(G_undirected, set(species_nodes))
        
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # LEFT: Bipartite layout
        ax1 = axes[0]
        pos = self.get_bipartite_layout(G)
        
        species_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'species']
        reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'reaction']
        
        # Color species
        species_colors = []
        for sp in species_nodes:
            if sp in self.reservoir.input_channels and sp in self.reservoir.output_channels:
                species_colors.append('#9467bd')
            elif sp in self.reservoir.input_channels:
                species_colors.append('#2ca02c')
            elif sp in self.reservoir.output_channels:
                species_colors.append('#d62728')
            else:
                species_colors.append('#1f77b4')
        
        nx.draw_networkx_nodes(G, pos, nodelist=species_nodes,
                              node_color=species_colors, node_size=800, ax=ax1)
        
        nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes,
                              node_color='#ff7f0e', node_size=500, node_shape='s', ax=ax1)
        
        substrate_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'substrate']
        product_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'product']
        
        nx.draw_networkx_edges(G, pos, edgelist=substrate_edges,
                              edge_color='#1f77b4', arrows=True, arrowsize=15,
                              connectionstyle='arc3,rad=0.1', ax=ax1)
        
        nx.draw_networkx_edges(G, pos, edgelist=product_edges,
                              edge_color='#d62728', style='dashed', arrows=True, arrowsize=15,
                              connectionstyle='arc3,rad=-0.1', ax=ax1)
        
        nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold', ax=ax1)
        
        ax1.set_title(f"{title}\n(Bipartite: {is_bipartite_graph})", fontsize=14, fontweight='bold')
        ax1.axis('off')
        
        # Legend
        legend_elements = [
            mpatches.Patch(color='#2ca02c', label='Input Species'),
            mpatches.Patch(color='#d62728', label='Output Species'),
            mpatches.Patch(color='#9467bd', label='Input & Output'),
            mpatches.Patch(color='#1f77b4', label='Internal Species'),
            mpatches.Patch(color='#ff7f0e', label='Reaction'),
            plt.Line2D([0], [0], color='#1f77b4', linestyle='-', label='Substrate'),
            plt.Line2D([0], [0], color='#d62728', linestyle='--', label='Product'),
        ]
        ax1.legend(handles=legend_elements, loc='lower left', fontsize=9)
        
        # RIGHT: Spring layout
        ax2 = axes[1]
        pos_spring = nx.spring_layout(G, k=2, iterations=50, seed=42)
        
        nx.draw_networkx_nodes(G, pos_spring, nodelist=species_nodes,
                              node_color=species_colors, node_size=800, ax=ax2)
        
        nx.draw_networkx_nodes(G, pos_spring, nodelist=reaction_nodes,
                              node_color='#ff7f0e', node_size=500, node_shape='s', ax=ax2)
        
        nx.draw_networkx_edges(G, pos_spring, edgelist=substrate_edges,
                              edge_color='#1f77b4', arrows=True, arrowsize=15, ax=ax2)
        
        nx.draw_networkx_edges(G, pos_spring, edgelist=product_edges,
                              edge_color='#d62728', style='dashed', arrows=True, arrowsize=15, ax=ax2)
        
        nx.draw_networkx_labels(G, pos_spring, font_size=10, font_weight='bold', ax=ax2)
        
        ax2.set_title("Spring Layout View", fontsize=14, fontweight='bold')
        ax2.axis('off')
        
        plt.tight_layout()
        save_figure(fig, filename, config)
        
        return fig
    
    def plot_detailed_bipartite(self,
                                filename: str = "reaction_network_detailed",
                                figsize: Tuple[float, float] = (16, 12),
                                config: OutputConfig = None) -> plt.Figure:
        """Create detailed bipartite visualization."""
        if config is None:
            config = get_output_config()
        
        G = self.build_directed_bipartite_graph()
        
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
        
        # TOP LEFT: Bipartite layout
        ax1 = fig.add_subplot(gs[0, 0])
        pos = self.get_bipartite_layout(G)
        
        species_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'species']
        reaction_nodes = [n for n, d in G.nodes(data=True) if d.get('node_type') == 'reaction']
        
        species_colors = []
        for sp in species_nodes:
            if sp in self.reservoir.input_channels and sp in self.reservoir.output_channels:
                species_colors.append('#9467bd')
            elif sp in self.reservoir.input_channels:
                species_colors.append('#2ca02c')
            elif sp in self.reservoir.output_channels:
                species_colors.append('#d62728')
            else:
                species_colors.append('#1f77b4')
        
        nx.draw_networkx_nodes(G, pos, nodelist=species_nodes,
                              node_color=species_colors, node_size=700, ax=ax1)
        nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes,
                              node_color='#ff7f0e', node_size=400, node_shape='s', ax=ax1)
        
        substrate_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'substrate']
        product_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'product']
        
        nx.draw_networkx_edges(G, pos, edgelist=substrate_edges,
                              edge_color='#1f77b4', arrows=True, arrowsize=12, ax=ax1)
        nx.draw_networkx_edges(G, pos, edgelist=product_edges,
                              edge_color='#d62728', style='dashed', arrows=True, arrowsize=12, ax=ax1)
        
        nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax1)
        
        ax1.set_title("Bipartite Representation", fontsize=12, fontweight='bold')
        ax1.axis('off')
        
        # TOP RIGHT: Stoichiometry matrix
        ax2 = fig.add_subplot(gs[0, 1])
        
        biadj = np.zeros((len(species_nodes), len(reaction_nodes)))
        species_idx = {sp: i for i, sp in enumerate(sorted(species_nodes))}
        reaction_idx = {rxn: i for i, rxn in enumerate(sorted(reaction_nodes))}
        
        for u, v, d in G.edges(data=True):
            if d.get('edge_type') == 'substrate':
                sp, rxn = u, v
                biadj[species_idx[sp], reaction_idx[rxn]] += d.get('stoich', 1)
            elif d.get('edge_type') == 'product':
                rxn, sp = u, v
                biadj[species_idx[sp], reaction_idx[rxn]] -= d.get('stoich', 1)
        
        im = ax2.imshow(biadj, cmap='RdBu', aspect='auto', vmin=-3, vmax=3)
        ax2.set_xticks(range(len(reaction_nodes)))
        ax2.set_xticklabels(sorted(reaction_nodes), rotation=45, ha='right', fontsize=9)
        ax2.set_yticks(range(len(species_nodes)))
        ax2.set_yticklabels(sorted(species_nodes), fontsize=9)
        ax2.set_xlabel('Reactions', fontsize=11)
        ax2.set_ylabel('Species', fontsize=11)
        ax2.set_title('Stoichiometry Matrix', fontsize=12, fontweight='bold')
        plt.colorbar(im, ax=ax2, label='Net Stoichiometry')
        
        # BOTTOM LEFT: Degree distribution
        ax3 = fig.add_subplot(gs[1, 0])
        
        species_degrees = [G.degree(sp) for sp in species_nodes]
        
        x_sp = np.arange(len(species_nodes))
        ax3.bar(x_sp, species_degrees, color='#1f77b4', alpha=0.7)
        ax3.set_xticks(x_sp)
        ax3.set_xticklabels(sorted(species_nodes), rotation=45, ha='right', fontsize=9)
        ax3.set_ylabel('Degree', fontsize=11)
        ax3.set_title('Species Connectivity', fontsize=12, fontweight='bold')
        ax3.grid(alpha=0.3, axis='y')
        
        # BOTTOM RIGHT: Statistics
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.axis('off')
        
        G_undirected = G.to_undirected()
        
        stats_text = (
            f"NETWORK STATISTICS\n"
            f"{'='*40}\n\n"
            f"Nodes:\n"
            f"  Species: {len(species_nodes)}\n"
            f"  Reactions: {len(reaction_nodes)}\n\n"
            f"Edges:\n"
            f"  Total: {G.number_of_edges()}\n"
            f"  Substrate: {len(substrate_edges)}\n"
            f"  Product: {len(product_edges)}\n\n"
            f"I/O Configuration:\n"
            f"  Input: {self.reservoir.input_channels}\n"
            f"  Output: {self.reservoir.output_channels}\n"
        )
        
        ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.suptitle('Chemical Reaction Network Analysis', fontsize=14, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, config)
        
        return fig
    
    def plot_reaction_equations(self,
                               filename: str = "reaction_equations",
                               figsize: Tuple[float, float] = (10, 8),
                               config: OutputConfig = None) -> plt.Figure:
        """Create figure showing chemical reaction equations."""
        if config is None:
            config = get_output_config()
        
        fig, ax = plt.subplots(figsize=figsize)
        ax.axis('off')
        
        reactions_text = "CHEMICAL REACTION NETWORK\n" + "=" * 50 + "\n\n"
        
        for rxn_idx in range(self.n_reactions):
            substrates = []
            products = []
            
            for sp_idx, sp in enumerate(self.species_names):
                substrate_stoich = self.substrate_matrix[rxn_idx, sp_idx]
                product_stoich = substrate_stoich + self.stoich_matrix[rxn_idx, sp_idx]
                
                if substrate_stoich > 0:
                    if substrate_stoich == 1:
                        substrates.append(sp)
                    else:
                        substrates.append(f"{int(substrate_stoich)}{sp}")
                
                if product_stoich > 0:
                    if product_stoich == 1:
                        products.append(sp)
                    else:
                        products.append(f"{int(product_stoich)}{sp}")
            
            substrate_str = " + ".join(substrates) if substrates else "∅"
            product_str = " + ".join(products) if products else "∅"
            rate = self.rate_constants[rxn_idx]
            
            rxn_str = f"R{rxn_idx + 1}:  {substrate_str}  →  {product_str}    (k = {rate:.4f})\n"
            reactions_text += rxn_str
        
        reactions_text += "\n" + "=" * 50 + "\n"
        reactions_text += f"\nInput Species: {', '.join(self.reservoir.input_channels)}\n"
        reactions_text += f"Output Species: {', '.join(self.reservoir.output_channels)}\n"
        
        ax.text(0.05, 0.95, reactions_text, transform=ax.transAxes,
               fontsize=11, verticalalignment='top', fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='white', edgecolor='black'))
        
        plt.title('Reaction Network Equations', fontsize=14, fontweight='bold')
        
        save_figure(fig, filename, config)
        
        return fig