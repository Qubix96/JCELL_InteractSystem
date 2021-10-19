

from kedro.pipeline import Pipeline, node

from .nodes import load_data, inter_share, inter_analysis, mca_dim_red, plot_dendrogram, calculate_k, k_estimator, preapre_data_java, k_mean, marker_selection, pathways_reactoma, pathways_plot, venn_interactions, network


def create_pipeline(**kwargs):
    return Pipeline(
        [
            node(
                load_data,
                ['interaction_data' , 'cell_data'],
                ['binary_data','d'],
                name="load_data",
            ),
            node(
                inter_share,
                ['binary_data', 'd'],
                ['df', 'dfr'],
                name="inter_share",
            ),
            node(
                inter_analysis,
                'df',
                'analysis_dict',
                name="inter_analysis",
            ),
            node(
                mca_dim_red,
                'analysis_dict',
                ['tmp', 'label'],
                name="mca_dim_red",
            ),
            node(
                plot_dendrogram,
                ['label', 'tmp'],
                'den',
                name="plot_dendrogram",
            ),
            node(
                calculate_k,
                ['tmp', 'den'],
                ['sse', 'tmp3'],
                name="calculate_k",
            ),
            node(
                k_estimator,
                'sse',
                'k',
                name='k_estimator',
            ),
            node(
                preapre_data_java,
                'tmp3',
                ['valid','train', 'df1'],
                name="preapre_data_java",
            ),
            node(
                k_mean,
                ['k', 'train', 'valid', 'df1', 'analysis_dict'],
                'pred',
                name="k_mean",
            ),
            node(
                marker_selection,
                ['analysis_dict', 'pred', 'parameters'],
                ['results', 'analysis_clusters', 'clusters'],
                name="marker_selection",
            ),
            node(
                pathways_reactoma,
                ['results', 'tmp3'],
                'final_results',
                name="pathways_reactoma",
            ),
             node(
                pathways_plot,
                ['analysis_clusters', 'clusters', 'final_results'],
                'cluster_df',
                name="pathways_plot",
            ),
            node(
                venn_interactions,
                ['clusters', 'cluster_df'],
                None,
                name="venn_interactions",
            ),
             node(
                network,
                ['df', 'parameters', 'cell_data'],
                None,
                name="network",
            )
        ]
    )
