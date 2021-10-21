

from kedro.pipeline import Pipeline, node

from .nodes import gene_relation, path_drug, plot


def create_pipeline(**kwargs):
    return Pipeline(
        [
            node(
                gene_relation,
                ['cell_data', 'parameters'],
                'first_relation_results',
                name="gene_relation",
            ),
            node(
                path_drug,
                'first_relation_results',
                'relation_results',
                name="path_drug",
            ),
            node(
                plot,
                ['relation_results','parameters'],
                None,
                name="plot",
            )
        ]
    )
