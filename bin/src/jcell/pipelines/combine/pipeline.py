

from kedro.pipeline import Pipeline, node

from .nodes import combine_data, venn, summary_path


def create_pipeline(**kwargs):
    return Pipeline(
        [
            node(
                combine_data,
                ['results_interaction' , 'results_relation', 'parameters'],
                'df_summary',
                name="combine_data",
            ),
            node(
                venn,
                ['results_interaction', 'results_relation', 'parameters'],
                None,
                name="venn",
            ),
            node(
                summary_path,
                ['df_summary', 'parameters'],
                None,
                name="summary_path",
            )
        ]
    )
