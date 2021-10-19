

from kedro.pipeline import Pipeline, node

from .nodes import marker_selector


def create_pipeline(**kwargs):
    return Pipeline(
        [
            node(
                marker_selector,
                'cell_data',
                'cell_markers',
                name="marker_selector",
            )
        ]
    )
