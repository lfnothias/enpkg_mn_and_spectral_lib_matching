import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import itertools

def plotter_count(dt_isdb_results_int, treemap_chemo_counted_results_path):

    """ Plotter functions
    Args:
        TODO
    Returns:
        TODO
    """
    dt_isdb_results_int = dt_isdb_results_int.replace({np.nan:'None'})
    fig = px.treemap(dt_isdb_results_int, path=['structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'],
                     color='structure_taxonomy_npclassifier_01pathway',
                     color_discrete_map={
                         'Terpenoids':'#E76F51',
                         'Alkaloids': '#264653',
                         'Amino acids and Peptides': '#287271',
                         'Polyketides': '#2A9D8F',
                         'Shikimates and Phenylpropanoids': '#E9C46A',
                         'Fatty acids': '#8AB17D',
                         'Carbohydrates': '#F4A261',})
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to number of annotations)")
    fig.update_annotations(font_size=12)
    fig.write_html(treemap_chemo_counted_results_path)


def plotter_intensity(dt_isdb_results_int, feature_table, treemap_chemo_intensity_results_path):

    """ Plotter functions
    Args:
        TODO
    Returns:
        TODO
    """
    dt_isdb_results_int = dt_isdb_results_int.replace({np.nan:'None'})
    dt_isdb_results_int = dt_isdb_results_int.merge(feature_table, left_on = 'feature_id', right_index=True, how='left')
    fig = px.treemap(dt_isdb_results_int, path=['structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class'],
                     color='intensity',
                     color_continuous_scale='RdBu_r',)
    fig.update_layout(margin = dict(t=50, l=25, r=25, b=25),
    title_text="Metabolite annotation overview (size proportional to the average features intensities)")
    fig.update_annotations(font_size=12)
    fig.write_html(treemap_chemo_intensity_results_path)

