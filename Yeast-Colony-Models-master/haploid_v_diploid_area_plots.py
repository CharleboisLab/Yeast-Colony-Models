import pandas as pd
import plotly.graph_objects as go
import os
from pathlib import Path

## creates a violin plot comparing the areas of haploid and diploid colonies in rich and low nutrient conditions

def main():

    # set font sizes
    axis_font_size = 38
    xaxis_tick_font_size = 38
    yaxis_tick_font_size = 30
    legend_font_size = 30

    # set image dimensions
    height_int = 1300
    width_int = 1700

    diploidPathName = 'python_output/time_all_results_diploid_not_unipolar.csv'
    haploidPathName = 'python_output/time_all_results_haploid_not_unipolar.csv'

    diploidFile = Path(diploidPathName)
    haploidFile = Path(haploidPathName)

    diploidDF = pd.read_csv(diploidFile)
    diploidDF = diploidDF[diploidDF.diffusion == 'with_diffusion']

    haploidDF = pd.read_csv(haploidFile)
    haploidDF = haploidDF[haploidDF.diffusion == 'with_diffusion']

    ploidyList = ['diploid'] * len(diploidDF)
    diploidDF['ploidy'] = ploidyList

    ploidyList = ['haploid'] * len(haploidDF)
    haploidDF['ploidy'] = ploidyList

    df = pd.concat([diploidDF, haploidDF])

    fig = go.Figure()

    fig.add_trace(go.Violin(x=diploidDF['concentrations'][(
        diploidDF['directionType'] == 'none')],
        y=diploidDF['areaList'][(
            diploidDF['directionType'] == 'none')],
        legendgroup='none',
        name='diploid',
        box_visible=True,
        meanline_visible=True))

    fig.add_trace(go.Violin(x=haploidDF['concentrations'][(
        haploidDF['directionType'] == 'none')],
        y=haploidDF['areaList'][(
            haploidDF['directionType'] == 'none')],
        legendgroup='none',
        name='haploid',
        box_visible=True,
        meanline_visible=True))

    fig.update_traces(meanline_visible=True, points=False)
    fig.update_layout(yaxis_title='area (pixels)',
                      xaxis_title='nutrient concentration',
                      violinmode='group',
                      paper_bgcolor='rgba(\
                        255,255,255,255)',
                      plot_bgcolor='rgba(255,255,255,255)',
                      legend=dict(
                          font=dict(size=legend_font_size)),
                      height=height_int,
                      width=width_int)
    fig.update_xaxes(title_font=dict(size=axis_font_size),
                     tickfont=dict(
        size=xaxis_tick_font_size),
        mirror=True,
        ticks='outside',
        showline=True)
    fig.update_yaxes(title_font=dict(size=axis_font_size),
                     tickfont=dict(
        size=yaxis_tick_font_size),
        mirror=True,
        ticks='outside',
        showline=True,
        gridcolor='Grey')

    directory_name = './figures_output/compare_ploidy_figures'
    if not os.path.isdir(directory_name):
        directory_name = directory_name[2:]
        os.mkdir(directory_name)
    file_name = 'figures_output/compare_ploidy_figures/haploid_v_diploid_area_violin_plot.pdf'
    fig.write_image(file_name)


main()
