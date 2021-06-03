import pandas as pd
import plotly.graph_objects as go
import os
from pathlib import Path

# Creates split violin plots for all measurements except area


def main():

    ploidy_list = ['haploid', 'diploid']
    budding_list = ['not_unipolar']
    diffusion_list = ['with_diffusion']
    end_list = ['time', 'count']
    
    y_axis_measurements = ['finalCellCounts',
                           'roundnessList',
                           'compactnessList',
                           'convexityList',
                           'boundary_fluctuations',
                           'perimeterList',
                           'solidityList',
                           'elongationList',
                           'areaList',
                           'timeList']

    y_axis_labels = dict([
        ('finalCellCounts', 'number of cells'),
        ('roundnessList', 'roundness'),
        ('compactnessList', 'compactness'),
        ('convexityList', 'convexity'),
        ('boundary_fluctuations', 'boundary fluctuations'),
        ('timeList', 'timestep'),
        ('perimeterList', 'perimeter  (pixels)'),
        ('solidityList', 'solidity'),
        ('elongationList', 'elongation'),
        ('areaList', 'colony area  (pixels)')
    ])

    y_axis_file_name = dict([
        ('finalCellCounts', 'cellno'),
        ('roundnessList', 'roundness'),
        ('compactnessList', 'compactness'),
        ('convexityList', 'convexity'),
        ('boundary_fluctuations', 'fluctuations'),
        ('timeList', 'time'),
        ('perimeterList', 'perimeter'),
        ('solidityList', 'solidity'),
        ('elongationList', 'elongation'),
        ('areaList', 'area')
    ])

    x_axis_labels = dict([
        ('no_diffusion', 'nutrient concentration (no diffusion)'),
        ('with_diffusion', 'nutrient concentration')
    ])

    # set font sizes
    axis_font_size = 38
    xaxis_tick_font_size = 38
    yaxis_tick_font_size = 30
    legend_font_size = 30

    # set image dimensions
    height_int = 1300
    width_int = 1700

    for end_condition in end_list:
        for ploidy in ploidy_list:
            for bud in budding_list:
                path_name = 'python_output/' + end_condition + \
                    '_all_results_'+ploidy+'_'+bud+'.csv'
                print(path_name)
                file = Path(path_name)

                for measurement in y_axis_measurements:
                    for diff in diffusion_list:
                        df = pd.read_csv(file)
                        df = df[df.diffusion == diff]

                        # make violin plot
                        fig = go.Figure()

                        # weak diagonal MF trace
                        fig.add_trace(go.Violin(x=df['concentrations'][(
                            df['directionType'] == 'diagonal') & (
                            df['strength'] == 'weak')],
                            y=df[measurement][(
                                df['directionType'] == 'diagonal') & (
                                df['strength'] == 'weak')],
                            legendgroup='diagonal',
                            name='weak diagonal MFs',
                            box_visible=True,
                            meanline_visible=True,
                            side='negative',
                            offsetgroup='diagonal'))

                        # strong diagonal MF trace
                        fig.add_trace(go.Violin(x=df['concentrations'][(
                            df['directionType'] == 'diagonal') & (
                            df['strength'] == 'strong')],
                            y=df[measurement][(
                                df['directionType'] == 'diagonal') & (
                                df['strength'] == 'strong')],
                            legendgroup='diagonal',
                            name='strong diagonal MFs',
                            box_visible=True,
                            meanline_visible=True,
                            side='positive',
                            offsetgroup='diagonal'))

                        # weak axial MF trace
                        fig.add_trace(go.Violin(x=df['concentrations'][(
                            df['directionType'] == 'axis') & (
                            df['strength'] == 'weak')],
                            y=df[measurement][(
                                df['directionType'] == 'axis') & (
                                df['strength'] == 'weak')],
                            legendgroup='axis',
                            name='weak axial MFs',
                            box_visible=True,
                            meanline_visible=True,
                            side='negative',
                            offsetgroup='axis'))

                        # strong axial MF trace
                        fig.add_trace(go.Violin(x=df['concentrations'][(
                            df['directionType'] == 'axis') & (
                            df['strength'] == 'strong')],
                            y=df[measurement][(
                                df['directionType'] == 'axis') & (
                                df['strength'] == 'strong')],
                            legendgroup='axis',
                            name='strong axial MFs',
                            box_visible=True,
                            meanline_visible=True,
                            side='positive',
                            offsetgroup='axis'))

                        # no MF trace
                        fig.add_trace(go.Violin(x=df['concentrations'][(
                            df['directionType'] == 'none')],
                            y=df[measurement][(
                                df['directionType'] == 'none')],
                            legendgroup='none',
                            name='no MFs',
                            box_visible=True,
                            meanline_visible=True,
                            offsetgroup='none'))

                        # formatting the figure
                        fig.update_traces(meanline_visible=True, points=False)
                        fig.update_layout(yaxis_title=y_axis_labels[
                                          measurement],
                                          xaxis_title=x_axis_labels[diff],
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

                        # save figure file
                        directory_name = './figures_output/' + ploidy + '_' + \
                                         bud + '_figures'
                        if not os.path.isdir(directory_name):
                            directory_name = directory_name[2:]
                            os.mkdir(directory_name)
                        file_name = 'figures_output/' + \
                            ploidy + '_' + \
                            bud + \
                            '_figures/' + \
                            end_condition + '_' + \
                            y_axis_file_name[measurement] + \
                            '_' + diff + \
                            '_violin_plot.pdf'
                        fig.write_image(file_name)


main()
