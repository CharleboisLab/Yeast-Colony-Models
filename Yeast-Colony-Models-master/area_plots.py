import pandas as pd
import plotly.graph_objects as go
import os
from pathlib import Path

# Creates violin plots for area measurements, one for each nutrient condition


def main():

    # set font sizes
    axis_font_size = 32
    xaxis_tick_font_size = 30
    yaxis_tick_font_size = 24
    legend_font_size = 30

    # set image dimensions
    height_int = 1250
    width_int = 1700

    # add 'diploid' to ploidy_list to make plots for diploid measurments
    ploidy_list = ['haploid']
    concentration_list = ['rich', 'low']
    diffusion_list = ['no_diffusion', 'with_diffusion']
    budding_list = ['not_unipolar']

    x_axis_labels = dict([
        ('no_diffusion', 'no diffusion'),
        ('with_diffusion', 'diffusion')
    ])

    for ploidy in ploidy_list:
        for bud in budding_list:
            path_name = "python_output/all_results_"+ploidy+'_'+bud+".csv"
            print(path_name)
            file = Path(path_name)

            for diff in diffusion_list:
                df = pd.read_csv(file)
                df = df[df.diffusion == diff]

                for concentration in concentration_list:

                    # make x-axis titles
                    x_axis_title = 'magnetic field direction applied to ' + \
                        'colonies in ' + concentration + ' nutrient conditions' + \
                        ' with ' + x_axis_labels[diff]

                    fig = go.Figure()

                    fig.add_trace(go.Violin(x=df['directionType'][
                        df['strength'] == 'weak'],
                        y=df['areaList'][(
                            df['strength'] == 'weak') & (
                            df['concentrations'] == concentration)],
                        legendgroup='Field Strength',
                        scalegroup='Field Strength',
                        name='Weak MFs',
                        side='negative',
                        line_color='blue')
                    )
                    fig.add_trace(go.Violin(x=df['directionType'][
                        df['strength'] == 'strong'],
                        y=df['areaList'][(
                            df['strength'] == 'strong') & (
                            df['concentrations'] == concentration)],
                        legendgroup='Field Strength',
                        scalegroup='Field Strength',
                        name='Strong MFs',
                        side='positive',
                        line_color='orange')
                    )
                    fig.add_trace(go.Violin(x=df['directionType'][
                        df['strength'] == 'no'],
                        y=df['areaList'][(
                            df['strength'] == 'no') & (
                            df['concentrations'] == concentration)],
                        legendgroup='Field Strength',
                        scalegroup='Field Strength',
                        name='No MFs',
                        line_color='purple')
                    )

                    # formatting the figure
                    fig.update_traces(meanline_visible=False, points=False, box_visible=True)
                    fig.update_layout(yaxis_title='area (pixels)',
                                      xaxis_title=x_axis_title,
                                      violinmode='overlay',
                                      paper_bgcolor='rgba(255,255,255,255)',
                                      plot_bgcolor='rgba(255,255,255,255)',
                                      legend=dict(
                                          font=dict(size=legend_font_size)),
                                      height=height_int,
                                      width=width_int)
                    fig.update_xaxes(title_font=dict(size=axis_font_size),
                                     tickfont=dict(size=xaxis_tick_font_size))
                    fig.update_yaxes(title_font=dict(size=axis_font_size),
                                     tickfont=dict(size=yaxis_tick_font_size),
                                     gridcolor='Grey')

                    # save figure file
                    directory_name = './figures_output/' + ploidy + '_' + bud + '_figures'
                    if not os.path.isdir(directory_name):
                        directory_name = directory_name[2:]
                        os.mkdir(directory_name)
                    file_name = 'figures_output/' + ploidy + '_' + bud + '_figures/' + \
                        'area_' + concentration + '_' + diff + '_violin_plot.png'
                    fig.write_image(file_name)
                    # fig.show()


main()
