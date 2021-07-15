import pandas as pd
from pathlib import Path
import plotly.graph_objects as go
import os

## Creates a polar bar chart of the budding angles

def main():

    ploidy_list = ['haploid','diploid']
    conc_list = ['rich', 'low']
    strengths_list = ['strong', 'weak', 'no']
    time_count = 'time'
    bud = 'not_unipolar'
    diffusion = 'with_diffusion'

    for ploidy in ploidy_list:
        for conc in conc_list:
            for strength in strengths_list:

                input_data_folder = Path("matlab_output")

                # assumes a folder titled "python_output" is present
                output_data_folder = Path("python_output")

                file_name_list = [['-1', '-1'], ['-1', '0'],
                                    ['-1', '1'], ['0', '-1'],
                                    ['0', '1'], ['1', '-1'],
                                    ['1', '0'], ['1', '1']]

                # assign magnetic field type to MF directions
                MF_dictionary = dict([('(-1, -1)', 'diagonal'),
                                        ('(-1, 0)', 'axis'),
                                        ('(-1, 1)', 'diagonal'),
                                        ('(0, -1)', 'axis'),
                                        ('(0, 1)', 'axis'),
                                        ('(1, -1)', 'diagonal'),
                                        ('(1, 0)', 'axis'),
                                        ('(1, 1)', 'diagonal')])

                results_list = []

                # combine data from all MF directions for the given
                # conditions:
                # ploidy, strength, concentration
                for MF_list in file_name_list:

                    x = MF_list[0]

                    y = MF_list[1]

                    MF = x+'_'+y

                    file_name = 'buddingAngles_' + time_count + '_' + bud + '_' + \
                        ploidy+'_'+conc+'_' + \
                        strength+'MF_('+MF+')_'+diffusion+'_file.xlsx'

                    folder1 = Path('matlab_output')
                    folder2 = Path('budding_angles')

                    file = folder1 / folder2 / file_name

                    results = pd.read_excel(file)

                    results = results[results.buddingAngles != -1]

                    zero = (results.buddingAngles.values == 0).sum()
                    fortyfive = (results.buddingAngles.values == 45).sum()
                    ninety = (results.buddingAngles.values == 90).sum()
                    onethirtyfive = (results.buddingAngles.values == 135).sum()
                    oneeighty = (results.buddingAngles.values == 180).sum()
                    negfortyfive = (results.buddingAngles.values == -45).sum()
                    negninety = (results.buddingAngles.values == -90).sum()
                    negonethirtyfive = (results.buddingAngles.values == -135).sum()

                    fig = go.Figure(go.Barpolar(
                        r=[zero, fortyfive, ninety, onethirtyfive,oneeighty,negfortyfive,negninety,negonethirtyfive],
                        theta=[0,45,90,135,180,-135,-90,-45],
                        width=[1,1,1,1,1,1,1,1,],
                        marker_color=["#E4FF87", '#709BFF', '#709BFF', '#FFAA70', '#FFAA70', '#FFDF70', '#B6FFB4'],
                        marker_line_color="black",
                        marker_line_width=2,
                        opacity=0.8
                    ))

                    fig.update_layout(
                        template=None,
                        polar = dict(
                            sector=list([-180,180]),
                            radialaxis = dict(showticklabels=True, ticks='outside',nticks=4),
                            angularaxis = dict(showticklabels=True, ticks='outside')
                        )
                    )

                    # save figure file
                    directory_name = './figures_output/budding_angle_figures'
                    if not os.path.isdir(directory_name):
                        directory_name = directory_name[2:]
                        os.mkdir(directory_name)
                    file_name = 'figures_output/budding_angle_figures/' + \
                        time_count + '_' + bud + '_' + \
                        ploidy+'_'+conc+'_' + \
                        strength+'MF_('+MF+')_'+diffusion+'_angles_plot.pdf'
                    fig.write_image(file_name)
main()
