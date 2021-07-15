import pandas as pd
from pathlib import Path

def main():

    budding_pattern_list = ['not_unipolar']
    ploidy_list = ['haploid']
    diffusion_steps_list = ['with_diffusion']
    conc_list = ['rich', 'low']
    strengths_list = ['strong', 'weak', 'no']
    time_count = 'time'

    output_data_folder = Path("python_output/bud_scar_angles")

    for bud in budding_pattern_list:
        for ploidy in ploidy_list:
            for diffusion in diffusion_steps_list:
                for conc in conc_list:
                    for strength in strengths_list:

                        input_data_folder = Path("matlab_output")

                        # assumes a folder titled "python_output" is present
                        output_data_folder = Path("python_output")

                        file_name_list = [['-1', '-1'], ['-1', '0'],
                                            ['-1', '1'], ['0', '-1'],
                                            ['0', '1'], ['1', '-1'],
                                            ['1', '0'], ['1', '1']]

                        for direction in file_name_list:
                            MF = '(' + direction[0] + '_' + direction[1] + ')'

                            name = time_count + '_' + bud + '_' + ploidy + '_' + conc + '_' + strength + 'MF_' + MF + '_' + diffusion + '_file'

                            file_name1 = 'budScarAngles_' + name + '.xlsx'
                            file_name2 = name + '.xlsx'

                            file1 = input_data_folder / file_name1
                            file2 = input_data_folder / file_name2

                            results1 = pd.read_excel(file1)
                            results2 = pd.read_excel(file2)

                            for run in range(200):
                                cell_count = results2['finalCellCounts'][run]
                                column_name = 'budScarAnglesCell' + str(run + 1)
                                results1[column_name] = results1[column_name][0:cell_count]

                            csv_file_name = 'budScarAngles_' + name + '.csv'

                            results1.to_csv(output_data_folder / csv_file_name, index=False)

main()

