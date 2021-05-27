import pandas as pd
from pathlib import Path

# Collects all MATLAB output files from matlab_output into a collection of
# files for each magnetic field direction and ploidy
# Requires "python_output" folder in directory


def main():

    budding_pattern_list = ['not_unipolar']
    ploidy_list = ['haploid', 'diploid']
    diffusion_steps_list = ['with_diffusion']
    conc_list = ['rich', 'low']
    strengths_list = ['strong', 'weak', 'no']
    time_count = 'count'

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
                        # ploidy, diffusion, strength, concentration
                        for MF_list in file_name_list:

                            x = MF_list[0]

                            y = MF_list[1]

                            MF = x+'_'+y

                            file_name = time_count + '_' + bud + '_' + \
                                ploidy+'_'+conc+'_' + \
                                strength+'MF_('+MF+')_'+diffusion+'_file.xlsx'

                            file = input_data_folder / file_name

                            results = pd.read_excel(file)

                            direction_string = '('+x+', '+y+')'

                            direction_list = [direction_string] * len(results)

                            concentration_list = [conc] * len(results)

                            strength_list = [strength] * len(results)

                            diffusion_list = [diffusion] * len(results)

                            results = results.drop(columns='densityList')

                            results['direction'] = direction_list

                            results['directionType'] = results[
                                    'direction'].map(MF_dictionary)

                            results['concentrations'] = concentration_list

                            results['strength'] = strength_list

                            results['diffusion'] = diffusion_list

                            results['boundary_fluctuations'] = \
                                results['stdDevLengthList'] \
                                / results['meanLengthsList']

                            results_list.append(results)

                        complete_table = results_list[0]
                        for table in results_list[1:]:
                            complete_table = pd.concat([complete_table, table])

                        output = time_count + '_' + bud + '_' + ploidy + \
                            '_' + conc + '_nutrient_' + diffusion + \
                            '_' + strength + '_MF.csv'

                        complete_table.to_csv(output_data_folder /
                                              output, index=False)


main()
