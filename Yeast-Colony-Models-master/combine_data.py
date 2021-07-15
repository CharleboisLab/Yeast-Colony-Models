import pandas as pd
from pathlib import Path

# Collects files from format_data.py into two large CSV files for each ploidy


def main():

    input_data_folder = Path("python_output")

    output_data_folder = Path("python_output")

    ploidy_list = ['diploid', 'haploid']

    concentration_list = ['rich', 'low']

    strength_list = ['extraStrong', 'strong', 'weak', 'no']

    diffusion_list = ['with_diffusion']

    budding_pattern_list = ['not_unipolar']

    time_count = 'time'

    for ploidy in ploidy_list:
        for bud in budding_pattern_list:
            results_list = []
            for diff in diffusion_list:
                for conc in concentration_list:
                    for stre in strength_list:

                        file_name = time_count + '_' + bud + '_' + ploidy + \
                            '_' + conc + '_nutrient_' + diff + \
                            '_' + stre + '_MF.csv'

                        file = input_data_folder / file_name

                        results = pd.read_csv(file)

                        results_list.append(results)

            complete_table = results_list[0]
            for table in results_list[1:]:
                complete_table = pd.concat([complete_table, table])

            complete_table.loc[(complete_table.strength == 'no'),
                               'directionType'] = 'none'

            output_file = time_count+'_all_results_'+ploidy+'_'+bud+'.csv'
            complete_table.to_csv(output_data_folder/output_file, index=False)


main()
