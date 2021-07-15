import researchpy as rp
import pandas as pd
from pathlib import Path
from scipy.stats import levene


def main():
    end_condition = 'time'
    ploidy = 'haploid'
    bud = 'not_unipolar'

    path_name = 'python_output/' + end_condition + \
        '_all_results_' + ploidy + '_' + bud + '.csv'
    file = Path(path_name)
    df = pd.read_csv(file)

    path_name = 'python_output/summary.txt'
    summary_file = Path(path_name)

    sum_file = open(summary_file, 'w')

    # # Compare area of colonies in rich nutrient conditions
    # # exposed to strong axial magnetic fields vs no magnetic fields
    # t_test(df, 'area', ['strong', 'rich'], ['axis', 'none'], sum_file)

    # # Compare area of colonies in low nutrient conditions
    # # exposed to strong axial magnetic fields vs no magnetic fields
    # t_test(df, 'area', ['strong', 'low'], ['axis', 'none'], sum_file)

    # # Compare area of colonies in rich nutrient conditions
    # # exposed to weak axial magnetic fields vs no magnetic fields
    # t_test(df, 'area', ['weak', 'rich'], ['axis', 'none'], sum_file)

    # # Compare area of colonies in low nutrient conditions
    # # exposed to strong diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'area', ['strong', 'low'], ['diagonal', 'none'], sum_file)

    # # Compare area of colonies in low nutrient conditions
    # # exposed to strong diagonal magnetic fields vs weak magnetic fields
    # t_test(df, 'area', ['weak', 'low'], ['diagonal', 'none'], sum_file)

    # # Compare compactness of colonies in low nutrient conditions
    # # exposed to strong diagonal magnetic fields vs weak magnetic fields
    # t_test(df, 'compactness', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # # Compare compactness of colonies in rich nutrient conditions
    # # exposed to strong diagonal magnetic fields vs weak magnetic fields
    # t_test(df, 'compactness', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # Compare compactness of colonies in low nutrient conditions
    # vs in rich nutrient conditions, no magnetic fields
    t_test(df, 'compactness', ['none', 'no'], ['rich', 'low'], sum_file)

    # Compare compactness of colonies in rich nutrient conditions
    # exposed to weak axial magnetic fields vs no magnetic fields
    t_test(df, 'compactness', ['strong', 'rich'],
           ['axis', 'none'], sum_file)

    # Compare compactness of colonies in low nutrient conditions
    # exposed to weak axial magnetic fields vs no magnetic fields
    t_test(df, 'compactness', ['strong', 'low'],
           ['axis', 'none'], sum_file)

    # # Compare compactness of colonies in rich nutrient conditions
    # # exposed to strong diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'compactness', ['strong', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare compactness of colonies in low nutrient conditions
    # # exposed to strong diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'compactness', ['strong', 'low'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare convexity of colonies in low nutrient conditions
    # # vs in rich nutrient conditions, no magnetic fields
    # t_test(df, 'convexity', ['none', 'no'], ['rich', 'low'], sum_file)

    # # Compare convexity of colonies in rich nutrient conditions
    # # exposed to weak axial magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['weak', 'rich'], ['axis', 'none'], sum_file)

    # # Compare convexity of colonies in low nutrient conditions
    # # exposed to weak axial magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['weak', 'low'], ['axis', 'none'], sum_file)

    # # Compare convexity of colonies in rich nutrient conditions
    # # exposed to weak diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['weak', 'rich'], ['diagonal', 'none'], sum_file)

    # # Compare convexity of colonies in low nutrient conditions
    # # exposed to weak diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['weak', 'low'], ['diagonal', 'none'], sum_file)

    # # Compare convexity of colonies in rich nutrient conditions
    # # exposed to strong diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['strong', 'rich'], ['diagonal', 'none'], sum_file)

    # # Compare convexity of colonies in low nutrient conditions
    # # exposed to strong diagonal magnetic fields vs no magnetic fields
    # t_test(df, 'convexity', ['strong', 'low'], ['diagonal', 'none'], sum_file)

    # # Compare elongation of colonies in rich colonies
    # # exposed to strong vs weak axial magnetic fields
    # t_test(df, 'elongation', ['axis', 'rich'], ['strong', 'weak'],
    #        sum_file)

    # # Compare elongation of colonies in low nutrient conditions
    # # exposed to strong vs weak axial magnetic fields
    # t_test(df, 'elongation', ['axis', 'low'], ['strong', 'weak'], sum_file)

    # # Compare elongation of colonies in low nutrient conditions
    # # vs in rich nutrient conditions, no magnetic fields
    # t_test(df, 'elongation', ['none', 'no'], ['rich', 'low'], sum_file)

    # # Compare boundary fluctuations of colonies in low nutrient conditions
    # # vs rich nutrient conditions, no magnetic fields
    # t_test(df, 'boundary fluctuations', ['none', 'no'], ['rich', 'low'],
    #        sum_file)

    # # Compare boundary fluctuations of colonies in rich nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'boundary fluctuations', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # # Compare boundary fluctuations of colonies in low nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'boundary fluctuations', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # # Compare boundary fluctuations of colonies in rich nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'boundary fluctuations', ['strong', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare boundary fluctuations of colonies in rich nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'boundary fluctuations', ['weak', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare boundary fluctuations of colonies in low nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'boundary fluctuations', ['strong', 'low'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare perimeter of colonies in rich nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'perimeter', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # # Compare perimeter of colonies in low nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'perimeter', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # # Compare perimeter of colonies in rich nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'perimeter', ['strong', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare perimeter of colonies in low nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'perimeter', ['strong', 'low'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare roundness of colonies in rich nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'roundness', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # # Compare roundness of colonies in low nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'roundness', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # Compare roundness of colonies in rich nutrient conditions
    # exposed to strong vs weak axial magnetic fields
    t_test(df, 'roundness', ['axis', 'rich'], ['strong', 'weak'], sum_file)

    # Compare roundness of colonies in low nutrient conditions
    # exposed to strong vs weak axial magnetic fields
    t_test(df, 'roundness', ['axis', 'low'], ['strong', 'weak'], sum_file)

    # # Compare roundness of colonies in rich nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'roundness', ['strong', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare roundness of colonies in low nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'roundness', ['strong', 'low'],
    #        ['diagonal', 'none'], sum_file)

    # Compare roundness of colonies in rich nutrient conditions
    # exposed to strong axial fields vs no magnetic fields
    t_test(df, 'roundness', ['strong', 'rich'], ['axis', 'none'], sum_file)

    # Compare roundness of colonies in low nutrient conditions
    # exposed to strong axial fields vs no magnetic fields
    t_test(df, 'roundness', ['strong', 'low'], ['axis', 'none'], sum_file)

    # Compare roundness of colonies in rich nutrient conditions
    # exposed to weak axial fields vs no magnetic fields
    t_test(df, 'roundness', ['weak', 'rich'], ['axis', 'none'], sum_file)

    # Compare roundness of colonies in low nutrient conditions
    # exposed to weak axial fields vs no magnetic fields
    t_test(df, 'roundness', ['weak', 'low'], ['axis', 'none'], sum_file)

    # # Compare solidity of colonies in rich nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'solidity', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # # Compare solidity of colonies in low nutrient conditions
    # # exposed to strong vs weak diagonal magnetic fields
    # t_test(df, 'solidity', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # # Compare solidity of colonies in rich nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'solidity', ['strong', 'rich'],
    #        ['diagonal', 'none'], sum_file)

    # # Compare solidity of colonies in low nutrient conditions
    # # exposed to strong diagonal fields vs no magnetic fields
    # t_test(df, 'solidity', ['strong', 'low'],
    #        ['diagonal', 'none'], sum_file)

    # # Switch to count data to test time

    # end_condition = 'count'

    # path_name = 'python_output/' + end_condition + \
    #     '_all_results_' + ploidy + '_' + bud + '.csv'
    # file = Path(path_name)
    # df = pd.read_csv(file)

    # # Compare time to reach 10 000 cells of colonies in
    # # rich nutrient conditions exposed to strong vs weak
    # # diagonal magnetic fields
    # t_test(df, 'time', ['diagonal', 'rich'],
    #        ['strong', 'weak'], sum_file)

    # # Compare time to reach 10 000 cells of colonies in
    # # low nutrient conditions exposed to strong vs weak
    # # diagonal magnetic fields
    # t_test(df, 'time', ['diagonal', 'low'],
    #        ['strong', 'weak'], sum_file)

    # Close file
    sum_file.close()


def t_test(df, measurement, constants, variable, sum_file):
    '''
    Runs two-sided t-test, prints detailed summary and creates summary file of
    all p-values.

    Measurement string can be one of:
        cellno,
        roundness,
        compactness,
        convexity,
        boundary fluctuations,
        time,
        perimeter,
        solidity,
        elongation,
        area
    Strings in the constants and variable lists can be a:
    Direction - axis, diagonal, or none
    MF Strength - strong, weak, or no
    Nutrient concentration - rich or low

    If the comparison is between magnetic fields and no magnetic fields,
    the strength of the magnetic field applied should be in constants

    Input: df (dataframe) - dataframe containing all results
           measurement (str) - name of measurement being compared
           constants (list) - strings of the conditions kept constant
           variable (list) - strings of the two conditions being compared,
                             should be of the same type
           sum_file (file) - summary file to add p-values to
    Output: N/A
    '''
    measure_names = dict([
        ('cellno', 'finalCellCounts'),
        ('roundness', 'roundnessList'),
        ('compactness', 'compactnessList'),
        ('convexity', 'convexityList'),
        ('boundary fluctuations', 'boundary_fluctuations'),
        ('time', 'timeList'),
        ('perimeter', 'perimeterList'),
        ('solidity', 'solidityList'),
        ('elongation', 'elongationList'),
        ('area', 'areaList')
    ])
    variable_types = dict([
        ('axis', 'directionType'),
        ('diagonal', 'directionType'),
        ('none', 'directionType'),
        ('strong', 'strength'),
        ('weak', 'strength'),
        ('no', 'strength'),
        ('rich', 'concentrations'),
        ('low', 'concentrations')
    ])

    title = build_title(measurement, constants, variable, variable_types)

    if variable[0] == 'none':
        for constant in constants:
            if variable_types[constant] == 'strength':
                strength_ind = constants.index(constant)
        old_constant = constants[strength_ind]
        constants[strength_ind] = 'no'

    group1 = df[measure_names[measurement]][(
        df[variable_types[constants[0]]] == constants[0]) & (
        df[variable_types[constants[1]]] == constants[1]) & (
        df[variable_types[variable[0]]] == variable[0])]

    if variable[1] == 'none':
        for constant in constants:
            if variable_types[constant] == 'strength':
                strength_ind = constants.index(constant)
        old_constant = constants[strength_ind]
        constants[strength_ind] = 'no'

    group2 = df[measure_names[measurement]][(
        df[variable_types[constants[0]]] == constants[0]) & (
        df[variable_types[constants[1]]] == constants[1]) & (
        df[variable_types[variable[1]]] == variable[1])]

    if variable[0] == 'none':
        constants[strength_ind] = old_constant

    print(title)
    print("Levene's test")
    stat, p = levene(group1, group2)
    print(p)
    if p < 0.05:
        equalVar = False
    else:
        equalVar = True
    summary, results = rp.ttest(group1=group1,
                                group2=group2,
                                group1_name=variable[0],
                                group2_name=variable[1],
                                equal_variances=equalVar)
    print(summary)
    print(results)

    if results['results'][3] < 0.05:
        significance = 'significant'
    else:
        significance = 'insignificant'
    print(significance)

    sum_file.write('{} \n'.format(title))
    sum_file.write('p = {} \n'.format(results['results'][3]))
    sum_file.write('{} \n \n'.format(significance))


def build_title(measurement, constants, variable, variable_types):
    if variable_types[variable[0]] == 'concentrations':
        variable_str = '{} vs {} Nutrients'.format(variable[0], variable[1])
        constant_str = 'NO MF'
    elif variable_types[variable[0]] == 'strength':
        for constant in constants:
            if variable_types[constant] == 'directionType':
                if constant == 'axis':
                    constant = 'axial'
                direction = constant
            elif variable_types[constant] == 'concentrations':
                concentration = constant
        variable_str = '{} vs {} {} MF'.format(
            variable[0], variable[1], direction)
        constant_str = '{} Nutrients'.format(concentration)
    elif 'none' in variable:
        index = variable.index('none')
        second_ind = (len(variable) - index - 1) % len(variable)
        if variable[second_ind] == 'axis':
            direction = 'axial'
        else:
            direction = variable[second_ind]
        for constant in constants:
            if variable_types[constant] == 'strength':
                strength = constant
            elif variable_types[constant] == 'concentrations':
                concentration = constant
        variable_str = '{} {} vs No MF'.format(strength, direction)
        constant_str = '{} Nutrients'.format(concentration)
    else:
        for constant in constants:
            if variable_types[constant] == 'strength':
                strength = constant
            elif variable_types[constant] == 'concentrations':
                concentration = constant
        variable_str = '{} {} vs {} MF'.format(
            strength, variable[0], variable[1])
        constant_str = '{} Nutrients'.format(concentration)

    title = ('{} - {}, {}'.format(measurement,
                                  variable_str,
                                  constant_str)).upper()

    return title


main()
