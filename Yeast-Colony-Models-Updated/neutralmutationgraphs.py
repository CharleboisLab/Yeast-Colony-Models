import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():

    # fileName = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults.csv')
    # results = pd.read_csv(fileName)

    # fileName1 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults1.csv')
    # results1 = pd.read_csv(fileName1)

    # fileName2 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults2.csv')
    # results2 = pd.read_csv(fileName2)

    # fileName3 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults3.csv')
    # results3 = pd.read_csv(fileName3)

    # fileName4 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults4.csv')
    # results4 = pd.read_csv(fileName4)

    # fileName5 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults5.csv')
    # results5 = pd.read_csv(fileName5)

    # fileName6 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults6.csv')
    # results6 = pd.read_csv(fileName6)

    # fileName7 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults7.csv')
    # results7 = pd.read_csv(fileName7)

    # fileName8 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults8.csv')
    # results8 = pd.read_csv(fileName8)

    # fileName9 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults9.csv')
    # results9 = pd.read_csv(fileName9)

    # fileName10 = Path('0.0008208_50wtfit_neutralmutationresults/neutralmutationresults10.csv')
    # results10 = pd.read_csv(fileName10)

    cells = np.zeros([2000,90])
    mutants = np.zeros([2000,90])
    subcolonies = np.zeros([2000,90])
    subcolonyareas = np.zeros([2000,90])

    n = 7   # plot height
    m = 10   # plot width

    for i in range(1,2001):
        no = str(i)
        name = '0.0008208_75wtfit_neutralmutationresults1000/neutralmutationresults' + no + '.csv'
        fileName = Path(name)
        results = pd.read_csv(fileName)

        cell = results.cellCount[:90]
        mutant = results.mutantCount[:90]
        subcolony = results.subColonyCount[:90]
        subcolonyarea = mutant/cell

        cells[i-1,:] = cell
        mutants[i-1,:] = mutant
        subcolonies[i-1,:] = subcolony
        subcolonyareas[i-1,:] = subcolonyarea

    plt.figure(figsize=(n,m))
    plt.plot(cells.T)

    plt.xlabel('days')
    plt.ylabel('total cell count')
    plt.figure(1)

    plt.figure(figsize=(n,m))
    plt.plot(mutants.T)

    plt.xlabel('days')
    plt.ylabel('mutant cell count')
    plt.figure(2)

    plt.figure(figsize=(n,m))
    plt.plot(subcolonies.T)

    plt.xlabel('days')
    plt.ylabel('mutant subcolony count')
    plt.figure(3)

    plt.figure(figsize=(n,m))
    plt.plot(subcolonyareas.T)

    plt.xlabel('days')
    plt.ylabel('mutant subcolony area fraction')
    plt.figure(4)

    # finalCounts = np.zeros(100)

    # for i in range(1,101):
    #     no = str(i)
    #     name = 'neutralmutationresults/neutralmutationresults' + no + '.csv'
    #     fileName = Path(name)
    #     results = pd.read_csv(fileName)
    #     finalCount = results.mutantCount[90]
    #     finalCounts[i-1] = finalCount

    # plt.hist(finalCounts,bins=30)
    # plt.xlabel('mutant cell counts')
    # plt.show()

    # n = 7   # plot height
    # m = 10   # plot width

    # plt.figure(1)
    # plt.figure(figsize=(n,m))

    # cells1 = results1.cellCount[:90]
    # cells2 = results2.cellCount[:90]
    # cells3 = results3.cellCount[:90]
    # cells4 = results4.cellCount[:90]
    # cells5 = results5.cellCount[:90]
    # cells6 = results6.cellCount[:90]
    # cells7 = results7.cellCount[:90]
    # cells8 = results8.cellCount[:90]
    # cells9 = results9.cellCount[:90]
    # cells10 = results10.cellCount[:90]

    # mutants1 = results1.mutantCount[:90]
    # mutants2 = results2.mutantCount[:90]
    # mutants3 = results3.mutantCount[:90]
    # mutants4 = results4.mutantCount[:90]
    # mutants5 = results5.mutantCount[:90]
    # mutants6 = results6.mutantCount[:90]
    # mutants7 = results7.mutantCount[:90]
    # mutants8 = results8.mutantCount[:90]
    # mutants9 = results9.mutantCount[:90]
    # mutants10 = results10.mutantCount[:90]

    # plt.plot(mutants1)
    # plt.plot(mutants2)
    # plt.plot(mutants3)
    # plt.plot(mutants4)
    # plt.plot(mutants5)
    # plt.plot(mutants6)
    # plt.plot(mutants7)
    # plt.plot(mutants8)
    # plt.plot(mutants9)
    # plt.plot(mutants10)

    # plt.xlabel('days')
    # plt.ylabel('mutant cell count')
    # plt.legend(['simulation run 1','simulation run 2','simulation run 3','simulation run 4','simulation run 5','simulation run 6','simulation run 7','simulation run 8','simulation run 9','simulation run 10'])
    # # plt.legend(['simulation run 1','simulation run 4','simulation run 5','simulation run 6','simulation run 9','simulation run 11','simulation run 12','simulation run 13','simulation run 14','simulation run 15'])

    # plt.figure(2)
    # plt.figure(figsize=(n,m))

    # subcolonies1 = results1.subColonyCount[:90]
    # subcolonies2 = results2.subColonyCount[:90]
    # subcolonies3 = results3.subColonyCount[:90]
    # subcolonies4 = results4.subColonyCount[:90]
    # subcolonies5 = results5.subColonyCount[:90]
    # subcolonies6 = results6.subColonyCount[:90]
    # subcolonies7 = results7.subColonyCount[:90]
    # subcolonies8 = results8.subColonyCount[:90]
    # subcolonies9 = results9.subColonyCount[:90]
    # subcolonies10 = results10.subColonyCount[:90]

    # plt.plot(subcolonies1)
    # plt.plot(subcolonies2)
    # plt.plot(subcolonies3)
    # plt.plot(subcolonies4)
    # plt.plot(subcolonies5)
    # plt.plot(subcolonies6)
    # plt.plot(subcolonies7)
    # plt.plot(subcolonies8)
    # plt.plot(subcolonies9)
    # plt.plot(subcolonies10)

    # plt.xlabel('days')
    # plt.ylabel('mutant subcolony count')
    # plt.legend(['simulation run 1','simulation run 2','simulation run 3','simulation run 4','simulation run 5','simulation run 6','simulation run 7','simulation run 8','simulation run 9','simulation run 10'])
    # # plt.legend(['simulation run 1','simulation run 4','simulation run 5','simulation run 6','simulation run 9','simulation run 11','simulation run 12','simulation run 13','simulation run 14','simulation run 15'])

    # plt.figure(3)
    # plt.figure(figsize=(n,m))

    # subcolonyarea1 = mutants1/cells1
    # subcolonyarea2 = mutants2/cells2
    # subcolonyarea3 = mutants3/cells3
    # subcolonyarea4 = mutants4/cells4
    # subcolonyarea5 = mutants5/cells5
    # subcolonyarea6 = mutants6/cells6
    # subcolonyarea7 = mutants7/cells7
    # subcolonyarea8 = mutants8/cells8
    # subcolonyarea9 = mutants9/cells9
    # subcolonyarea10 = mutants10/cells10

    # plt.plot(subcolonyarea1)
    # plt.plot(subcolonyarea2)
    # plt.plot(subcolonyarea3)
    # plt.plot(subcolonyarea4)
    # plt.plot(subcolonyarea5)
    # plt.plot(subcolonyarea6)
    # plt.plot(subcolonyarea7)
    # plt.plot(subcolonyarea8)
    # plt.plot(subcolonyarea9)
    # plt.plot(subcolonyarea10)

    # plt.xlabel('days')
    # plt.ylabel('mutant subcolony area fraction')
    # plt.legend(['simulation run 1','simulation run 2','simulation run 3','simulation run 4','simulation run 5','simulation run 6','simulation run 7','simulation run 8','simulation run 9','simulation run 10'])
    
    plt.show()

main()