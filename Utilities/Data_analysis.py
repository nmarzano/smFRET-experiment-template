import pandas as pd


def filter_TDP(data_frame, thresh = 0.3):
    filtered_mol = []
    for treatment, df in data_frame.groupby("treatment_name"):
        mol_list = df[(df["FRET before transition"] <= thresh)|(df["FRET after transition"] <= thresh)].Molecule.unique().tolist()
        filtered = df[df["Molecule"].isin(mol_list)]
        filtered_mol.append(filtered)
    filtered_mol = pd.concat(filtered_mol)
    return filtered_mol






def remove_outliers(compiled, plot_type, data_type = "raw"):
    if plot_type == 'hist':
        if data_type == "raw":
            rawFRET = compiled[(compiled[3] > -0.5) & (compiled[3] < 1.5)].copy()
            return rawFRET
        if data_type == "idealized":
            idealizedFRET = compiled[(compiled[4] > -0.5) & (compiled[4] < 1.5)].copy()
            return idealizedFRET
    elif plot_type == 'TDP':
        outliers = compiled[(compiled["FRET before transition"] < -0.5)|(compiled["FRET before transition"] > 1.5)|(compiled["FRET after transition"] < -0.5) | (compiled["FRET after transition"] > 1.5)].index
        compiled.drop(outliers, inplace = True)
        return compiled
    else:
        print('invalid plot type, please set plot_type as "hist" or "TDP" - you idiot')