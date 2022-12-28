import pandas as pd
import xarray as xr

def mif_to_xarray(filepath):
    """

    """
    dataframe = pd.read_csv(
        filepath,
        sep=";",
        index_col=["Region", "Variable", "Unit"],
        encoding="latin-1",
    ).drop(columns=["Model", "Scenario"])

    # if new sub-European regions are present, we remove EUR and NEU
    if any(
        x in dataframe.index.get_level_values("Region").unique()
        for x in ["ESC", "DEU", "NEN"]
    ):
        dataframe = dataframe.loc[
            ~dataframe.index.get_level_values("Region").isin(["EUR", "NEU"])
        ]

    if len(dataframe.columns == 20):
        dataframe.drop(columns=dataframe.columns[-1], inplace=True)

    dataframe.columns = dataframe.columns.astype(int)
    dataframe = dataframe.reset_index()

    # dataframe = dataframe.loc[dataframe["Variable"].isin(variables)]

    dataframe = dataframe.rename(
        columns={"Region": "region", "Variable": "variables", "Unit": "unit"}
    )

    array = (
        dataframe.melt(
            id_vars=["region", "variables", "unit"],
            var_name="year",
            value_name="value",
        )[["region", "variables", "year", "value"]]
        .groupby(["region", "variables", "year"])["value"]
        .mean()
        .to_xarray()
    )

    return array