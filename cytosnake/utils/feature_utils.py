def infer_dp_features(dp_profile) -> list[str]:
    """Returns a list of Deep profiler features found within the single cell
    dataframe

    Parameters
    ----------
    dp_profile : pd.DataFrame
        dataframe features captured from deep profiler

    Returns
    -------
    list[str]
        list of deep profiler features

    Raises
    ------
    ValueError
        Raised if no Deep profiler features are found within the given DataFrame
    """
    metadata_model = dp_profile["Metadata_Model"].unique().tolist()
    dp_features = [
        column
        for column in dp_profile.columns.tolist()
        if any(column.startswith(f"{meta_model}_") for meta_model in metadata_model)
    ]
    if len(dp_features) <= 0:
        raise ValueError(
            "No DP features found."
            "Are you sure that this dataframe is from DeepProfiler?"
        )

    return dp_features
